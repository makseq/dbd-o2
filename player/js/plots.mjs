// plots.mjs — графики плеера 2D-ДБР (canvas 2D, тёмная тема, ретина).
//
// Классы:
//   WaveformPlot       U_app/U_gap (кВ, левая ось) + ток (мА, правая, лог/лин)
//   LissajousPlot      Q–V параллелограмм, автофит наклонов -> C_diel, C_cell, W, P
//   RadialProfilePlot  sigma(r) на обеих поверхностях диэлектрика
//   AxialProfilePlot   лог-профили плотностей вдоль оси, легенда-переключатель
//   PhotoPanel         три фотоканала против ударной ионизации, лог-ось
// Утилиты:
//   cellVolumes, integrateField, axialColumn, PhotoSeriesBuilder, SPECIES, CHANNELS
//
// Контракт и примеры — PLOTS_API.md (рядом).
//
// ПРАВИЛА, КОТОРЫЕ СОБЛЮДАЮТСЯ ЖЁСТКО:
//  1. Ноль аллокаций в кадре. Все буферы (декимация, тики, подписи) выделяются
//     в _onResize/при смене диапазона, render() только пишет в них.
//  2. Лог-оси не врут: значения <= 0 не рисуются (разрыв пути), их количество
//     показывается в подписи оси. Ток в лог-режиме — это |I| плюс отдельная
//     полоса-индикатор знака, а не молча выброшенная половина данных.
//  3. Ни одна извлечённая величина не показывается без сравнения с аналитикой
//     там, где аналитика есть (ёмкости), и без бейджа качества, если фит плохой.
//
// Зависимости: только ./metrics.mjs. Модуль не знает про index.html/app.mjs.

import {
  THEME, fmtSI, fmtExp, fmtFixed, fmtPct, fmtTime, fmtPow10,
  analyticCapacitances,
} from './metrics.mjs';

const RAF = typeof requestAnimationFrame === 'function'
  ? requestAnimationFrame
  : (fn) => setTimeout(fn, 16);

// ─────────────────────────────────────────────────────────── Scale

/** Линейная или логарифмическая шкала «значение -> пиксель». */
class Scale {
  constructor() {
    this.min = 0; this.max = 1; this.log = false;
    this.p0 = 0; this.p1 = 1; this.k = 1; this.lmin = 0; this.lspan = 1;
  }
  setDomain(min, max, log) {
    this.log = !!log;
    if (this.log) {
      if (!(min > 0)) min = (max > 0 ? max : 1) * 1e-12;
      if (!(max > min)) max = min * 10;
    } else if (!(max > min)) {
      const c = Number.isFinite(min) ? min : 0;
      min = c - 1; max = c + 1;
    }
    this.min = min; this.max = max;
    this._upd();
    return this;
  }
  setPixels(p0, p1) { this.p0 = p0; this.p1 = p1; this._upd(); return this; }
  _upd() {
    if (this.log) {
      this.lmin = Math.log10(this.min);
      this.lspan = Math.log10(this.max) - this.lmin;
      if (!(this.lspan > 0)) this.lspan = 1;
      this.k = (this.p1 - this.p0) / this.lspan;
    } else {
      this.k = (this.p1 - this.p0) / (this.max - this.min);
    }
  }
  /** значение -> пиксель. Для лог-шкалы v<=0 даёт NaN (путь обязан рваться). */
  to(v) {
    if (this.log) {
      if (!(v > 0)) return NaN;
      return this.p0 + (Math.log10(v) - this.lmin) * this.k;
    }
    return this.p0 + (v - this.min) * this.k;
  }
  /** пиксель -> значение. */
  from(px) {
    if (this.log) return Math.pow(10, this.lmin + (px - this.p0) / this.k);
    return this.min + (px - this.p0) / this.k;
  }
}

// ─────────────────────────────────────────────────────────── Ticks

/**
 * Тики с кэшем подписей: строки пересобираются ТОЛЬКО когда изменился
 * диапазон/шаг. В установившемся окне render() не аллоцирует ничего.
 */
class Ticks {
  constructor(cap = 32) {
    this.v = new Float64Array(cap);
    this.lab = new Array(cap).fill('');
    this.minor = new Float64Array(cap * 8);
    this.n = 0; this.nMinor = 0;
    this._a = NaN; this._b = NaN; this._t = 0; this._log = null;
  }
  _same(a, b, target, log) {
    return this._a === a && this._b === b && this._t === target && this._log === log;
  }
  _stamp(a, b, target, log) { this._a = a; this._b = b; this._t = target; this._log = log; }

  linear(min, max, target, fmt) {
    if (this._same(min, max, target, false)) return this;
    this._stamp(min, max, target, false);
    this.n = 0; this.nMinor = 0;
    const span = max - min;
    if (!(span > 0) || !Number.isFinite(span)) return this;
    const raw = span / Math.max(2, target);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const nn = raw / mag;
    const step = mag * (nn < 1.5 ? 1 : nn < 3 ? 2 : nn < 7 ? 5 : 10);
    const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    let x = Math.ceil(min / step - 1e-9) * step;
    while (x <= max + step * 1e-9 && this.n < this.v.length) {
      const val = Math.abs(x) < step * 1e-6 ? 0 : x;
      this.v[this.n] = val;
      this.lab[this.n] = fmt ? fmt(val, dec) : fmtFixed(val, dec);
      this.n++;
      x += step;
    }
    return this;
  }

  log(min, max, target) {
    if (this._same(min, max, target, true)) return this;
    this._stamp(min, max, target, true);
    this.n = 0; this.nMinor = 0;
    if (!(min > 0) || !(max > min)) return this;
    const e0 = Math.floor(Math.log10(min) + 1e-9);
    const e1 = Math.ceil(Math.log10(max) - 1e-9);
    const dec = e1 - e0;
    const stride = Math.max(1, Math.ceil(dec / Math.max(2, target)));
    for (let e = e0; e <= e1 && this.n < this.v.length; e += stride) {
      const val = Math.pow(10, e);
      if (val < min * 0.999 || val > max * 1.001) continue;
      this.v[this.n] = val;
      this.lab[this.n] = fmtPow10(e);
      this.n++;
    }
    if (dec <= 4) {
      for (let e = e0; e <= e1; e++) {
        for (let m = 2; m <= 9; m++) {
          const val = m * Math.pow(10, e);
          if (val < min || val > max || this.nMinor >= this.minor.length) continue;
          this.minor[this.nMinor++] = val;
        }
      }
    }
    return this;
  }
}

// ───────────────────────────────────────────────────────── BasePlot

/**
 * Общий каркас: ретина-размер, ResizeObserver, поля/оси/сетка, курсор,
 * отложенная перерисовка. Наследник реализует draw().
 */
class BasePlot {
  constructor(canvas, opts = {}) {
    if (!canvas) throw new Error('plots: не передан canvas');
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.opts = opts;
    this.pad = Object.assign({ l: 56, r: 54, t: 14, b: 26 }, opts.pad);
    this.title = opts.title || '';
    this.w = 0; this.h = 0; this.dpr = 1;
    this.x = new Scale(); this.y = new Scale();
    this.xt = new Ticks(); this.yt = new Ticks();
    this.cursorT = NaN;
    this._rafPending = false;
    this._boundFrame = () => { this._rafPending = false; this.render(); };
    this._boundResize = () => { this._measure(); this.onResize(); this.requestRender(); };
    this._ro = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(this._boundResize) : null;
    if (this._ro) this._ro.observe(canvas);
    this._measure();
  }

  _measure() {
    const c = this.canvas;
    const dpr = (typeof devicePixelRatio === 'number' && devicePixelRatio > 0) ? devicePixelRatio : 1;
    const w = c.clientWidth || c.width || 320;
    const h = c.clientHeight || c.height || 180;
    const W = Math.max(1, Math.round(w * dpr));
    const H = Math.max(1, Math.round(h * dpr));
    const changed = (this.w !== w || this.h !== h || this.dpr !== dpr);
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    this.w = w; this.h = h; this.dpr = dpr;
    this.px0 = this.pad.l; this.px1 = Math.max(this.pad.l + 1, w - this.pad.r);
    this.py1 = this.pad.t; this.py0 = Math.max(this.pad.t + 1, h - this.pad.b);
    this.x.setPixels(this.px0, this.px1);
    this.y.setPixels(this.py0, this.py1);
    return changed;
  }

  /** Переопределяется наследником: пересоздать пиксельные буферы. */
  onResize() {}

  requestRender() {
    if (this._rafPending) return;
    this._rafPending = true;
    RAF(this._boundFrame);
  }

  /** Синхронная полная перерисовка. */
  render() {
    this._measure();
    if (this.w < 8 || this.h < 8) return;
    const g = this.ctx;
    g.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    g.fillStyle = THEME.surface1;
    g.fillRect(0, 0, this.w, this.h);
    this.draw(g);
  }

  draw() {}

  destroy() {
    if (this._ro) { this._ro.disconnect(); this._ro = null; }
  }

  // ── примитивы

  _clipPlot(g) {
    g.save();
    g.beginPath();
    g.rect(this.px0, this.py1, this.px1 - this.px0, this.py0 - this.py1);
    g.clip();
  }

  _grid(g, sx, sy, xt, yt) {
    g.lineWidth = 1;
    g.strokeStyle = THEME.grid;
    g.beginPath();
    if (yt) {
      for (let i = 0; i < yt.nMinor; i++) {
        const py = Math.round(sy.to(yt.minor[i])) + 0.5;
        if (py < this.py1 || py > this.py0) continue;
        g.moveTo(this.px0, py); g.lineTo(this.px1, py);
      }
      for (let i = 0; i < yt.n; i++) {
        const py = Math.round(sy.to(yt.v[i])) + 0.5;
        if (py < this.py1 || py > this.py0) continue;
        g.moveTo(this.px0, py); g.lineTo(this.px1, py);
      }
    }
    if (xt) {
      for (let i = 0; i < xt.n; i++) {
        const px = Math.round(sx.to(xt.v[i])) + 0.5;
        if (px < this.px0 || px > this.px1) continue;
        g.moveTo(px, this.py1); g.lineTo(px, this.py0);
      }
    }
    g.stroke();
    g.strokeStyle = THEME.axis;
    g.beginPath();
    g.rect(this.px0 + 0.5, this.py1 + 0.5, this.px1 - this.px0 - 1, this.py0 - this.py1 - 1);
    g.stroke();
  }

  _xLabels(g, sx, xt, title) {
    g.fillStyle = THEME.textMuted;
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'center';
    g.textBaseline = 'top';
    for (let i = 0; i < xt.n; i++) {
      const px = sx.to(xt.v[i]);
      if (px < this.px0 - 1 || px > this.px1 + 1) continue;
      g.fillText(xt.lab[i], px, this.py0 + 4);
    }
    if (title) {
      g.textAlign = 'right';
      g.fillStyle = THEME.textSecondary;
      g.font = THEME.fontLabel;
      g.fillText(title, this.px1, this.py0 + 13);
    }
  }

  _yLabels(g, sy, yt, title, side, color) {
    g.fillStyle = color || THEME.textMuted;
    g.font = THEME.fontMonoSmall;
    g.textAlign = side === 'right' ? 'left' : 'right';
    g.textBaseline = 'middle';
    const xr = side === 'right' ? this.px1 + 5 : this.px0 - 5;
    for (let i = 0; i < yt.n; i++) {
      const py = sy.to(yt.v[i]);
      if (py < this.py1 - 1 || py > this.py0 + 1) continue;
      g.fillText(yt.lab[i], xr, py);
    }
    if (title) {
      g.save();
      // rotate(-90°) переводит локальную ось +y в экранную +x, поэтому обе подписи
      // рисуются с baseline 'bottom' (глифы уходят влево от точки привязки).
      g.translate(side === 'right' ? this.w - 2 : 10, (this.py0 + this.py1) / 2);
      g.rotate(-Math.PI / 2);
      g.textAlign = 'center';
      g.textBaseline = 'bottom';
      g.fillStyle = color || THEME.textSecondary;
      g.font = THEME.fontLabel;
      g.fillText(title, 0, 0);
      g.restore();
    }
  }

  _titleBar(g, text, right) {
    g.font = THEME.fontLabel;
    g.textBaseline = 'top';
    if (text) {
      g.fillStyle = THEME.textSecondary;
      g.textAlign = 'left';
      g.fillText(text, this.px0, 2);
    }
    if (right) {
      g.fillStyle = THEME.textMuted;
      g.textAlign = 'right';
      g.font = THEME.fontMonoSmall;
      g.fillText(right, this.px1, 2);
    }
  }

  _badge(g, text, color) {
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    // бейдж не имеет права вылезать за поле графика — обрезаем с многоточием.
    // Кэш: строка предупреждения меняется редко, в кадре ничего не считается.
    const avail = this.px1 - this.px0 - 12;
    if (this._badgeSrc !== text || this._badgeW !== avail) {
      this._badgeSrc = text; this._badgeW = avail;
      let s = text;
      if (g.measureText(s).width > avail) {
        while (s.length > 4 && g.measureText(s + '…').width > avail) s = s.slice(0, -1);
        s += '…';
      }
      this._badgeTxt = s;
    }
    text = this._badgeTxt;
    const wpx = g.measureText(text).width + 10;
    const x = this.px0 + 4, y = this.py1 + 4;
    g.fillStyle = 'rgba(10,13,18,0.85)';
    g.fillRect(x, y, wpx, 15);
    g.strokeStyle = color;
    g.lineWidth = 1;
    g.strokeRect(x + 0.5, y + 0.5, wpx - 1, 14);
    g.fillStyle = color;
    g.fillText(text, x + 5, y + 3);
  }
}

// ───────────────────────────────────── бинарный поиск по времени

function lowerBound(arr, n, t) {
  let lo = 0, hi = n;
  while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < t) lo = m + 1; else hi = m; }
  return lo;
}

// ══════════════════════════════════════════════════════ WaveformPlot

const WF_TRACES = [
  { key: 'Uapp', label: 'U_app', color: THEME.Uapp, width: 1.6, dash: null, axis: 'U' },
  { key: 'Ugap', label: 'U_gap', color: THEME.Ugap, width: 1.4, dash: [6, 3], axis: 'U' },
  { key: 'Itot', label: 'I_total', color: THEME.Itot, width: 1.0, dash: [3, 3], axis: 'I' },
  { key: 'Icond', label: 'I_cond', color: THEME.Icond, width: 1.6, dash: null, axis: 'I' },
];

/**
 * Осциллограмма U/I на общей оси времени.
 *
 * Левая ось — напряжение в кВ (линейная, симметричная).
 * Правая ось — ток в мА (lin | log по opts.logCurrent).
 *
 * Декимация — min/max по пиксельным колонкам (НЕ подвыборка): импульс
 * длительностью 20 нс в окне периода 100 мкс не имеет права исчезнуть.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {object} series  колоночные ряды из loader.mjs (t, Uapp, Ugap, Icond, Idisp, Itot, ...)
 * @param {object} opts    { logCurrent, showUgap, showItot, showIcond, pad, onSeek }
 */
export class WaveformPlot extends BasePlot {
  constructor(canvas, series, opts = {}) {
    super(canvas, Object.assign({ pad: { l: 52, r: 56, t: 16, b: 28 } }, opts));
    this.yI = new Scale();
    this.ytI = new Ticks();
    this.xs = new Scale(); this.ys = new Scale(); this.yIs = new Scale();
    this.logCurrent = opts.logCurrent ?? true;
    this.show = {
      Uapp: true,
      Ugap: opts.showUgap ?? true,
      Itot: opts.showItot ?? true,
      Icond: opts.showIcond ?? true,
    };
    this.decades = opts.decades || 6;
    this._minB = null; this._maxB = null; this._negB = null; this._colBuf = 0;
    this.setSeries(series);
    if (opts.onSeek) {
      this._onClick = (ev) => {
        const r = this.canvas.getBoundingClientRect();
        const px = ev.clientX - r.left;
        if (px < this.px0 || px > this.px1) return;
        opts.onSeek(this.x.from(px));
      };
      canvas.addEventListener('click', this._onClick);
    }
    this.onResize();
  }

  setSeries(series) {
    this.series = series || null;
    this._scKey = '';
    this.t = series ? series.t : null;
    this.nS = this.t ? this.t.length : 0;
    if (this.nS > 0) this.setWindow(this.t[0], this.t[this.nS - 1]);
    return this;
  }

  onResize() {
    const cols = Math.max(2, Math.ceil(this.px1 - this.px0) + 2);
    if (cols !== this._colBuf) { this._colBuf = cols; this._scKey = ''; }
  }

  /** Окно времени (с). Значения зажимаются в доступный диапазон ряда. */
  setWindow(t0, t1) {
    if (!this.nS) return this;
    const lo = this.t[0], hi = this.t[this.nS - 1];
    if (!(t1 > t0)) { t0 = lo; t1 = hi; }
    this.tw0 = Math.max(lo, Math.min(t0, hi));
    this.tw1 = Math.min(hi, Math.max(t1, this.tw0 + 1e-12));
    this.requestRender();
    return this;
  }

  /** Курсор времени (с). NaN — скрыть. */
  setCursor(t) { this.cursorT = t; this.requestRender(); return this; }

  setLogCurrent(on) { this.logCurrent = !!on; this.requestRender(); return this; }
  toggle(key, on) { this.show[key] = on ?? !this.show[key]; this.requestRender(); return this; }

  /** Пиксель -> время (для host-обработчиков). */
  timeAt(px) { return this.x.from(px); }

  /**
   * Декимация окна во ВСЕ показываемые кривые сразу, с кэшем.
   *
   * Раньше каждый вызов draw() делал 8 проходов по ряду (4 на домены + 4 на
   * колонки). На `run-default` в полных данных ряд — 491 653 отсчёта, и один
   * кадр стоил 72 мс: движение курсора (кадр не менялся!) роняло плеер до
   * 13 fps. Теперь проход один на кривую и только при смене окна/размера/
   * набора кривых; курсор пересчёта не вызывает. Сама декимация min/max —
   * прежняя, бит в бит (импульс 20 нс в окне 100 мкс не теряется).
   */
  _prepare(i0, i1) {
    const cols = this._colBuf;
    let showKey = '';
    for (const tr of WF_TRACES) showKey += (this.show[tr.key] && this.series[tr.key]) ? '1' : '0';
    const key = `${this.tw0}|${this.tw1}|${cols}|${this.logCurrent ? 1 : 0}|${showKey}|${this.nS}`;
    if (this._scKey === key) return this._scStats;
    this._scKey = key;
    if (!this._sc) this._sc = new Map();
    let uMax = 0, iMax = 0, iMinPos = Infinity;
    for (const tr of WF_TRACES) {
      if (!this.show[tr.key] || !this.series[tr.key]) continue;
      let e = this._sc.get(tr.key);
      if (!e || e.min.length !== cols) {
        e = { min: new Float32Array(cols), max: new Float32Array(cols), neg: new Uint8Array(cols) };
        this._sc.set(tr.key, e);
      }
      const arr = this.series[tr.key];
      const abs = tr.axis === 'I' && this.logCurrent;
      e.min.fill(Infinity); e.max.fill(-Infinity); e.neg.fill(0);
      e.any = false; e.gMax = 0; e.gMinPos = Infinity;
      const sx = (cols - 1) / (this.tw1 - this.tw0);
      for (let k = i0; k <= i1; k++) {
        let c = ((this.t[k] - this.tw0) * sx) | 0;
        if (c < 0) c = 0; else if (c >= cols) c = cols - 1;
        let v = arr[k];
        if (!Number.isFinite(v)) continue;
        const a = v < 0 ? -v : v;
        if (a > e.gMax) e.gMax = a;
        if (a > 0 && a < e.gMinPos) e.gMinPos = a;
        if (abs) { if (v < 0) { e.neg[c] = 1; v = a; } }
        if (v < e.min[c]) e.min[c] = v;
        if (v > e.max[c]) e.max[c] = v;
        e.any = true;
      }
      if (tr.axis === 'U') { if (e.gMax > uMax) uMax = e.gMax; }
      else { if (e.gMax > iMax) iMax = e.gMax; if (e.gMinPos < iMinPos) iMinPos = e.gMinPos; }
    }
    this._scStats = { uMax, iMax, iMinPos };
    return this._scStats;
  }

  _drawEnvelope(g, sy, color, width, dash, e) {
    const min = e.min, max = e.max, cols = this._colBuf;
    const sx = (this.px1 - this.px0) / (cols - 1);
    g.beginPath();
    let open = false;
    for (let c = 0; c < cols; c++) {
      if (min[c] === Infinity) { open = false; continue; }
      const px = this.px0 + c * sx;
      const a = sy.to(max[c]);
      let b = sy.to(min[c]);
      // лог-ось: нулевой минимум колонки не должен убивать валидный максимум
      if (!Number.isFinite(a)) { open = false; continue; }
      if (!Number.isFinite(b)) b = a;
      if (!open) { g.moveTo(px, a); open = true; }
      else g.lineTo(px, a);
      if (b !== a) g.lineTo(px, b);
    }
    g.strokeStyle = color;
    g.lineWidth = width;
    if (dash) g.setLineDash(dash); else g.setLineDash(EMPTY_DASH);
    g.stroke();
    g.setLineDash(EMPTY_DASH);
  }

  _drawNegStrip(g, color, e) {
    const neg = e.neg, min = e.min, cols = this._colBuf;
    const sx = (this.px1 - this.px0) / (cols - 1);
    g.fillStyle = color;
    for (let c = 0; c < cols; c++) {
      if (!neg[c] || min[c] === Infinity) continue;
      g.fillRect(this.px0 + c * sx, this.py0 - 3, Math.max(1, sx), 3);
    }
  }

  draw(g) {
    if (!this.nS) { drawEmpty(g, this, 'нет series.bin'); return; }
    const i0 = Math.max(0, lowerBound(this.t, this.nS, this.tw0) - 1);
    const i1 = Math.min(this.nS - 1, lowerBound(this.t, this.nS, this.tw1));

    // ── масштабы
    // yI — вторая ось, _measure про неё не знает: пиксели ставим здесь,
    // ДО setDomain (setDomain пересчитывает k по текущим p0/p1).
    this.yI.setPixels(this.py0, this.py1);
    const tScale = pickTimeScale(this.tw1 - this.tw0);
    this.x.setDomain(this.tw0, this.tw1, false);

    const st = this._prepare(i0, i1);
    let uMax = st.uMax;
    if (!(uMax > 0)) uMax = 1;
    this.y.setDomain(-uMax * 1.08, uMax * 1.08, false);

    let iMax = st.iMax;
    const iAbsMinPos = st.iMinPos;
    if (!(iMax > 0)) iMax = 1e-6;
    if (this.logCurrent) {
      const hi = Math.pow(10, Math.ceil(Math.log10(iMax)));
      const loFloor = hi * Math.pow(10, -this.decades);
      const lo = Math.min(loFloor, Math.max(iAbsMinPos, hi * 1e-12));
      this.yI.setDomain(Math.max(lo, hi * 1e-12), hi, true);
    } else {
      this.yI.setDomain(-iMax * 1.08, iMax * 1.08, false);
    }

    this.xt.linear(this.tw0 / tScale.k, this.tw1 / tScale.k, Math.max(3, (this.px1 - this.px0) / 90));
    const xs = dispScale(this.xs, this.x, tScale.k);
    this.yt.linear(this.y.min / 1e3, this.y.max / 1e3, Math.max(3, (this.py0 - this.py1) / 34));
    const ys = dispScale(this.ys, this.y, 1e3);
    if (this.logCurrent) this.ytI.log(this.yI.min, this.yI.max, Math.max(2, (this.py0 - this.py1) / 34));
    else this.ytI.linear(this.yI.min * 1e3, this.yI.max * 1e3, Math.max(3, (this.py0 - this.py1) / 34));

    this._grid(g, xs, ys, this.xt, this.yt);

    // ── нулевая линия напряжения
    const y0 = Math.round(this.y.to(0)) + 0.5;
    if (y0 > this.py1 && y0 < this.py0) {
      g.strokeStyle = THEME.axis;
      g.setLineDash(DASH_2_2);
      g.beginPath(); g.moveTo(this.px0, y0); g.lineTo(this.px1, y0); g.stroke();
      g.setLineDash(EMPTY_DASH);
    }

    this._clipPlot(g);
    // ток рисуем первым — напряжение важнее и должно быть сверху
    for (const tr of WF_TRACES) {
      if (tr.axis !== 'I' || !this.show[tr.key] || !this.series[tr.key]) continue;
      const e = this._sc.get(tr.key);
      if (!e || !e.any) continue;
      this._drawEnvelope(g, this.yI, tr.color, tr.width, tr.dash, e);
      if (this.logCurrent && tr.key === 'Icond') this._drawNegStrip(g, 'rgba(57,135,229,0.85)', e);
    }
    for (const tr of WF_TRACES) {
      if (tr.axis !== 'U' || !this.show[tr.key] || !this.series[tr.key]) continue;
      const e = this._sc.get(tr.key);
      if (!e || !e.any) continue;
      this._drawEnvelope(g, this.y, tr.color, tr.width, tr.dash, e);
    }
    g.restore();

    // ── оси
    this._xLabels(g, xs, this.xt, `t, ${tScale.unit}`);
    this._yLabels(g, ys, this.yt, 'U, kV', 'left', THEME.textSecondary);
    const yIs = this.logCurrent ? this.yI : dispScale(this.yIs, this.yI, 1e3);
    this._yLabels(g, yIs, this.ytI, this.logCurrent ? '|I|, A (log)' : 'I, mA', 'right', THEME.Icond);

    // ── курсор + считывание
    if (Number.isFinite(this.cursorT) && this.cursorT >= this.tw0 && this.cursorT <= this.tw1) {
      const px = Math.round(this.x.to(this.cursorT)) + 0.5;
      g.strokeStyle = THEME.cursor;
      g.lineWidth = 1;
      g.beginPath(); g.moveTo(px, this.py1); g.lineTo(px, this.py0); g.stroke();
      const k = Math.min(this.nS - 1, lowerBound(this.t, this.nS, this.cursorT));
      this._readout(g, k);
    }
    this._legend(g);
    if (this.logCurrent) {
      g.font = THEME.fontMonoSmall;
      g.fillStyle = THEME.textMuted;
      g.textAlign = 'left'; g.textBaseline = 'bottom';
      g.fillText('log: |I|, синяя полоса = I < 0', this.px0 + 3, this.py0 - 4);
    }
  }

  _readout(g, k) {
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'right';
    g.textBaseline = 'top';
    let y = this.py1 + 3;
    g.fillStyle = THEME.textMuted;
    g.fillText(fmtTime(this.t[k]), this.px1 - 3, y); y += 12;
    for (const tr of WF_TRACES) {
      if (!this.show[tr.key] || !this.series[tr.key]) continue;
      g.fillStyle = tr.color;
      const v = this.series[tr.key][k];
      g.fillText(`${tr.label} ${tr.axis === 'U' ? fmtSI(v, 'V', 4) : fmtSI(v, 'A', 3)}`, this.px1 - 3, y);
      y += 12;
    }
  }

  _legend(g) {
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    if (!this._legW) {                       // ширины подписей меряются один раз
      this._legW = new Float64Array(WF_TRACES.length);
      for (let i = 0; i < WF_TRACES.length; i++) this._legW[i] = g.measureText(WF_TRACES[i].label).width;
    }
    let x = this.px0 + 3;
    const y = this.py1 + 3;
    for (let i = 0; i < WF_TRACES.length; i++) {
      const tr = WF_TRACES[i];
      if (!this.show[tr.key] || !this.series[tr.key]) continue;
      g.strokeStyle = tr.color;
      g.lineWidth = tr.width;
      if (tr.dash) g.setLineDash(tr.dash); else g.setLineDash(EMPTY_DASH);
      g.beginPath(); g.moveTo(x, y + 6); g.lineTo(x + 14, y + 6); g.stroke();
      g.setLineDash(EMPTY_DASH);
      g.fillStyle = tr.color;
      g.fillText(tr.label, x + 17, y);
      x += 27 + this._legW[i];
    }
  }

  destroy() {
    if (this._onClick) this.canvas.removeEventListener('click', this._onClick);
    super.destroy();
  }
}

const EMPTY_DASH = [];
const DASH_2_2 = [2, 2];
const DASH_4_4 = [4, 4];
const DASH_6_3 = [6, 3];

/**
 * «Псевдоскейл» для подписей: те же пиксели, но домен поделён на k (В->кВ, с->мкс).
 * Пишет в переданный экземпляр — ни одной аллокации в кадре.
 */
function dispScale(dst, src, k) {
  dst.log = src.log;
  dst.min = src.min / k; dst.max = src.max / k;
  dst.p0 = src.p0; dst.p1 = src.p1;
  dst._upd();
  return dst;
}

function pickTimeScale(span) {
  if (span >= 1e-1) return TS_S;
  if (span >= 1e-4) return TS_MS;
  if (span >= 1e-7) return TS_US;
  return TS_NS;
}
const TS_S = { k: 1, unit: 's' };
const TS_MS = { k: 1e-3, unit: 'ms' };
const TS_US = { k: 1e-6, unit: 'µs' };
const TS_NS = { k: 1e-9, unit: 'ns' };

function drawEmpty(g, plot, msg) {
  g.fillStyle = THEME.textMuted;
  g.font = THEME.fontLabel;
  g.textAlign = 'center';
  g.textBaseline = 'middle';
  g.fillText(msg, (plot.px0 + plot.px1) / 2, (plot.py0 + plot.py1) / 2);
}

// ═════════════════════════════════════════════════════ LissajousPlot

/**
 * Фигура Лиссажу Q–U и извлечение из неё ёмкостей (метод Мэнли).
 *
 * ЧТО ИМЕННО ДЕЛАЕТСЯ (важно, потому что «наклон Лиссажу» часто считают неверно):
 *  1. Локальные наклоны dQ/dU по соседним отсчётам, вес |dU| (короткие сегменты
 *     около разворота напряжения не должны тянуть фит).
 *  2. Кластеризация наклонов на ДВА кластера (k-means, 1D, старт с 20/80 перцентилей
 *     — НЕ с аналитических значений, иначе фит перестал бы быть независимым).
 *  3. Разделение каждого кластера на две ветви по знаку dU/dt -> четыре стороны.
 *  4. Взвешенный МНК по каждой стороне отдельно; C_cell — среднее наклонов
 *     «ёмкостных» сторон, C_diel — «разрядных».
 *  5. Сравнение с аналитикой из геометрии (metrics.analyticCapacitances) и вывод
 *     расхождения в процентах. Аналитика: C_cell 0.797 пФ/см², C_diel 7.97 пФ/см²
 *     (ERRATA §C V4), пересчёт на площадь ячейки pi*R^2 — внутри.
 *  6. Бейдж качества: если петля не замкнута, сторон меньше четырёх, R² < 0.98
 *     или кластеры не разделены (C_diel/C_cell < 1.5) — наклоны помечаются
 *     как недостоверные. Молча выдать C_diel по кривой петле хуже, чем не выдать.
 *
 * Энергия за период W = ∮U dQ считается ЧИСЛЕННО по данным (не по фиту),
 * мощность P = W/T. U_min получается обращением формулы Мэнли
 * W = 4·C_diel·U_min·(U0 − U_min) и помечается как производная величина.
 *
 * @param {object} opts { freqHz, geometry:{epsR,gapMM,dielMM,radiusMM}, analytic:{C_cell,C_diel},
 *                        mode:'period'|'last8'|'all', ghosts:number }
 */
export class LissajousPlot extends BasePlot {
  constructor(canvas, series, opts = {}) {
    super(canvas, Object.assign({ pad: { l: 58, r: 14, t: 16, b: 26 } }, opts));
    this.series = series || null;
    this.t = series ? series.t : null;
    this.nS = this.t ? this.t.length : 0;
    this.freqHz = opts.freqHz || 1e4;
    this.T = 1 / this.freqHz;
    this.analytic = opts.analytic
      || (opts.geometry ? analyticCapacitances(opts.geometry) : null);
    this.mode = opts.mode || 'period';
    this.ghosts = opts.ghosts ?? 7;
    this._fits = new Map();     // periodIndex -> fit (аллокация раз на период)
    this.xs = new Scale(); this.ys = new Scale();
    this._rowY = 0; this._rowX = 0;   // состояние построчного вывода отчёта
    this.cursorT = this.nS ? this.t[this.nS - 1] : NaN;
  }

  setCursor(t) { this.cursorT = t; this.requestRender(); return this; }
  setMode(m) { this.mode = m; this.requestRender(); return this; }
  /** Результат фита текущего периода (или null). Для MetricsPanel. */
  get fit() { return this._fitFor(this._periodOf(this.cursorT)); }

  _periodOf(t) {
    if (!Number.isFinite(t)) return 0;
    return Math.floor(t / this.T);
  }

  _rangeOfPeriod(p) {
    const t0 = p * this.T, t1 = (p + 1) * this.T;
    const i0 = lowerBound(this.t, this.nS, t0);
    const i1 = Math.min(this.nS - 1, lowerBound(this.t, this.nS, t1));
    return { i0, i1, t0, t1 };
  }

  _fitFor(p) {
    if (!this.nS) return null;
    const cached = this._fits.get(p);
    if (cached) return cached;
    const { i0, i1 } = this._rangeOfPeriod(p);
    const f = fitLissajous(this.series.Uapp, this.series.Q, i0, i1, this.T, this.analytic);
    this._fits.set(p, f);
    if (this._fits.size > 32) this._fits.delete(this._fits.keys().next().value);
    return f;
  }

  draw(g) {
    if (!this.nS) { drawEmpty(g, this, 'нет series.bin'); return; }
    if (!this.series.Q || !this.series.Uapp) { drawEmpty(g, this, 'в series нет Q/Uapp'); return; }
    const p = this._periodOf(this.cursorT);
    const fit = this._fitFor(p);

    // квадратная область построения: наклоны читаются верно только при 1:1
    const availW = this.px1 - this.px0, availH = this.py0 - this.py1;
    const side = Math.min(availW, availH);
    const cx0 = this.px0, cx1 = this.px0 + side;
    const cy0 = this.py0, cy1 = this.py0 - side;
    const savedPx1 = this.px1, savedPy1 = this.py1;
    this.px1 = cx1; this.py1 = cy1;
    this.x.setPixels(cx0, cx1);
    this.y.setPixels(cy0, cy1);

    // домены: по всем показываемым периодам
    const pFirst = this.mode === 'all' ? 0
      : this.mode === 'last8' ? Math.max(0, p - this.ghosts) : p;
    let uMax = 0, qMin = Infinity, qMax = -Infinity;
    for (let q = pFirst; q <= p; q++) {
      const r = this._rangeOfPeriod(q);
      for (let k = r.i0; k <= r.i1; k++) {
        const u = Math.abs(this.series.Uapp[k]); if (u > uMax) uMax = u;
        const qq = this.series.Q[k];
        if (qq < qMin) qMin = qq;
        if (qq > qMax) qMax = qq;
      }
    }
    if (!(uMax > 0)) uMax = 1;
    if (!(qMax > qMin)) { qMax = qMin + 1e-15; }
    const qPad = (qMax - qMin) * 0.08;
    this.x.setDomain(-uMax * 1.08, uMax * 1.08, false);
    this.y.setDomain(qMin - qPad, qMax + qPad, false);

    const uk = 1e3;                                  // В -> кВ
    const qk = pickChargeScale(Math.max(Math.abs(qMin), Math.abs(qMax)));
    this.xt.linear(this.x.min / uk, this.x.max / uk, 5);
    this.yt.linear(this.y.min / qk.k, this.y.max / qk.k, Math.max(3, side / 36));
    const xs = dispScale(this.xs, this.x, uk);
    const ys = dispScale(this.ys, this.y, qk.k);
    this._grid(g, xs, ys, this.xt, this.yt);

    this._clipPlot(g);
    // призраки предыдущих периодов — без цвета, только прозрачность (UI_SPEC §2.3)
    for (let q = pFirst; q < p; q++) {
      const a = 0.10 + 0.25 * (1 - (p - q) / Math.max(1, p - pFirst + 1));
      this._path(g, q, `rgba(232,234,237,${a.toFixed(3)})`, 1);
    }
    this._path(g, p, THEME.text, 2);

    // фит-линии
    if (fit && fit.sides) {
      g.setLineDash(DASH_4_4);
      g.lineWidth = 1;
      for (let s = 0; s < 4; s++) {
        const sd = fit.sides[s];
        if (!sd || !sd.ok) continue;
        g.strokeStyle = sd.cluster === 1 ? 'rgba(255,209,102,0.75)' : 'rgba(139,147,161,0.9)';
        g.beginPath();
        const u0 = this.x.min, u1 = this.x.max;
        g.moveTo(this.x.to(u0), this.y.to(sd.a * u0 + sd.b));
        g.lineTo(this.x.to(u1), this.y.to(sd.a * u1 + sd.b));
        g.stroke();
      }
      g.setLineDash(EMPTY_DASH);
      if (fit.vertices) {
        g.fillStyle = THEME.accent;
        for (let v = 0; v < fit.vertices.length; v++) {
          const V = fit.vertices[v];
          if (!V) continue;
          g.beginPath();
          g.arc(this.x.to(V.u), this.y.to(V.q), 3, 0, 6.283185307179586);
          g.fill();
        }
      }
    }

    // движущаяся точка
    if (Number.isFinite(this.cursorT)) {
      const k = Math.min(this.nS - 1, lowerBound(this.t, this.nS, this.cursorT));
      const px = this.x.to(this.series.Uapp[k]), py = this.y.to(this.series.Q[k]);
      g.fillStyle = THEME.cursor;
      g.beginPath(); g.arc(px, py, 4, 0, 6.283185307179586); g.fill();
      g.strokeStyle = THEME.bg; g.lineWidth = 1.5; g.stroke();
    }
    g.restore();

    this._xLabels(g, xs, this.xt, 'U_app, kV');
    this._yLabels(g, ys, this.yt, `Q, ${qk.unit}`, 'left', THEME.textSecondary);

    // ── правая колонка чисел
    this._report(g, fit, cx1 + 8, p);

    this.px1 = savedPx1; this.py1 = savedPy1;
    this.x.setPixels(this.px0, this.px1);
    this.y.setPixels(this.py0, this.py1);

    if (fit && fit.warn) this._badge(g, '⚠ ' + fit.warn, THEME.warning);
  }

  /**
   * Ломаная периода с ПИКСЕЛЬНОЙ декимацией: подряд идущие отсчёты, попадающие
   * в тот же полупиксель, в путь не добавляются.
   *
   * Это не подвыборка по индексу (она срезала бы углы фигуры), а отбрасывание
   * геометрически неразличимых точек: соседняя точка либо сдвинулась хотя бы на
   * полпикселя — и рисуется, либо не сдвинулась — и не видна вовсе. Последняя
   * точка периода добавляется всегда, иначе ветвь визуально не доходит до
   * вершины. Мера: на `run-default` (полные данные, 491 653 отсчёта) один
   * render стоил 1020 мс — плеер шёл 1.3 fps; после декимации ~6 мс.
   */
  _path(g, p, color, width) {
    const r = this._rangeOfPeriod(p);
    if (r.i1 <= r.i0) return;
    const U = this.series.Uapp, Q = this.series.Q;
    g.beginPath();
    let lx = this.x.to(U[r.i0]), ly = this.y.to(Q[r.i0]);
    g.moveTo(lx, ly);
    for (let k = r.i0 + 1; k <= r.i1; k++) {
      const px = this.x.to(U[k]), py = this.y.to(Q[k]);
      const dx = px - lx, dy = py - ly;
      if (k !== r.i1 && dx > -0.5 && dx < 0.5 && dy > -0.5 && dy < 0.5) continue;
      g.lineTo(px, py);
      lx = px; ly = py;
    }
    g.strokeStyle = color;
    g.lineWidth = width;
    g.stroke();
  }

  _line(g, label, value, color) {
    g.fillStyle = THEME.textMuted;
    g.fillText(label, this._rowX, this._rowY);
    g.fillStyle = color || THEME.text;
    g.fillText(value, this._rowX, this._rowY + 11);
    this._rowY += 25;
  }

  _report(g, fit, x, p) {
    if (x > this.w - 64) return;   // узкая панель — числа отдаёт MetricsPanel
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    this._rowX = x;
    this._rowY = this.py1 + 2;
    this._line(g, 'period #', numStr(p));
    const an = this.analytic;
    if (fit && fit.ok) {
      const dCell = an ? 100 * (fit.C_cell / an.C_cell - 1) : NaN;
      const dDiel = an ? 100 * (fit.C_diel / an.C_diel - 1) : NaN;
      this._line(g, 'C_cell fit', fmtSI(fit.C_cell, 'F'), devColor(dCell, 4, 25));
      if (an) this._line(g, '  vs аналит.', fmtSI(an.C_cell, 'F') + '  ' + fmtPct(dCell), devColor(dCell, 4, 25));
      this._line(g, 'C_diel fit', fmtSI(fit.C_diel, 'F'), devColor(dDiel, 4, 25));
      if (an) this._line(g, '  vs аналит.', fmtSI(an.C_diel, 'F') + '  ' + fmtPct(dDiel), devColor(dDiel, 4, 25));
      this._line(g, 'R² сторон', fmtFixed(fit.r2min, 3), fit.r2min < 0.98 ? THEME.warning : THEME.text);
    } else {
      this._line(g, 'C_cell fit', '—', THEME.textMuted);
      this._line(g, 'C_diel fit', '—', THEME.textMuted);
    }
    if (fit) {
      this._line(g, 'W = ∮U dQ', fmtSI(fit.W, 'J'));
      this._line(g, 'P = f·W', fmtSI(fit.P, 'W'));
      this._line(g, 'ΔQ', fit.dQ != null ? fmtSI(fit.dQ, 'C') : '—');
      this._line(g, 'U_min (Manley)', fit.Umin != null ? fmtSI(fit.Umin, 'V') : '—',
        fit.Umin == null ? THEME.textMuted : THEME.text);
      this._line(g, 'замыкание', fit.closure != null ? fmtPct(100 * fit.closure, 2) : '—',
        fit.closure != null && Math.abs(fit.closure) > 0.05 ? THEME.warning : THEME.text);
    }
  }
}

function devColor(pct, warn, crit) {
  if (!Number.isFinite(pct)) return THEME.textMuted;
  const a = Math.abs(pct);
  if (a > crit) return THEME.critical;
  if (a > warn) return THEME.warning;
  return THEME.good;
}

/** Кэш строк маленьких целых — чтобы String(n) не аллоцировал каждый кадр. */
const NUM_CACHE = [];
function numStr(n) {
  if (n >= 0 && n < 4096 && (n | 0) === n) {
    let s = NUM_CACHE[n];
    if (s === undefined) { s = String(n); NUM_CACHE[n] = s; }
    return s;
  }
  return String(n);
}

function pickChargeScale(q) {
  if (!(q > 0)) return { k: 1e-12, unit: 'pC' };
  if (q >= 1e-6) return { k: 1e-6, unit: 'µC' };
  if (q >= 1e-9) return { k: 1e-9, unit: 'nC' };
  if (q >= 1e-12) return { k: 1e-12, unit: 'pC' };
  return { k: 1e-15, unit: 'fC' };
}

/** Взвешенная медиана значений v[idx[0..n)]. */
function weightedMedian(v, w, idx, n) {
  const ord = Array.prototype.slice.call(idx.subarray(0, n));
  ord.sort((a, b) => v[a] - v[b]);
  let W = 0;
  for (let i = 0; i < n; i++) W += w[ord[i]];
  let acc = 0;
  for (let i = 0; i < n; i++) { acc += w[ord[i]]; if (acc >= 0.5 * W) return v[ord[i]]; }
  return v[ord[n - 1]];
}

/** Взвешенная медиана выражения (q - a*u) — свободный член прямой заданного наклона. */
function weightedMedianExpr(q, u, a, w, idx, n) {
  const tmp = new Float64Array(n);
  const ord = new Int32Array(n);
  for (let i = 0; i < n; i++) { tmp[i] = q[idx[i]] - a * u[idx[i]]; ord[i] = i; }
  const arr = Array.prototype.slice.call(ord);
  arr.sort((x, y) => tmp[x] - tmp[y]);
  let W = 0;
  for (let i = 0; i < n; i++) W += w[idx[i]];
  let acc = 0;
  for (let i = 0; i < n; i++) { acc += w[idx[arr[i]]]; if (acc >= 0.5 * W) return tmp[arr[i]]; }
  return tmp[arr[n - 1]];
}

/** Размах остатков (95-й минус 5-й перцентиль) — устойчивая мера «ветвь не прямая». */
function residualSpread(q, u, a, b, idx, n) {
  const r = new Float64Array(n);
  for (let i = 0; i < n; i++) r[i] = q[idx[i]] - (a * u[idx[i]] + b);
  r.sort();
  return r[Math.floor(0.95 * (n - 1))] - r[Math.floor(0.05 * (n - 1))];
}

/**
 * Фит фигуры Лиссажу. Чистая функция — тестируется отдельно от canvas.
 * @returns {{ok,C_cell,C_diel,sides,vertices,W,P,dQ,Umin,r2min,closure,warn}}
 */
export function fitLissajous(U, Q, i0, i1, T, analytic) {
  const out = {
    ok: false, C_cell: NaN, C_diel: NaN, sides: null, vertices: null,
    W: NaN, P: NaN, dQ: null, Umin: null, r2min: NaN, cellSpread: NaN, closure: null,
    clusterSep: NaN, warn: null, note: null,
  };
  const n = i1 - i0 + 1;
  if (n < 16) { out.warn = 'период не набран: точек < 16'; return out; }

  // энергия и замыкание — по данным, без фита
  let W = 0, uMax = 0, qMin = Infinity, qMax = -Infinity;
  for (let k = i0; k < i1; k++) {
    W += 0.5 * (U[k] + U[k + 1]) * (Q[k + 1] - Q[k]);
    const u = Math.abs(U[k]); if (u > uMax) uMax = u;
    if (Q[k] < qMin) qMin = Q[k];
    if (Q[k] > qMax) qMax = Q[k];
  }
  out.W = Math.abs(W);
  out.P = out.W / T;
  const qSpan = qMax - qMin;
  out.closure = qSpan > 0 ? (Q[i1] - Q[i0]) / qSpan : null;

  // локальные наклоны.
  // Вес сегмента — его ДЛИНА В НОРМИРОВАННЫХ КООРДИНАТАХ (dU/Uspan, dQ/Qspan).
  // Вес |dU| был бы систематически несправедлив к разрядной ветви: она проходится
  // почти при постоянном U, и её сегменты получили бы вес ~0.
  const m = n - 1;
  const s = new Float64Array(m), w = new Float64Array(m), dir = new Int8Array(m);
  const uMid = new Float64Array(m), qMid = new Float64Array(m);
  const eps = Math.max(1e-12, uMax * 1e-5);
  const uSpan = Math.max(2 * uMax, 1e-30), qs = Math.max(qSpan, 1e-30);
  let cnt = 0;
  for (let k = i0; k < i1; k++) {
    const du = U[k + 1] - U[k];
    if (Math.abs(du) < eps) continue;
    const dq = Q[k + 1] - Q[k];
    const sl = dq / du;
    if (!Number.isFinite(sl) || sl <= 0) continue;  // отрицательный наклон — не сторона параллелограмма
    s[cnt] = sl;
    w[cnt] = Math.hypot(du / uSpan, dq / qs);
    dir[cnt] = du > 0 ? 1 : -1;
    uMid[cnt] = 0.5 * (U[k] + U[k + 1]);
    qMid[cnt] = 0.5 * (Q[k] + Q[k + 1]);
    cnt++;
  }
  if (cnt < 16) { out.warn = 'нет пригодных сегментов dQ/dU'; return out; }

  // k-means ПО ЛОГАРИФМУ наклона: ёмкостная и разрядная ветви различаются
  // в 10..100 раз, в линейной метрике доминирующий кластер съедает второй
  // (проверено на run-synth: линейный k-means давал C_cell ×7.8).
  // Старт — 5/95 перцентили, аналитика в инициализацию НЕ входит.
  const ls = new Float64Array(cnt);
  for (let k = 0; k < cnt; k++) ls[k] = Math.log10(s[k]);
  const sorted = Float64Array.from(ls).sort();
  let c0 = sorted[Math.floor(0.05 * (cnt - 1))], c1 = sorted[Math.floor(0.95 * (cnt - 1))];
  if (!(c1 > c0)) c1 = c0 + 0.05;
  const lab = new Uint8Array(cnt);
  for (let it = 0; it < 50; it++) {
    let s0 = 0, n0 = 0, s1 = 0, n1 = 0, moved = 0;
    for (let k = 0; k < cnt; k++) {
      const l = Math.abs(ls[k] - c1) < Math.abs(ls[k] - c0) ? 1 : 0;
      if (l !== lab[k]) { lab[k] = l; moved++; }
      if (l) { s1 += ls[k]; n1++; } else { s0 += ls[k]; n0++; }
    }
    if (n0 > 0) c0 = s0 / n0;
    if (n1 > 0) c1 = s1 / n1;
    if (!moved && it > 1) break;
  }
  if (!(c1 > c0)) { const t2 = c0; c0 = c1; c1 = t2; for (let k = 0; k < cnt; k++) lab[k] ^= 1; }
  out.clusterSep = Math.pow(10, c1 - c0);

  // ── четыре стороны: кластер × направление dU.
  //
  // ДВА РАЗНЫХ ОЦЕНЩИКА, И ЭТО НЕ ЛЕНЬ, А ФИЗИКА:
  //  * C_cell (ёмкостная ветвь) — ЛОКАЛЬНЫЙ наклон: взвешенная медиана dQ/dU.
  //    Глобальный МНК по точкам этой ветви даёт ерунду, если за полупериод
  //    прошло много микроразрядов: ветвь превращается в лестницу, и прямая,
  //    проведённая по её ступеням, имеет наклон ЛЕСТНИЦЫ (~C_diel), а не ступени.
  //    На run-synth (52 микроразряда за период) глобальный МНК давал 48 фФ
  //    вместо 6.26 фФ — ровно та ошибка, о которой предупреждает VALIDATION §V4b.
  //  * C_diel (разрядная ветвь) — ГЛОБАЛЬНЫЙ взвешенный МНК: здесь нужен
  //    именно наклон огибающей ветви, а не наклон отдельного скачка.
  //
  // Для ёмкостных сторон дополнительно меряется разброс свободного члена b:
  // если ветвь не одна прямая (лестница), это видно как spreadB/Qspan.
  const sides = new Array(4);
  const bufIdx = new Int32Array(cnt);
  for (let cl = 0; cl < 2; cl++) {
    for (let dd = 0; dd < 2; dd++) {
      const want = dd === 0 ? 1 : -1;
      const idx = cl * 2 + dd;
      let np = 0;
      for (let k = 0; k < cnt; k++) if (lab[k] === cl && dir[k] === want) bufIdx[np++] = k;
      if (np < 6) { sides[idx] = { ok: false, cluster: cl, dir: want, n: np }; continue; }
      let a, b, spreadB = 0;
      if (cl === 0) {
        a = weightedMedian(s, w, bufIdx, np);
        b = weightedMedianExpr(qMid, uMid, a, w, bufIdx, np);
        spreadB = residualSpread(qMid, uMid, a, b, bufIdx, np) / Math.max(qSpan, 1e-30);
      } else {
        let sw = 0, su = 0, sq = 0, suu = 0, suq = 0;
        for (let t2 = 0; t2 < np; t2++) {
          const k = bufIdx[t2], ww = w[k];
          sw += ww; su += ww * uMid[k]; sq += ww * qMid[k];
          suu += ww * uMid[k] * uMid[k]; suq += ww * uMid[k] * qMid[k];
        }
        const den = sw * suu - su * su;
        if (!(Math.abs(den) > 0)) { sides[idx] = { ok: false, cluster: cl, dir: want, n: np }; continue; }
        a = (sw * suq - su * sq) / den;
        b = (sq - a * su) / sw;
      }
      // взвешенный R² относительно ИТОГОВОЙ прямой (a,b)
      let ssTot = 0, ssRes = 0, sw2 = 0, sq2 = 0;
      for (let t2 = 0; t2 < np; t2++) { const k = bufIdx[t2]; sw2 += w[k]; sq2 += w[k] * qMid[k]; }
      const qBar = sq2 / sw2;
      for (let t2 = 0; t2 < np; t2++) {
        const k = bufIdx[t2];
        const r = qMid[k] - (a * uMid[k] + b);
        ssRes += w[k] * r * r;
        const d = qMid[k] - qBar;
        ssTot += w[k] * d * d;
      }
      const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
      sides[idx] = { ok: true, cluster: cl, dir: want, a, b, r2, n: np, spreadB };
    }
  }
  out.sides = sides;

  const cellSides = [sides[0], sides[1]].filter((x) => x && x.ok);
  const dielSides = [sides[2], sides[3]].filter((x) => x && x.ok);
  if (!cellSides.length || !dielSides.length) {
    out.warn = 'нет четырёх сторон: петля не параллелограмм';
    return out;
  }
  out.C_cell = cellSides.reduce((acc, x) => acc + x.a, 0) / cellSides.length;
  out.C_diel = dielSides.reduce((acc, x) => acc + x.a, 0) / dielSides.length;
  // R² предъявляется к РАЗРЯДНЫМ сторонам (там он и означает «сторона — прямая»).
  // Качество ёмкостной ветви меряется разбросом свободного члена (см. выше).
  out.r2min = Math.min(...dielSides.map((x) => x.r2));
  out.cellSpread = Math.max(...cellSides.map((x) => x.spreadB || 0));

  // вершины: пересечения соседних сторон
  const inter = (A, B) => {
    if (!A || !B || !A.ok || !B.ok || Math.abs(A.a - B.a) < 1e-30) return null;
    const u = (B.b - A.b) / (A.a - B.a);
    return { u, q: A.a * u + A.b };
  };
  out.vertices = [inter(sides[0], sides[2]), inter(sides[2], sides[1]),
    inter(sides[1], sides[3]), inter(sides[3], sides[0])];

  // перенесённый заряд: размах Q на разрядной ветви между её вершинами
  const v = out.vertices;
  if (v[0] && v[1]) out.dQ = Math.abs(v[1].q - v[0].q);

  // U_min обращением формулы Мэнли: W = 4 C_diel U_min (U0 - U_min)
  const disc = uMax * uMax - out.W / out.C_diel;
  out.Umin = disc >= 0 ? 0.5 * (uMax - Math.sqrt(disc)) : null;

  out.ok = true;
  const warns = [], notes = [];
  if (out.closure != null && Math.abs(out.closure) > 0.05) warns.push('петля не замкнута');
  if (out.r2min < 0.98) warns.push(`R² разрядных сторон = ${out.r2min.toFixed(3)} < 0.98`);
  // информационно, НЕ бракует наклон: медиана локального dQ/dU от дробления не страдает
  if (out.cellSpread > 0.05) {
    notes.push(`ёмкостная ветвь раздроблена микроразрядами (разброс b ${(100 * out.cellSpread).toFixed(0)} % от ΔQ)`);
  }
  if (!(out.clusterSep > 1.5)) warns.push(`ветви не разделены (×${out.clusterSep.toFixed(2)})`);
  if (cellSides.length < 2 || dielSides.length < 2) warns.push('найдено < 4 сторон');
  if (analytic) {
    const d = Math.max(Math.abs(out.C_cell / analytic.C_cell - 1), Math.abs(out.C_diel / analytic.C_diel - 1));
    if (d > 0.25) warns.push(`расхождение с аналитикой ${(100 * d).toFixed(0)} %`);
  }
  if (warns.length) out.warn = warns.join(' · ') + ' — наклоны Мэнли недостоверны';
  out.note = notes.length ? notes.join(' · ') : null;
  return out;
}

// ════════════════════════════════════════════════ RadialProfilePlot

/**
 * sigma(r) на обеих поверхностях диэлектрика.
 * Данные подаются host'ом покадрово: setData(r, sigmaL, sigmaR) — ссылки, без копий.
 * Ось Y: 'lin' (симметричная) либо 'asinh' (знаковая «лог-подобная»), нули корректны.
 */
export class RadialProfilePlot extends BasePlot {
  constructor(canvas, opts = {}) {
    super(canvas, Object.assign({ pad: { l: 62, r: 14, t: 16, b: 28 } }, opts));
    this.r = null; this.sL = null; this.sR = null;
    this.scaleMode = opts.scale || 'lin';       // 'lin' | 'asinh'
    this.unitK = 1e-5;                          // Кл/м² -> нКл/см²  (1 Кл/м² = 1e5 нКл/см²)
    this._smooth = 0;
    this._v0 = 1; this._alim = 1;
    // связанные один раз — в кадре не создаются
    this._asinhTickFmt = (v) => {
      const val = Math.sinh(v * this._alim) * this._v0 / this.unitK;
      return fmtFixed(val, Math.abs(val) < 10 ? 2 : 0);
    };
  }

  /** Значение σ (Кл/м²) -> пиксель Y, с учётом режима шкалы. */
  _yOf(v) {
    return this.scaleMode === 'asinh'
      ? this.y.to(Math.asinh(v / this._v0) / this._alim)
      : this.y.to(v / this.unitK);
  }

  /** @param {Float64Array|number[]} r мм, @param {Float32Array} sigmaL,sigmaR Кл/м² */
  setData(r, sigmaL, sigmaR) {
    this.r = r; this.sL = sigmaL; this.sR = sigmaR;
    this.requestRender();
    return this;
  }
  setScaleMode(m) { this.scaleMode = m; this.requestRender(); return this; }

  draw(g) {
    if (!this.r || !this.sL) { drawEmpty(g, this, 'нет кадра'); return; }
    const n = Math.min(this.r.length, this.sL.length);
    let mx = 0;
    for (let i = 0; i < n; i++) {
      const a = Math.abs(this.sL[i]); if (a > mx) mx = a;
      const b = Math.abs(this.sR ? this.sR[i] : 0); if (b > mx) mx = b;
    }
    if (!(mx > 0)) mx = 1e-9;
    // сглаживание верхней границы, чтобы шкала не мерцала покадрово
    this._smooth = this._smooth > 0 ? Math.max(mx, this._smooth * 0.92) : mx;
    const lim = this._smooth * 1.1;

    this.x.setDomain(this.r[0], this.r[n - 1], false);
    const K = this.unitK;
    if (this.scaleMode === 'asinh') {
      this.y.setDomain(-1, 1, false);           // домен в asinh-координатах
      this._v0 = lim * 1e-3;
      this._alim = Math.asinh(lim / this._v0);
    } else {
      this.y.setDomain(-lim / K, lim / K, false);
    }

    this.xt.linear(this.x.min, this.x.max, Math.max(3, (this.px1 - this.px0) / 70));
    if (this.scaleMode === 'asinh') {
      this.yt._a = NaN;                                   // домен постоянный, подписи — нет
      this.yt.linear(-1, 1, 5, this._asinhTickFmt);
    } else {
      this.yt.linear(this.y.min, this.y.max, Math.max(3, (this.py0 - this.py1) / 34));
    }
    this._grid(g, this.x, this.y, this.xt, this.yt);

    const y0 = Math.round(this._yOf(0)) + 0.5;
    g.strokeStyle = THEME.axis; g.setLineDash(DASH_2_2);
    g.beginPath(); g.moveTo(this.px0, y0); g.lineTo(this.px1, y0); g.stroke();
    g.setLineDash(EMPTY_DASH);

    this._clipPlot(g);
    this._line(g, this.sL, n, THEME.sigmaL, 2, null);
    if (this.sR) this._line(g, this.sR, n, THEME.sigmaR, 2, DASH_6_3);
    g.restore();

    this._xLabels(g, this.x, this.xt, 'r, mm');
    this._yLabels(g, this.y, this.yt, 'σ, nC/cm²', 'left', THEME.textSecondary);

    // легенда + числа на оси
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left'; g.textBaseline = 'top';
    g.fillStyle = THEME.sigmaL;
    g.fillText(`σ_lo (z=0.5 мм)  ось: ${fmtFixed(this.sL[0] / K, 2)} nC/cm²`, this.px0 + 4, this.py1 + 3);
    if (this.sR) {
      g.fillStyle = THEME.sigmaR;
      g.fillText(`σ_hi (z=1.5 мм)  ось: ${fmtFixed(this.sR[0] / K, 2)} nC/cm²`, this.px0 + 4, this.py1 + 15);
    }
    if (this.scaleMode === 'asinh') {
      g.fillStyle = THEME.textMuted;
      g.textAlign = 'right';
      g.fillText('asinh-ось: линейна у нуля', this.px1 - 4, this.py1 + 3);
    }
  }

  _line(g, arr, n, color, width, dash) {
    g.beginPath();
    let open = false;
    for (let i = 0; i < n; i++) {
      const px = this.x.to(this.r[i]), py = this._yOf(arr[i]);
      if (!Number.isFinite(py)) { open = false; continue; }
      if (!open) { g.moveTo(px, py); open = true; } else g.lineTo(px, py);
    }
    g.strokeStyle = color; g.lineWidth = width;
    if (dash) g.setLineDash(dash); else g.setLineDash(EMPTY_DASH);
    g.stroke();
    g.setLineDash(EMPTY_DASH);
  }
}

// ═════════════════════════════════════════════════ AxialProfilePlot

/**
 * Реестр сортов. Поле `field` — имя в frames.bin (recorder.mjs defaultFields).
 *
 * ЧЕСТНО: recorder пишет ТОЛЬКО n_e, n_O3m, n_O3. Остальные шесть сортов в
 * контейнере отсутствуют — они помечаются в легенде как «нет в прогоне» и НЕ
 * рисуются. Дорисовывать их из воздуха плеер не имеет права.
 *
 * Цвета 1–6 — валидированная палитра UI_SPEC §3.2. Цвета для O₄⁺, O₃⁻, O₂(a)
 * добавлены сверх неё и через validate_palette.js НЕ проходили: для них
 * обязательны штрих и прямая подпись (они и так обязательны для всех).
 */
export const SPECIES = [
  { key: 'e', field: 'n_e', label: 'e⁻', color: '#3987e5', dash: null, width: 2.5 },
  { key: 'O2p', field: 'n_O2p', label: 'O₂⁺', color: '#e66767', dash: null, width: 2 },
  { key: 'O4p', field: 'n_O4p', label: 'O₄⁺', color: '#ff8b3d', dash: [5, 2], width: 2 },
  { key: 'Om', field: 'n_Om', label: 'O⁻', color: '#d55181', dash: [6, 3], width: 2 },
  { key: 'O2m', field: 'n_O2m', label: 'O₂⁻', color: '#c98500', dash: [2, 3], width: 2 },
  { key: 'O3m', field: 'n_O3m', label: 'O₃⁻', color: '#45b8c9', dash: [4, 2, 1, 2], width: 2 },
  { key: 'O', field: 'n_O', label: 'O', color: '#9085e9', dash: [8, 3, 2, 3], width: 2 },
  { key: 'O3', field: 'n_O3', label: 'O₃', color: '#199e70', dash: [1, 3], width: 2 },
  { key: 'O2a', field: 'n_O2a', label: 'O₂(a)', color: '#b0b8c4', dash: [3, 3], width: 2 },
];

/**
 * Лог-профили плотностей вдоль оси симметрии (или средних по r).
 *
 * Данные: setZ(zCenters) один раз + setSpecies(key, Float32Array|null) покадрово.
 * Удобная обёртка — updateFromPlayback() ниже.
 *
 * Легенда — кликабельная (canvas hit-test): клик = solo, alt+клик = скрыть,
 * повторный клик = вернуть всё.
 */
export class AxialProfilePlot extends BasePlot {
  constructor(canvas, opts = {}) {
    super(canvas, Object.assign({ pad: { l: 58, r: 58, t: 16, b: 50 } }, opts));
    this.z = null;
    this.data = new Map();      // key -> Float32Array | null
    this.present = new Map();   // key -> bool
    this.hidden = new Set();
    this.solo = null;
    this.minDecades = opts.minDecades || 6;
    this.floor = opts.floor || 1e10;
    // боксы легенды создаются ОДИН раз и только мутируются (ноль аллокаций в кадре)
    this._legendBoxes = SPECIES.map((sp) => ({ key: sp.key, x: 0, y: 0, w: 0, h: 12, label: sp.label, tw: 0 }));
    this._hitDirty = true;
    this._onClickBound = (ev) => this._onClick(ev);
    canvas.addEventListener('click', this._onClickBound);
    canvas.style.cursor = 'default';
  }

  setZ(z) { this.z = z; this.requestRender(); return this; }

  /** @param {string} key ключ из SPECIES; values=null => «нет в прогоне» */
  setSpecies(key, values) {
    const has = !!values;
    this.data.set(key, values || null);
    if (this.present.get(key) !== has) { this.present.set(key, has); this._hitDirty = true; }
    this.requestRender();
    return this;
  }

  _visible(sp) {
    if (!this.present.get(sp.key)) return false;
    if (this.solo) return this.solo === sp.key;
    return !this.hidden.has(sp.key);
  }

  _onClick(ev) {
    const r = this.canvas.getBoundingClientRect();
    const mx = ev.clientX - r.left, my = ev.clientY - r.top;
    for (let i = 0; i < this._legendBoxes.length; i++) {
      const b = this._legendBoxes[i];
      if (mx < b.x || mx > b.x + b.w || my < b.y || my > b.y + b.h) continue;
      if (!this.present.get(b.key)) return;
      if (ev.altKey) {
        if (this.hidden.has(b.key)) this.hidden.delete(b.key); else this.hidden.add(b.key);
      } else {
        this.solo = this.solo === b.key ? null : b.key;
      }
      this.render();
      return;
    }
  }

  draw(g) {
    if (!this.z) { drawEmpty(g, this, 'нет кадра'); return; }
    const nz = this.z.length;
    let vmax = 0;
    for (const sp of SPECIES) {
      if (!this._visible(sp)) continue;
      const a = this.data.get(sp.key);
      for (let j = 0; j < nz; j++) { const v = a[j]; if (v > vmax) vmax = v; }
    }
    if (!(vmax > 0)) vmax = 1e16;
    const hi = Math.pow(10, Math.ceil(Math.log10(vmax)));
    const lo = Math.max(this.floor, hi * Math.pow(10, -Math.max(this.minDecades, 6)));

    this.x.setDomain(this.z[0], this.z[nz - 1], false);
    this.y.setDomain(lo, hi, true);
    this.xt.linear(this.x.min, this.x.max, Math.max(3, (this.px1 - this.px0) / 70));
    this.yt.log(lo, hi, Math.max(3, (this.py0 - this.py1) / 30));
    this._grid(g, this.x, this.y, this.xt, this.yt);

    this._clipPlot(g);
    let clipped = 0;
    for (const sp of SPECIES) {
      if (!this._visible(sp)) continue;
      const a = this.data.get(sp.key);
      g.beginPath();
      let open = false;
      for (let j = 0; j < nz; j++) {
        const v = a[j];
        const py = this.y.to(v);
        if (!Number.isFinite(py)) { open = false; clipped++; continue; }
        const px = this.x.to(this.z[j]);
        if (!open) { g.moveTo(px, py); open = true; } else g.lineTo(px, py);
      }
      g.strokeStyle = sp.color;
      g.lineWidth = sp.width;
      if (sp.dash) g.setLineDash(sp.dash); else g.setLineDash(EMPTY_DASH);
      g.stroke();
      g.setLineDash(EMPTY_DASH);
    }
    g.restore();

    // прямые подписи у правого конца кривой (secondary encoding, UI_SPEC §3.2)
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left'; g.textBaseline = 'middle';
    for (const sp of SPECIES) {
      if (!this._visible(sp)) continue;
      const a = this.data.get(sp.key);
      let j = nz - 1;
      while (j > 0 && !(a[j] > 0)) j--;
      const py = this.y.to(a[j]);
      if (!Number.isFinite(py)) continue;
      g.fillStyle = sp.color;
      g.fillText(sp.label, this.px1 + 3, Math.max(this.py1 + 5, Math.min(this.py0 - 5, py)));
    }

    this._xLabels(g, this.x, this.xt, 'z, mm');
    this._yLabels(g, this.y, this.yt, 'n, m⁻³', 'left', THEME.textSecondary);
    if (clipped) {
      g.fillStyle = THEME.textMuted;
      g.font = THEME.fontMonoSmall;
      g.textAlign = 'right'; g.textBaseline = 'bottom';
      g.fillText(`${clipped} точек ≤ 0 или ниже ${fmtExp(lo)} — не показаны`, this.px1, this.py0 - 3);
    }
    this._drawLegend(g);
  }

  _drawLegend(g) {
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left'; g.textBaseline = 'top';
    const y = this.py0 + 24;
    let x = this.px0;
    if (this._hitDirty) {                    // подписи и их ширины — только при смене набора
      for (let i = 0; i < SPECIES.length; i++) {
        const sp = SPECIES[i], b = this._legendBoxes[i];
        b.label = sp.label + (this.present.get(sp.key) ? '' : ' ×');
        b.tw = g.measureText(b.label).width;
        b.w = 26 + b.tw;
      }
      this._hitDirty = false;
    }
    for (let i = 0; i < SPECIES.length; i++) {
      const sp = SPECIES[i], b = this._legendBoxes[i];
      const has = !!this.present.get(sp.key);
      const on = this._visible(sp);
      g.globalAlpha = has ? (on ? 1 : 0.35) : 0.28;
      g.strokeStyle = sp.color; g.lineWidth = sp.width;
      if (sp.dash) g.setLineDash(sp.dash); else g.setLineDash(EMPTY_DASH);
      g.beginPath(); g.moveTo(x, y + 6); g.lineTo(x + 13, y + 6); g.stroke();
      g.setLineDash(EMPTY_DASH);
      g.fillStyle = has ? sp.color : THEME.textMuted;
      g.fillText(b.label, x + 16, y);
      g.globalAlpha = 1;
      b.x = x; b.y = y;
      x += b.w;
    }
    g.fillStyle = THEME.textMuted;
    g.fillText('клик — solo · alt+клик — скрыть · × — сорт не записан в прогоне', this.px0, y + 13);
  }

  destroy() {
    this.canvas.removeEventListener('click', this._onClickBound);
    super.destroy();
  }
}

// ══════════════════════════════════════════════════════ PhotoPanel

/**
 * Каналы фотопанели. Три фотопроцесса + две реперные ударные реакции.
 * `source`: 'field' — объёмный интеграл поля кадра; 'series' — готовый ряд.
 */
export const CHANNELS = [
  { key: 'photoIon', label: 'фотоионизация', color: '#ffd166', dash: null, width: 2,
    source: 'field', field: 'photoIonRate', unit: 's⁻¹' },
  { key: 'photoEmit', label: 'фотоэмиссия (пов.)', color: '#4da3ff', dash: [5, 2], width: 2,
    source: 'series', series: ['photoEmitTotalL', 'photoEmitTotalR'], unit: 's⁻¹' },
  { key: 'photoDetach', label: 'фотоотлипание', color: '#9085e9', dash: [8, 3, 2, 3], width: 2,
    source: 'field', field: 'photoDetachRate', unit: 's⁻¹' },
  { key: 'impactIon', label: 'ударная ионизация', color: '#e66767', dash: null, width: 1.6,
    source: 'field', field: 'ionizRate', unit: 's⁻¹', reference: true },
  { key: 'collDetach', label: 'столкн. отлипание', color: '#199e70', dash: [1, 3], width: 1.6,
    source: 'none', unit: 's⁻¹', reference: true,
    note: 'не записано в контейнере (recorder не пишет detachRate)' },
];

/**
 * Объёмы выходных ячеек, м³. Сетка кадра — ТОЛЬКО газовый зазор
 * (recorder вырезает jGas0..jGas1), поэтому интеграл честно объёмный по газу.
 * @param {object} grid manifest.grid (rFaces/zFaces в мм)
 * @returns {Float64Array} длиной nrOut*nzOut, раскладка i*nzOut+j
 */
export function cellVolumes(grid) {
  const nr = grid.nrOut, nz = grid.nzOut;
  const rf = grid.rFaces, zf = grid.zFaces;
  const V = new Float64Array(nr * nz);
  for (let i = 0; i < nr; i++) {
    // мм² -> м²: 1e-6;  мм -> м: 1e-3
    const A = Math.PI * (rf[i + 1] * rf[i + 1] - rf[i] * rf[i]) * 1e-6;
    for (let j = 0; j < nz; j++) V[i * nz + j] = A * (zf[j + 1] - zf[j]) * 1e-3;
  }
  return V;
}

/** Σ v_k·V_k — интеграл поля по объёму газа. Поля прорежены площадным усреднением. */
export function integrateField(vals, vol) {
  let s = 0;
  const n = Math.min(vals.length, vol.length);
  for (let k = 0; k < n; k++) { const v = vals[k]; if (v > 0) s += v * vol[k]; }
  return s;
}

/**
 * Столбец поля вдоль оси (i=0) или средний по r.
 * @param {Float32Array} vals поле кадра (nrOut*nzOut)
 * @param {object} grid manifest.grid
 * @param {Float32Array} out длиной nzOut (переиспользуется — ноль аллокаций)
 * @param {'axis'|'mean'} mode
 */
export function axialColumn(vals, grid, out, mode = 'axis', vol = null) {
  const nr = grid.nrOut, nz = grid.nzOut;
  if (mode === 'axis') {
    for (let j = 0; j < nz; j++) out[j] = vals[j];
    return out;
  }
  for (let j = 0; j < nz; j++) {
    let s = 0, w = 0;
    for (let i = 0; i < nr; i++) {
      const ww = vol ? vol[i * nz + j] : 1;
      s += vals[i * nz + j] * ww; w += ww;
    }
    out[j] = w > 0 ? s / w : 0;
  }
  return out;
}

/**
 * Однократная сборка временных рядов фотоканалов из кадров.
 * Считается ОДИН раз при загрузке прогона (не в кадре анимации): объёмное
 * интегрирование 827 кадров × 3 поля — сотни мс, в rAF этому места нет.
 *
 *   const ph = await PhotoSeriesBuilder.build(playback, { onProgress });
 *   photoPanel.setData(ph);
 */
export class PhotoSeriesBuilder {
  /**
   * @param {import('./loader.mjs').Playback} pb
   * @param {{onProgress?:(done:number,total:number)=>void, chunk?:number}} opts
   * @returns {Promise<{t:Float64Array, chan:Object<string,Float32Array>, missing:string[]}>}
   */
  static async build(pb, opts = {}) {
    const grid = pb.grid;
    const vol = cellVolumes(grid);
    const nF = pb.frameCount;
    const t = Float64Array.from(pb.frameTimes);
    const chan = {};
    const missing = [];
    const available = new Set(pb.fields);
    for (const c of CHANNELS) {
      if (c.source === 'field') {
        if (!available.has(c.field)) { missing.push(c.key); continue; }
      } else if (c.source === 'series') {
        if (!pb.series || !c.series.every((s) => pb.series[s])) { missing.push(c.key); continue; }
      } else { missing.push(c.key); continue; }
      chan[c.key] = new Float32Array(nF);
    }
    const chunk = opts.chunk || 32;
    for (let i = 0; i < nF; i++) {
      for (const c of CHANNELS) {
        if (c.source !== 'field' || !chan[c.key]) continue;
        chan[c.key][i] = integrateField(await pb.getFrameAsync(i, c.field), vol);
      }
      if (i % chunk === chunk - 1) {
        if (opts.onProgress) opts.onProgress(i + 1, nF);
        await new Promise((res) => setTimeout(res, 0));   // не блокировать UI
      }
    }
    // поверхностные каналы: ближайший отсчёт series к времени кадра
    if (pb.series) {
      for (const c of CHANNELS) {
        if (c.source !== 'series' || !chan[c.key]) continue;
        const ts = pb.series.t, n = ts.length;
        for (let i = 0; i < nF; i++) {
          let k = lowerBound(ts, n, t[i]);
          if (k >= n) k = n - 1;
          let v = 0;
          for (const s of c.series) v += pb.series[s][k];
          chan[c.key][i] = v;
        }
      }
    }
    if (opts.onProgress) opts.onProgress(nF, nF);
    return { t, chan, missing };
  }
}

/**
 * Вклад фотоканалов во времени против ударных процессов, лог-ось.
 *
 * Смысл панели (PHOTO_PROCESSES.md): фотопроцессы почти всегда на 3–6 порядков
 * СЛАБЕЕ ударной ионизации по темпу — и при этом определяют, где и когда
 * зародится филамент. Панель показывает именно это: разрыв в темпе + отношение
 * фотоионизация/ударная ионизация в точке курсора.
 */
export class PhotoPanel extends BasePlot {
  constructor(canvas, opts = {}) {
    super(canvas, Object.assign({ pad: { l: 62, r: 14, t: 16, b: 50 } }, opts));
    this.data = null;
    this.hidden = new Set();
    this.decades = opts.decades || 12;
    this.xs = new Scale();
    this._legendBoxes = CHANNELS.map((c) => ({ key: c.key, x: 0, y: 0, w: 0, h: 12, label: c.label, tw: 0 }));
    this._legendDirty = true;
    this._onClickBound = (ev) => this._onClick(ev);
    canvas.addEventListener('click', this._onClickBound);
  }

  /** @param {{t:Float64Array, chan:object, missing:string[]}} d — из PhotoSeriesBuilder */
  setData(d) { this.data = d; this._legendDirty = true; this.requestRender(); return this; }
  setCursor(t) { this.cursorT = t; this.requestRender(); return this; }
  setWindow(t0, t1) { this.tw0 = t0; this.tw1 = t1; this.requestRender(); return this; }

  _onClick(ev) {
    const r = this.canvas.getBoundingClientRect();
    const mx = ev.clientX - r.left, my = ev.clientY - r.top;
    for (const b of this._legendBoxes) {
      if (mx < b.x || mx > b.x + b.w || my < b.y || my > b.y + b.h) continue;
      if (!this.data || !this.data.chan[b.key]) return;
      if (this.hidden.has(b.key)) this.hidden.delete(b.key); else this.hidden.add(b.key);
      this.render();
      return;
    }
  }

  draw(g) {
    const d = this.data;
    if (!d) { drawEmpty(g, this, 'фоторяды не построены'); return; }
    const n = d.t.length;
    const t0 = Number.isFinite(this.tw0) ? this.tw0 : d.t[0];
    const t1 = Number.isFinite(this.tw1) ? this.tw1 : d.t[n - 1];
    const i0 = Math.max(0, lowerBound(d.t, n, t0) - 1);
    const i1 = Math.min(n - 1, lowerBound(d.t, n, t1));

    let vmax = 0;
    for (const c of CHANNELS) {
      const a = d.chan[c.key];
      if (!a || this.hidden.has(c.key)) continue;
      for (let i = i0; i <= i1; i++) if (a[i] > vmax) vmax = a[i];
    }
    if (!(vmax > 0)) vmax = 1;
    const hi = Math.pow(10, Math.ceil(Math.log10(vmax)));
    const lo = hi * Math.pow(10, -this.decades);

    const ts = pickTimeScale(t1 - t0);
    this.x.setDomain(t0, t1, false);
    this.y.setDomain(lo, hi, true);
    this.xt.linear(t0 / ts.k, t1 / ts.k, Math.max(3, (this.px1 - this.px0) / 90));
    this.yt.log(lo, hi, Math.max(3, (this.py0 - this.py1) / 26));
    const xs = dispScale(this.xs, this.x, ts.k);
    this._grid(g, xs, this.y, this.xt, this.yt);

    this._clipPlot(g);
    let zeros = 0;
    for (const c of CHANNELS) {
      const a = d.chan[c.key];
      if (!a || this.hidden.has(c.key)) continue;
      g.beginPath();
      let open = false;
      for (let i = i0; i <= i1; i++) {
        const py = this.y.to(a[i]);
        if (!Number.isFinite(py)) { open = false; zeros++; continue; }
        const px = this.x.to(d.t[i]);
        if (!open) { g.moveTo(px, py); open = true; } else g.lineTo(px, py);
      }
      g.strokeStyle = c.color;
      g.lineWidth = c.width;
      if (c.dash) g.setLineDash(c.dash); else g.setLineDash(EMPTY_DASH);
      g.stroke();
      g.setLineDash(EMPTY_DASH);
    }
    if (Number.isFinite(this.cursorT) && this.cursorT >= t0 && this.cursorT <= t1) {
      const px = Math.round(this.x.to(this.cursorT)) + 0.5;
      g.strokeStyle = THEME.cursor; g.lineWidth = 1;
      g.beginPath(); g.moveTo(px, this.py1); g.lineTo(px, this.py0); g.stroke();
    }
    g.restore();

    this._xLabels(g, xs, this.xt, `t, ${ts.unit}`);
    this._yLabels(g, this.y, this.yt, 'темп, событий/с (по объёму газа)', 'left', THEME.textSecondary);

    // отношение в точке курсора — главный смысл панели
    if (Number.isFinite(this.cursorT) && d.chan.photoIon && d.chan.impactIon) {
      const k = Math.min(n - 1, lowerBound(d.t, n, this.cursorT));
      const pi = d.chan.photoIon[k], ii = d.chan.impactIon[k];
      const ratio = ii > 0 ? pi / ii : NaN;
      g.font = THEME.fontMonoSmall;
      g.textAlign = 'right'; g.textBaseline = 'top';
      g.fillStyle = THEME.textMuted;
      g.fillText(`S_ph/S_imp = ${Number.isFinite(ratio) ? fmtExp(ratio, '', 2) : '—'}`, this.px1 - 3, this.py1 + 3);
      g.fillText(`фотоионизация ${fmtExp(pi, 's⁻¹', 2)}`, this.px1 - 3, this.py1 + 15);
    }
    if (zeros) {
      g.font = THEME.fontMonoSmall;
      g.textAlign = 'right'; g.textBaseline = 'bottom';
      g.fillStyle = THEME.textMuted;
      g.fillText(`${zeros} нулей — на лог-оси разрыв, а не ноль`, this.px1 - 3, this.py0 - 3);
    }
    this._drawLegend(g);
  }

  _drawLegend(g) {
    g.font = THEME.fontMonoSmall;
    g.textAlign = 'left'; g.textBaseline = 'top';
    const y = this.py0 + 24;
    let x = this.px0;
    if (this._legendDirty) {
      for (let i = 0; i < CHANNELS.length; i++) {
        const c = CHANNELS[i], b = this._legendBoxes[i];
        b.label = c.label + (this.data && this.data.chan[c.key] ? '' : ' ×');
        b.tw = g.measureText(b.label).width;
        b.w = 26 + b.tw;
      }
      this._legendDirty = false;
    }
    for (let i = 0; i < CHANNELS.length; i++) {
      const c = CHANNELS[i], b = this._legendBoxes[i];
      const has = !!(this.data && this.data.chan[c.key]);
      const on = has && !this.hidden.has(c.key);
      g.globalAlpha = has ? (on ? 1 : 0.35) : 0.28;
      g.strokeStyle = c.color; g.lineWidth = c.width;
      if (c.dash) g.setLineDash(c.dash); else g.setLineDash(EMPTY_DASH);
      g.beginPath(); g.moveTo(x, y + 6); g.lineTo(x + 13, y + 6); g.stroke();
      g.setLineDash(EMPTY_DASH);
      g.fillStyle = has ? c.color : THEME.textMuted;
      g.fillText(b.label, x + 16, y);
      g.globalAlpha = 1;
      b.x = x; b.y = y;
      x += b.w;
    }
    g.fillStyle = THEME.textMuted;
    g.fillText('× — канал не записан в контейнере · клик — скрыть/показать', this.px0, y + 13);
  }

  destroy() {
    this.canvas.removeEventListener('click', this._onClickBound);
    super.destroy();
  }
}

// ──────────────────────────────────────────── удобные обёртки

/**
 * Заполнить AxialProfilePlot из кадра плейбэка. Буферы столбцов создаются
 * один раз и переиспользуются (ноль аллокаций в кадре).
 * @param {AxialProfilePlot} plot
 * @param {import('./loader.mjs').Playback} pb
 * @param {number} frameIndex
 * @param {{mode:'axis'|'mean', vol?:Float64Array, buffers?:Map}} o
 */
export function updateAxialFromPlayback(plot, pb, frameIndex, o = {}) {
  const grid = pb.grid;
  if (!plot.z) plot.setZ(grid.zCenters);
  if (!plot._colBufs) plot._colBufs = new Map();
  if (!plot._avail || plot._availOf !== pb) {   // Set строится один раз на прогон
    plot._avail = new Set(pb.fields);
    plot._availOf = pb;
  }
  const avail = plot._avail;
  for (const sp of SPECIES) {
    if (!avail.has(sp.field)) {
      if (plot.present.get(sp.key) !== false) plot.setSpecies(sp.key, null);
      continue;
    }
    let buf = plot._colBufs.get(sp.key);
    if (!buf || buf.length !== grid.nzOut) {
      buf = new Float32Array(grid.nzOut);
      plot._colBufs.set(sp.key, buf);
    }
    axialColumn(pb.getFrame(frameIndex, sp.field), grid, buf, o.mode || 'axis', o.vol);
    plot.setSpecies(sp.key, buf);
  }
  return plot;
}

export default {
  WaveformPlot, LissajousPlot, RadialProfilePlot, AxialProfilePlot, PhotoPanel,
  PhotoSeriesBuilder, SPECIES, CHANNELS, cellVolumes, integrateField, axialColumn,
  updateAxialFromPlayback, fitLissajous,
};
