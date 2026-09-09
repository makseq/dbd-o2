// scene.js — главная сцена: осесимметричный срез (r,z), зеркально отражённый
// относительно оси r = 0, с явной геометрией «металл / диэлектрик / газ»,
// накладками sigma(r) на обеих поверхностях и цветовой шкалой.
//
// Геометрия (мм, из manifest.grid):
//   z: 0 .. 0.5 диэлектрик | 0.5 .. 1.5 ГАЗ | 1.5 .. 2.0 диэлектрик
//   r: 0 .. R (данные) -> на экране -R .. +R (зеркало)
// Пропорции соблюдены строго: 1 px = const мм по обеим осям.
//
// Поля кадра записаны ТОЛЬКО по газовой части сетки (nrOut x nzOut, i-major).
// Сетка НЕРАВНОМЕРНАЯ (сгущение к оси и к стенкам), поэтому картинка не
// «растягивается» тупым imageSmoothing по индексам: строится билинейная
// выборка по ФИЗИЧЕСКИМ координатам через таблицы дробных индексов, а уже
// её результат масштабируется сглаживанием.

import { getLUT, isDiverging } from './colormaps.js';
import { t as tr, fieldLabel, fieldUnit } from './i18n.js';


const OFF_W = 288;          // разрешение внутреннего буфера газового окна (z)
const OFF_H = 288;          // (r, зеркально)  — газовое окно физически квадратное

const PAD = { l: 58, r: 100, t: 50, b: 56 };

const CSS = {
  bg: '#0a0d12',
  panel: '#10141a',
  metal: '#8d97a6',
  metalDark: '#5a626e',
  diel: '#232c3a',
  dielLine: '#3a475a',
  border: '#2a3442',
  text: '#e8eaed',
  dim: '#8b93a1',
  axis: '#39424f',
};

function fmtSci(v, digits = 2) {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) {
    const s = a >= 100 ? v.toFixed(0) : a >= 1 ? v.toFixed(2) : v.toFixed(3);
    return s.replace(/\.?0+$/, '');
  }
  const e = Math.floor(Math.log10(a));
  const m = v / Math.pow(10, e);
  return `${m.toFixed(digits)}e${e}`;
}

function makeHatch(color, bg, step = 7) {
  const c = document.createElement('canvas');
  c.width = step; c.height = step;
  const g = c.getContext('2d');
  if (bg) { g.fillStyle = bg; g.fillRect(0, 0, step, step); }
  g.strokeStyle = color; g.lineWidth = 1;
  g.beginPath();
  g.moveTo(-step, step); g.lineTo(step, -step);
  g.moveTo(0, 2 * step); g.lineTo(2 * step, 0);
  g.stroke();
  return c;
}

/** Дробный индекс x в массиве возрастающих центров (с зажимом на концах). */
function fracIndex(centers, x) {
  const n = centers.length;
  if (n === 1) return 0;
  if (x <= centers[0]) return 0;
  if (x >= centers[n - 1]) return n - 1;
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (centers[m] <= x) lo = m; else hi = m; }
  return lo + (x - centers[lo]) / (centers[hi] - centers[lo]);
}

export class Scene {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.off = document.createElement('canvas');
    this.off.width = OFF_W; this.off.height = OFF_H;
    this.offCtx = this.off.getContext('2d');
    this.img = this.offCtx.createImageData(OFF_W, OFF_H);
    this.acc = new Float32Array(OFF_W * OFF_H);      // послесвечение
    this.accActive = false;
    this.pb = null;
    this.frame = 0;
    this.field = 'n_e';
    this.opts = {
      scaleMode: 'frame',      // 'frame' | 'global'
      logScale: true,
      cmap: 'inferno',
      afterglow: false,
      decay: 0.86,
      showSigma: true,
      showProbe: true,
    };
    this.title = '';
    this.subtitle = '';
    this.pointer = null;
    this._imgKey = '';
    this._hatchDiel = makeHatch('rgba(120,140,170,0.16)', null, 8);
    this._hatchMetal = makeHatch('rgba(10,14,20,0.45)', null, 6);
    this.error = null;
  }

  // ── данные ───────────────────────────────────────────────────────────────
  setPlayback(pb) {
    this.pb = pb;
    this.acc.fill(0);
    this._imgKey = '';
    this._lut = null;
    if (!pb) return;
    const g = pb.grid;
    this.nr = g.nrOut; this.nz = g.nzOut;
    this.rCenters = g.rCenters; this.zCenters = g.zCenters;
    this.R = g.rFaces[g.rFaces.length - 1];
    this.zGas0 = g.zFaces[0]; this.zGas1 = g.zFaces[g.zFaces.length - 1];
    const diel = g.dielectrics || [];
    this.Lz = diel.length ? diel[diel.length - 1].z1 : this.zGas1;
    this.dielectrics = diel;
    this.norm = new Float32Array(this.nr * this.nz);
    this._buildMaps();
    this._globalRanges = this._computeGlobalRanges();
    this._sigmaMax = this._computeSigmaMax();
  }

  _buildMaps() {
    // столбцы буфера -> дробный индекс по z; строки -> по |r| (зеркало)
    const zi = new Float32Array(OFF_W), ri = new Float32Array(OFF_H);
    for (let px = 0; px < OFF_W; px++) {
      const z = this.zGas0 + ((px + 0.5) / OFF_W) * (this.zGas1 - this.zGas0);
      zi[px] = fracIndex(this.zCenters, z);
    }
    for (let py = 0; py < OFF_H; py++) {
      const r = -this.R + ((py + 0.5) / OFF_H) * (2 * this.R);
      ri[py] = fracIndex(this.rCenters, Math.abs(r));
    }
    this.zIdx = zi; this.rIdx = ri;
  }

  _computeGlobalRanges() {
    const out = {};
    const frames = this.pb.manifest.frames;
    for (const f of this.pb.manifest.fields) {
      let mx = 0, mn = 0;
      for (let i = 0; i < frames.length; i++) {
        const c = frames[i].fields[f.name];
        if (!c) continue;
        if (Number.isFinite(c.max) && c.max > mx) mx = c.max;
        if (Number.isFinite(c.min) && c.min < mn) mn = c.min;
      }
      out[f.name] = { max: mx, min: mn };
    }
    return out;
  }

  _computeSigmaMax() {
    let m = 0;
    try {
      const n = this.pb.frameCount;
      const step = n > 600 ? Math.ceil(n / 600) : 1;
      for (let i = 0; i < n; i += step) {
        const s = this.pb.getSurface(i);
        for (let k = 0; k < s.sigmaL.length; k++) {
          const a = Math.abs(s.sigmaL[k]); if (a > m) m = a;
          const b = Math.abs(s.sigmaR[k]); if (b > m) m = b;
        }
      }
    } catch (e) { /* режим range: посчитаем по кадру */ }
    return m > 0 ? m : 1e-12;
  }

  setField(name) { if (name !== this.field) { this.field = name; this.acc.fill(0); this._imgKey = ''; } }
  setFrameIndex(i) { this.frame = i; }
  setOptions(o) {
    const prevGlow = this.opts.afterglow;
    Object.assign(this.opts, o);
    if (!prevGlow && this.opts.afterglow) this.acc.fill(0);
    this._imgKey = '';
  }
  resetAfterglow() { this.acc.fill(0); }
  setPointer(p) { this.pointer = p; }

  /** Принудительный верх шкалы (синхронизация двух сцен в режиме сравнения). */
  setForceMax(v) {
    const nv = v > 0 ? v : 0;
    if (nv !== (this._forceMax || 0)) { this._forceMax = nv; this._imgKey = ''; }
  }
  globalMax(field) {
    const g = this._globalRanges && this._globalRanges[field];
    return g ? Math.max(Math.abs(g.max), Math.abs(g.min)) : 0;
  }

  fieldSpec(name) { return this.pb ? this.pb.fieldByName.get(name || this.field) : null; }

  // ── нормировка значений поля в t ∈ [0,1] ─────────────────────────────────
  _normalize() {
    const pb = this.pb, name = this.field;
    const spec = pb.fieldByName.get(name);
    const codec = pb.getCodec(this.frame, name) || {};
    const g = this._globalRanges[name] || { max: 0, min: 0 };
    const useGlobal = this.opts.scaleMode === 'global';
    const dec = pb.getFrame(this.frame, name);
    const out = this.norm;
    const n = out.length;

    const isSigned = spec.map === 'asinh';
    let range;

    const forced = this._forceMax || 0;

    if (isSigned) {
      const amax = forced > 0 ? forced : (useGlobal
        ? Math.max(Math.abs(g.max), Math.abs(g.min))
        : Math.max(Math.abs(codec.max || 0), Math.abs(codec.min || 0)));
      // v0 кодека (1e-4 от максимума) для показа слишком мал: asinh с ним
      // вырождается в знак и вся картина заливается двумя цветами.
      const v0 = Math.max((codec.v0 && codec.v0 > 0) ? codec.v0 : amax * 1e-4, amax * 3e-3, 1e-30);
      const den = Math.asinh(amax / v0) || 1;
      for (let k = 0; k < n; k++) {
        const v = dec[k];
        out[k] = amax > 0 ? 0.5 + 0.5 * Math.asinh(v / v0) / den : 0.5;
      }
      range = { kind: 'sym', vmax: amax, v0, unit: fieldUnit(spec) };
    } else {
      const vmax = forced > 0 ? forced : (useGlobal ? g.max : (codec.max || 0));
      const decades = spec.decades || 6;
      if (this.opts.logScale) {
        const vmin = vmax > 0 ? vmax * Math.pow(10, -decades) : 0;
        if (!(vmax > 0)) { out.fill(0); range = { kind: 'log', vmin: 0, vmax: 0, unit: fieldUnit(spec) }; }
        else {
          const lmin = Math.log10(vmin), inv = 1 / (Math.log10(vmax) - lmin);
          for (let k = 0; k < n; k++) {
            const v = dec[k];
            out[k] = v > 0 ? Math.min(1, Math.max(0, (Math.log10(v) - lmin) * inv)) : 0;
          }
          range = { kind: 'log', vmin, vmax, unit: fieldUnit(spec) };
        }
      } else {
        if (!(vmax > 0)) { out.fill(0); range = { kind: 'lin', vmin: 0, vmax: 0, unit: fieldUnit(spec) }; }
        else {
          const inv = 1 / vmax;
          for (let k = 0; k < n; k++) out[k] = Math.min(1, Math.max(0, dec[k] * inv));
          range = { kind: 'lin', vmin: 0, vmax, unit: fieldUnit(spec) };
        }
      }
    }
    this.range = range;
    this.decoded = dec;
    return range;
  }

  _rasterize() {
    const key = `${this.frame}|${this.field}|${this.opts.scaleMode}|${this.opts.logScale}|${this.opts.cmap}`;
    const glow = this.opts.afterglow;
    if (key === this._imgKey && !glow) return;
    this._imgKey = key;

    this._normalize();
    const lut = getLUT(this.opts.cmap);
    const nz = this.nz;
    const norm = this.norm, data = this.img.data;
    const zi = this.zIdx, ri = this.rIdx;
    const acc = this.acc, decay = this.opts.decay;
    const diverg = isDiverging(this.opts.cmap);

    for (let py = 0; py < OFF_H; py++) {
      const fi = ri[py];
      const i0 = fi | 0, i1 = Math.min(this.nr - 1, i0 + 1), wi = fi - i0;
      const rowA = i0 * nz, rowB = i1 * nz;
      let o = py * OFF_W * 4;
      const accRow = py * OFF_W;
      for (let px = 0; px < OFF_W; px++) {
        const fj = zi[px];
        const j0 = fj | 0, j1 = Math.min(nz - 1, j0 + 1), wj = fj - j0;
        const a = norm[rowA + j0] + (norm[rowA + j1] - norm[rowA + j0]) * wj;
        const b = norm[rowB + j0] + (norm[rowB + j1] - norm[rowB + j0]) * wj;
        let t = a + (b - a) * wi;
        if (glow && !diverg) {
          const prev = acc[accRow + px] * decay;
          t = t > prev ? t : prev;
          acc[accRow + px] = t;
        }
        const q = (t <= 0 ? 0 : t >= 1 ? 255 : (t * 255) | 0) << 2;
        data[o++] = lut[q]; data[o++] = lut[q + 1]; data[o++] = lut[q + 2]; data[o++] = 255;
      }
    }
    this.offCtx.putImageData(this.img, 0, 0);
  }

  // ── компоновка ───────────────────────────────────────────────────────────
  _layout(w, h) {
    const availW = Math.max(40, w - PAD.l - PAD.r);
    const availH = Math.max(30, h - PAD.t - PAD.b);
    const aspect = this.Lz / (2 * this.R);
    let W = availW, H = W / aspect;
    if (H > availH) { H = availH; W = H * aspect; }
    const x0 = PAD.l + (availW - W) / 2;
    // Свободную высоту (типично в режиме сравнения) отдаём осевому профилю,
    // сцену прижимаем к верху; иначе центрируем.
    const spare = availH - H;
    let prof = null, y0;
    if (spare > 110) {
      y0 = PAD.t;
      prof = { y: PAD.t + H + 48, h: Math.min(210, spare - 58) };
    } else {
      y0 = PAD.t + spare / 2;
    }
    return { x0, y0, W, H, aspect, prof };
  }

  /** Профиль поля вдоль оси r = 0 (i = 0), выровненный по оси z сцены. */
  _drawProfile(ctx, L) {
    const p = L.prof;
    if (!p || !this.decoded) return;
    const xa = this.X(L, this.zGas0), xb = this.X(L, this.zGas1);
    const y0 = p.y, h = p.h;
    ctx.fillStyle = '#0b0f15';
    ctx.fillRect(xa, y0, xb - xa, h);
    ctx.strokeStyle = CSS.border; ctx.lineWidth = 1;
    ctx.strokeRect(xa + 0.5, y0 + 0.5, xb - xa - 1, h - 1);

    const nz = this.nz, dec = this.decoded, norm = this.norm;
    const lut = getLUT(this.opts.cmap);

    // сетка по вертикали (в единицах нормировки)
    for (let k = 1; k < 4; k++) {
      const y = y0 + (k / 4) * h;
      ctx.strokeStyle = 'rgba(255,255,255,0.05)';
      ctx.beginPath(); ctx.moveTo(xa, y + 0.5); ctx.lineTo(xb, y + 0.5); ctx.stroke();
    }

    // заливка под кривой цветом палитры
    ctx.beginPath();
    ctx.moveTo(xa, y0 + h);
    let jMax = 0, vMax = -Infinity;
    for (let j = 0; j < nz; j++) {
      const x = this.X(L, this.zCenters[j]);
      const t = Math.max(0, Math.min(1, norm[j]));
      ctx.lineTo(x, y0 + h - t * (h - 4));
      if (dec[j] > vMax) { vMax = dec[j]; jMax = j; }
    }
    ctx.lineTo(xb, y0 + h);
    ctx.closePath();
    const q = 200 * 4;
    ctx.fillStyle = `rgba(${lut[q]},${lut[q + 1]},${lut[q + 2]},0.18)`;
    ctx.fill();

    ctx.beginPath();
    for (let j = 0; j < nz; j++) {
      const x = this.X(L, this.zCenters[j]);
      const t = Math.max(0, Math.min(1, norm[j]));
      const y = y0 + h - t * (h - 4);
      if (j === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.strokeStyle = `rgb(${lut[q]},${lut[q + 1]},${lut[q + 2]})`;
    ctx.lineWidth = 1.6;
    ctx.stroke();

    // радиальная полуширина в сечении максимума
    let rHalf = NaN;
    if (vMax > 0) {
      const half = vMax / 2;
      let prevR = 0, prevV = dec[jMax];
      for (let i = 1; i < this.nr; i++) {
        const v = dec[i * nz + jMax];
        if (v <= half) {
          const rr = this.rCenters[i - 1] + (this.rCenters[i] - this.rCenters[i - 1])
            * ((prevV - half) / Math.max(1e-300, prevV - v));
          rHalf = rr; break;
        }
        prevR = this.rCenters[i]; prevV = v;
      }
      if (!Number.isFinite(rHalf)) rHalf = this.R;
    }

    ctx.font = '9.5px "JetBrains Mono", ui-monospace, monospace';
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    const rightTxt = tr('scene.maxAt', { v: fmtSci(vMax), z: this.zCenters[jMax].toFixed(3) });
    const leftTxt = tr('scene.profileOnAxis');
    ctx.fillStyle = '#8b93a1';
    if (ctx.measureText(leftTxt).width + ctx.measureText(rightTxt).width + 24 < xb - xa) {
      ctx.fillText(leftTxt, xa + 5, y0 + 4);
    }
    ctx.textAlign = 'right';
    ctx.fillStyle = '#c8ced8';
    ctx.fillText(rightTxt, xb - 5, y0 + 4);
    if (Number.isFinite(rHalf)) {
      ctx.fillStyle = '#fab219';
      ctx.fillText(tr('scene.rHalf', { v: (rHalf * 1000).toFixed(1) }), xb - 5, y0 + 17);
    }
    // отметка максимума
    const xm = this.X(L, this.zCenters[jMax]);
    ctx.strokeStyle = 'rgba(250,178,25,0.6)';
    ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(xm + 0.5, y0); ctx.lineTo(xm + 0.5, y0 + h); ctx.stroke();
    ctx.setLineDash([]);
  }

  X(L, z) { return L.x0 + (z / this.Lz) * L.W; }
  Y(L, r) { return L.y0 + ((r + this.R) / (2 * this.R)) * L.H; }

  /** Физические координаты (z,r) в мм по позиции курсора, либо null. */
  probeAt(cssX, cssY) {
    if (!this.pb || !this._L) return null;
    const L = this._L;
    const z = ((cssX - L.x0) / L.W) * this.Lz;
    const r = ((cssY - L.y0) / L.H) * (2 * this.R) - this.R;
    if (z < 0 || z > this.Lz || Math.abs(r) > this.R) return null;
    return { z, r };
  }

  // ── отрисовка ────────────────────────────────────────────────────────────
  render() {
    const cv = this.canvas, ctx = this.ctx;
    const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
    const w = cv.clientWidth, h = cv.clientHeight;
    if (!w || !h) return;
    const pw = Math.round(w * dpr), ph = Math.round(h * dpr);
    if (cv.width !== pw || cv.height !== ph) { cv.width = pw; cv.height = ph; }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = CSS.bg;
    ctx.fillRect(0, 0, w, h);

    if (!this.pb) {
      ctx.fillStyle = CSS.dim;
      ctx.font = '13px Inter, system-ui, sans-serif';
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(this.error || tr('scene.noFrameData'), w / 2, h / 2);
      return;
    }

    const L = this._layout(w, h);
    this._L = L;
    this._rasterize();

    this._drawSolids(ctx, L);
    this._drawGas(ctx, L);
    this._drawSigma(ctx, L);
    this._drawFrameAndAxes(ctx, L);
    this._drawProfile(ctx, L);
    this._drawColorbar(ctx, L, w, h);
    this._drawHUD(ctx, L, w, h);
    this._drawProbe(ctx, L);
  }

  _drawSolids(ctx, L) {
    const yT = L.y0, yB = L.y0 + L.H;
    // диэлектрики
    for (const d of this.dielectrics) {
      if (d.gas) continue;
      const xa = this.X(L, d.z0), xb = this.X(L, d.z1);
      ctx.fillStyle = CSS.diel;
      ctx.fillRect(xa, yT, xb - xa, L.H);
      const pat = ctx.createPattern(this._hatchDiel, 'repeat');
      ctx.fillStyle = pat;
      ctx.fillRect(xa, yT, xb - xa, L.H);
      ctx.strokeStyle = CSS.dielLine; ctx.lineWidth = 1;
      ctx.strokeRect(xa + 0.5, yT + 0.5, xb - xa - 1, L.H - 1);
      // подпись
      ctx.save();
      ctx.fillStyle = 'rgba(232,234,237,0.82)';
      ctx.font = '600 11px Inter, system-ui, sans-serif';
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      const cx = (xa + xb) / 2;
      if (xb - xa > 74) {
        ctx.fillText('Al₂O₃', cx, yT + L.H / 2 - 9);
        ctx.font = '10px Inter, system-ui, sans-serif';
        ctx.fillStyle = 'rgba(169,176,187,0.85)';
        ctx.fillText('εᵣ = 9', cx, yT + L.H / 2 + 6);
        ctx.fillText(tr('scene.dielThick', { v: (d.z1 - d.z0).toFixed(1) }), cx, yT + L.H / 2 + 20);
      }
      ctx.restore();
    }
    // электроды (плоскости z = 0 и z = Lz)
    const bar = 11;
    const drawMetal = (x, dir) => {
      const gx = ctx.createLinearGradient(x, 0, x + dir * bar, 0);
      gx.addColorStop(0, CSS.metal); gx.addColorStop(1, CSS.metalDark);
      ctx.fillStyle = gx;
      const xa = dir > 0 ? x : x - bar;
      ctx.fillRect(xa, yT - 4, bar, L.H + 8);
      const pat = ctx.createPattern(this._hatchMetal, 'repeat');
      ctx.fillStyle = pat; ctx.fillRect(xa, yT - 4, bar, L.H + 8);
      ctx.strokeStyle = '#aeb7c4'; ctx.lineWidth = 1;
      ctx.strokeRect(xa + 0.5, yT - 3.5, bar - 1, L.H + 7);
    };
    drawMetal(this.X(L, 0), -1);
    drawMetal(this.X(L, this.Lz), 1);

    ctx.font = '600 10px Inter, system-ui, sans-serif';
    ctx.fillStyle = '#c8ced8';
    ctx.textBaseline = 'alphabetic';
    ctx.textAlign = 'left';
    ctx.fillText(tr('scene.metalHV'), this.X(L, 0) - bar, yT - 10);
    ctx.textAlign = 'right';
    ctx.fillText(tr('scene.metalGnd'), this.X(L, this.Lz) + bar, yT - 10);
  }

  _drawGas(ctx, L) {
    const xa = this.X(L, this.zGas0), xb = this.X(L, this.zGas1);
    ctx.save();
    ctx.beginPath(); ctx.rect(xa, L.y0, xb - xa, L.H); ctx.clip();
    ctx.imageSmoothingEnabled = true;
    if ('imageSmoothingQuality' in ctx) ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(this.off, 0, 0, OFF_W, OFF_H, xa, L.y0, xb - xa, L.H);
    // ось r = 0
    ctx.strokeStyle = 'rgba(232,234,237,0.28)';
    ctx.setLineDash([5, 5]); ctx.lineWidth = 1;
    const y = this.Y(L, 0);
    ctx.beginPath(); ctx.moveTo(xa, y + 0.5); ctx.lineTo(xb, y + 0.5); ctx.stroke();
    ctx.setLineDash([]);
    ctx.restore();
    ctx.strokeStyle = 'rgba(90,110,140,0.75)'; ctx.lineWidth = 1;
    ctx.strokeRect(xa + 0.5, L.y0 + 0.5, xb - xa - 1, L.H - 1);
    ctx.fillStyle = 'rgba(169,176,187,0.9)';
    ctx.font = '10px Inter, system-ui, sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.fillText(tr('scene.gas'), (xa + xb) / 2, L.y0 - 10);
    ctx.save();
    ctx.fillStyle = 'rgba(232,234,237,0.5)';
    ctx.font = '9px Inter, system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(tr('scene.axis0'), xa + 5, this.Y(L, 0) - 4);
    ctx.restore();
  }

  _drawSigma(ctx, L) {
    if (!this.opts.showSigma) return;
    let s;
    try { s = this.pb.getSurface(this.frame); } catch (e) { return; }
    const lut = getLUT('sigma');
    const amax = this._sigmaMax;
    const wStrip = 8;
    const draw = (arr, xEdge, dir) => {
      const x = dir > 0 ? xEdge : xEdge - wStrip;
      const n = arr.length;
      const H = L.H;
      const steps = Math.max(48, Math.min(240, Math.round(H)));
      for (let k = 0; k < steps; k++) {
        const r = -this.R + ((k + 0.5) / steps) * 2 * this.R;
        const fi = fracIndex(this.rCenters, Math.abs(r));
        const i0 = fi | 0, i1 = Math.min(n - 1, i0 + 1), wi = fi - i0;
        const v = arr[i0] + (arr[i1] - arr[i0]) * wi;
        const t = 0.5 + 0.5 * Math.max(-1, Math.min(1, v / amax));
        const q = Math.max(0, Math.min(255, (t * 255) | 0)) * 4;
        ctx.fillStyle = `rgb(${lut[q]},${lut[q + 1]},${lut[q + 2]})`;
        const y0 = L.y0 + (k / steps) * H, y1 = L.y0 + ((k + 1) / steps) * H;
        ctx.fillRect(x, y0, wStrip, y1 - y0 + 0.6);
      }
      ctx.strokeStyle = 'rgba(200,210,225,0.55)'; ctx.lineWidth = 1;
      ctx.strokeRect(x + 0.5, L.y0 + 0.5, wStrip - 1, L.H - 1);
    };
    draw(s.sigmaL, this.X(L, this.zGas0), -1);
    draw(s.sigmaR, this.X(L, this.zGas1), 1);

    ctx.fillStyle = 'rgba(169,176,187,0.95)';
    ctx.font = '9px "JetBrains Mono", ui-monospace, monospace';
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.fillText('σ(r)', this.X(L, this.zGas0) - wStrip / 2, L.y0 + L.H + 5);
    ctx.fillText('σ(r)', this.X(L, this.zGas1) + wStrip / 2, L.y0 + L.H + 5);
    this._sigmaNow = s;
  }

  _drawFrameAndAxes(ctx, L) {
    ctx.strokeStyle = CSS.border; ctx.lineWidth = 1;
    ctx.strokeRect(L.x0 + 0.5, L.y0 + 0.5, L.W - 1, L.H - 1);

    ctx.fillStyle = CSS.dim;
    ctx.font = '10px "JetBrains Mono", ui-monospace, monospace';
    // ось z
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const z of [0, 0.5, 1.0, 1.5, 2.0]) {
      if (z > this.Lz + 1e-9) continue;
      const x = this.X(L, z);
      ctx.strokeStyle = CSS.axis;
      ctx.beginPath(); ctx.moveTo(x, L.y0 + L.H); ctx.lineTo(x, L.y0 + L.H + 4); ctx.stroke();
      ctx.fillText(z.toFixed(1), x, L.y0 + L.H + 17);
    }
    ctx.textAlign = 'center';
    ctx.fillText(tr('scene.axisZ'), L.x0 + L.W / 2, L.y0 + L.H + 30);
    // ось r
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    const rt = [-this.R, -this.R / 2, 0, this.R / 2, this.R];
    for (const r of rt) {
      const y = this.Y(L, r);
      ctx.strokeStyle = CSS.axis;
      ctx.beginPath(); ctx.moveTo(L.x0 - 16, y); ctx.lineTo(L.x0 - 13, y); ctx.stroke();
      ctx.fillStyle = Math.abs(r) < 1e-9 ? '#c8ced8' : CSS.dim;
      ctx.fillText(r.toFixed(2), L.x0 - 20, y);
    }
    ctx.fillStyle = CSS.dim;
    ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
    ctx.fillText(tr('scene.axisR'), 6, L.y0 - 10);
  }

  _drawColorbar(ctx, L, w, h) {
    const range = this.range || { kind: 'lin', vmin: 0, vmax: 0, unit: '' };
    const bx = Math.min(w - 74, L.x0 + L.W + 34);
    const by = L.y0, bh = L.H, bw = 15;
    const lut = getLUT(this.opts.cmap);
    for (let k = 0; k < bh; k++) {
      const t = 1 - k / (bh - 1 || 1);
      const q = Math.max(0, Math.min(255, (t * 255) | 0)) * 4;
      ctx.fillStyle = `rgb(${lut[q]},${lut[q + 1]},${lut[q + 2]})`;
      ctx.fillRect(bx, by + k, bw, 1.2);
    }
    ctx.strokeStyle = CSS.border; ctx.lineWidth = 1;
    ctx.strokeRect(bx + 0.5, by + 0.5, bw - 1, bh - 1);

    ctx.fillStyle = CSS.dim;
    ctx.font = '9.5px "JetBrains Mono", ui-monospace, monospace';
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    const tick = (t, label, strong) => {
      const y = by + (1 - t) * bh;
      ctx.strokeStyle = 'rgba(232,234,237,0.5)';
      ctx.beginPath(); ctx.moveTo(bx + bw, y); ctx.lineTo(bx + bw + 4, y); ctx.stroke();
      ctx.fillStyle = strong ? '#c8ced8' : CSS.dim;
      ctx.fillText(label, bx + bw + 7, y);
    };
    if (range.kind === 'sym') {
      tick(1, `+${fmtSci(range.vmax)}`, true);
      tick(0.5, '0', true);
      tick(0, `−${fmtSci(range.vmax)}`, true);
      tick(0.75, `+${fmtSci(range.v0 * Math.sinh(Math.asinh(range.vmax / range.v0) * 0.5))}`);
      tick(0.25, `−${fmtSci(range.v0 * Math.sinh(Math.asinh(range.vmax / range.v0) * 0.5))}`);
    } else if (range.kind === 'log' && range.vmax > 0) {
      const l0 = Math.log10(range.vmin), l1 = Math.log10(range.vmax);
      const span = l1 - l0;
      const stepDec = span > 8 ? 2 : 1;
      const start = Math.ceil(l0 / stepDec) * stepDec;
      for (let e = start; e <= l1 + 1e-9; e += stepDec) {
        tick((e - l0) / span, `1e${e}`, false);
      }
      tick(1, fmtSci(range.vmax), true);
      ctx.fillStyle = CSS.dim;
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText('≤0', bx + bw + 7, by + bh + 8);
    } else {
      for (let k = 0; k <= 4; k++) tick(k / 4, fmtSci(range.vmax * k / 4), k === 4);
    }
    // единицы — ПОД шкалой, чтобы не сталкиваться с подписью электрода
    ctx.save();
    ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#a9b0bb';
    ctx.font = '10px Inter, system-ui, sans-serif';
    ctx.fillText(range.unit || '', bx - 2, by + bh + 22);
    ctx.fillStyle = '#59616e';
    ctx.font = '9px Inter, system-ui, sans-serif';
    const modeTxt = (range.kind === 'log' ? tr('scale.log') : range.kind === 'sym' ? 'asinh ±' : tr('scale.lin'))
      + (this._forceMax ? tr('scene.scaleSync') : '');
    ctx.fillText(modeTxt, bx - 2, by + bh + 34);
    ctx.restore();
  }

  _drawHUD(ctx, L, w, h) {
    const spec = this.fieldSpec();
    const fr = this.pb.manifest.frames[this.frame];
    ctx.font = '600 12px Inter, system-ui, sans-serif';
    ctx.fillStyle = CSS.text;
    ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
    const specLabel = fieldLabel(spec);
    const head = this.title ? `${this.title} · ${specLabel}` : specLabel;
    ctx.fillText(head, 8, 16);
    ctx.font = '10px "JetBrains Mono", ui-monospace, monospace';
    ctx.fillStyle = CSS.dim;
    ctx.textAlign = 'right';
    const t = fr ? fr.t : 0;
    ctx.fillText(tr('scene.hudTime', { t: (t * 1e6).toFixed(4), i: this.frame + 1, n: this.pb.frameCount }), w - 8, 16);
    if (this.opts.afterglow) {
      ctx.textAlign = 'left';
      ctx.fillStyle = 'rgba(250,178,25,0.9)';
      ctx.fillText(tr('scene.afterglow'), 8, 30);
    }
    if (this.subtitle) {
      ctx.textAlign = 'right'; ctx.fillStyle = CSS.dim;
      ctx.fillText(this.subtitle, w - 8, 30);
    }
  }

  _drawProbe(ctx, L) {
    if (!this.opts.showProbe || !this.pointer) return;
    const p = this.probeAt(this.pointer.x, this.pointer.y);
    if (!p) return;
    const x = this.X(L, p.z), y = this.Y(L, p.r);
    ctx.save();
    ctx.strokeStyle = 'rgba(232,234,237,0.35)'; ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(L.x0, y + 0.5); ctx.lineTo(L.x0 + L.W, y + 0.5);
    ctx.moveTo(x + 0.5, L.y0); ctx.lineTo(x + 0.5, L.y0 + L.H);
    ctx.stroke();

    let txt = tr('scene.probePos', { r: p.r.toFixed(3), z: p.z.toFixed(3) });
    let val = null;
    if (p.z >= this.zGas0 && p.z <= this.zGas1 && this.decoded) {
      const i = Math.round(fracIndex(this.rCenters, Math.abs(p.r)));
      const j = Math.round(fracIndex(this.zCenters, p.z));
      val = this.decoded[i * this.nz + j];
    } else if (this._sigmaNow) {
      const arr = p.z < this.zGas0 ? this._sigmaNow.sigmaL : this._sigmaNow.sigmaR;
      const i = Math.round(fracIndex(this.rCenters, Math.abs(p.r)));
      val = arr[Math.min(arr.length - 1, i)];
      txt += tr('scene.probeSigma', { v: fmtSci(val) });
      val = null;
    }
    if (val !== null) txt += `   ${fmtSci(val, 3)} ${fieldUnit(this.fieldSpec())}`;
    ctx.font = '10px "JetBrains Mono", ui-monospace, monospace';
    const tw = ctx.measureText(txt).width + 12;
    let bxx = x + 10, byy = y - 24;
    if (bxx + tw > L.x0 + L.W) bxx = x - 10 - tw;
    if (byy < L.y0 + 2) byy = y + 10;
    ctx.fillStyle = 'rgba(10,13,18,0.88)';
    ctx.strokeStyle = 'rgba(90,110,140,0.6)';
    ctx.beginPath(); ctx.rect(bxx, byy, tw, 18); ctx.fill(); ctx.stroke();
    ctx.fillStyle = CSS.text; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    ctx.fillText(txt, bxx + 6, byy + 9);
    ctx.restore();
  }
}

export default Scene;
