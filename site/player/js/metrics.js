// metrics.js — форматирование инженерных величин + панель метрик (DOM).
//
// Зона ответственности:
//   * fmtSI / fmtEng / fmtTime / fmtPct — единственное место в плеере, где число
//     превращается в строку с единицей. Никакого toFixed по месту в графиках.
//   * MetricsPanel — плитки метрик (UI_SPEC §2.5) + строка статуса солвера.
//
// Модуль не трогает css/app.css: свои стили инжектит один раз в <head> под
// префиксом .mp- (см. injectStyles). Классы уникальны и с app.css не пересекаются.
//
// Контракт: docs/PLOTS_API.md (в этом же каталоге — PLOTS_API.md).

import { t as tr } from './i18n.js';

// ───────────────────────────────────────────────────────── тема
// Значения обязаны совпадать с UI_SPEC §3.1. Здесь они продублированы числами,
// потому что canvas не умеет читать CSS-переменные без getComputedStyle на кадр
// (это аллокация строки каждый кадр — запрещено).
export const THEME = {
  bg: '#0a0d12',
  surface1: '#10141a',
  surface2: '#161b23',
  surface3: '#1c2430',
  border: '#242c38',
  borderStrong: '#333f4f',
  text: '#e8eaed',
  textSecondary: '#a9b0bb',
  textMuted: '#8b93a1',
  grid: '#1e2530',
  axis: '#39424f',
  accent: '#4da3ff',
  cursor: '#c8ced8',
  good: '#0ca30c',
  warning: '#fab219',
  serious: '#ec835a',
  critical: '#d03b3b',
  // вспомогательные серии (UI_SPEC §3.2, «Вспомогательные»)
  Uapp: '#e8eaed',
  Ugap: '#a9b0bb',
  Icond: '#ffd166',
  Itot: '#ffd166',
  Idisp: '#7f8794',
  sigmaL: '#3987e5',
  sigmaR: '#e66767',
  fontMono: '11px ui-monospace, "JetBrains Mono", SFMono-Regular, Menlo, monospace',
  fontMonoSmall: '10px ui-monospace, "JetBrains Mono", SFMono-Regular, Menlo, monospace',
  fontLabel: '11px Inter, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
};

// ───────────────────────────────────────────── форматирование чисел

const PREFIX = ['y', 'z', 'a', 'f', 'p', 'n', 'µ', 'm', '', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'];
const PREFIX_ZERO = 8; // индекс пустого префикса

/** Значащие цифры без экспоненты, если это возможно без вранья. */
function sigStr(v, sig) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  const d = Math.max(0, sig - 1 - Math.floor(Math.log10(a)));
  return v.toFixed(Math.min(20, d));
}

/**
 * Инженерный формат с приставкой СИ: 6.26e-15 -> "6.26 fF".
 * Ноль -> "0 <unit>"; не-число -> em-dash. Никогда не бросает.
 * @param {number} v
 * @param {string} unit  единица без приставки ('A', 'F', 'W', 'V', 'm^-3')
 * @param {number} sig   значащих цифр (по умолчанию 3)
 */
export function fmtSI(v, unit = '', sig = 3) {
  if (v == null || !Number.isFinite(v)) return '—';
  if (v === 0) return unit ? `0 ${unit}` : '0';
  const a = Math.abs(v);
  let e = Math.floor(Math.log10(a) / 3);
  e = Math.max(-8, Math.min(8, e));
  const idx = PREFIX_ZERO + e;
  const m = v / Math.pow(1000, e);
  const p = PREFIX[idx];
  const num = sigStr(m, sig);
  return unit || p ? `${num} ${p}${unit}` : num;
}

/**
 * Явная экспонента — для величин, у которых приставка бессмысленна
 * (плотности 1e21 м^-3, скорости 1e29 м^-3 с^-1).
 */
export function fmtExp(v, unit = '', sig = 3) {
  if (v == null || !Number.isFinite(v)) return '—';
  if (v === 0) return unit ? `0 ${unit}` : '0';
  const e = Math.floor(Math.log10(Math.abs(v)));
  const m = v / Math.pow(10, e);
  return `${sigStr(m, sig)}e${e >= 0 ? '+' : '−'}${String(Math.abs(e)).padStart(2, '0')}${unit ? ' ' + unit : ''}`;
}

/** Фиксированное число знаков после запятой, с разделителем тысяч тонким пробелом. */
export function fmtFixed(v, digits = 2, unit = '') {
  if (v == null || !Number.isFinite(v)) return '—';
  let s = v.toFixed(digits);
  const neg = s[0] === '-';
  if (neg) s = s.slice(1);
  const dot = s.indexOf('.');
  const intPart = dot < 0 ? s : s.slice(0, dot);
  const rest = dot < 0 ? '' : s.slice(dot);
  let grouped = '';
  for (let i = 0; i < intPart.length; i++) {
    if (i > 0 && (intPart.length - i) % 3 === 0) grouped += ' ';
    grouped += intPart[i];
  }
  return (neg ? '−' : '') + grouped + rest + (unit ? ' ' + unit : '');
}

/** Время в секундах -> удобная приставка, 4 знака мантиссы: 8.2984e-6 -> "8.2984 µs". */
export function fmtTime(t, sig = 5) {
  if (t == null || !Number.isFinite(t)) return '—';
  return fmtSI(t, 's', sig);
}

/** Относительное расхождение в процентах со знаком. */
export function fmtPct(x, digits = 1) {
  if (x == null || !Number.isFinite(x)) return '—';
  const s = x >= 0 ? '+' : '−';
  return `${s}${Math.abs(x).toFixed(digits)} %`;
}

/** Подпись для лог-тика: 1e21 -> "10²¹". */
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
export function fmtPow10(exp) {
  let s = '';
  const str = String(exp);
  for (let i = 0; i < str.length; i++) s += SUP[str[i]] ?? str[i];
  return '10' + s;
}

// ───────────────────────────────────────── статус по порогам

/**
 * Классификация значения по порогам. Возвращает 'ok'|'warning'|'serious'|'critical'.
 * Пороги задаются массивом [warn, serious, critical] в порядке возрастания «плохости»;
 * если dir === -1, плохо — когда МЕНЬШЕ порога.
 */
export function classify(v, thresholds, dir = 1) {
  if (v == null || !Number.isFinite(v) || !thresholds) return 'ok';
  const names = ['warning', 'serious', 'critical'];
  let out = 'ok';
  for (let i = 0; i < thresholds.length && i < 3; i++) {
    const th = thresholds[i];
    if (th == null) continue;
    if (dir > 0 ? v >= th : v <= th) out = names[i];
  }
  return out;
}

export const STATUS_COLOR = {
  ok: THEME.text, good: THEME.good, warning: THEME.warning,
  serious: THEME.serious, critical: THEME.critical,
};

// ───────────────────────────────────────────── аналитические ёмкости

/**
 * Аналитика ёмкостей ячейки (ERRATA §C V4, NUMERICS_2D §7.5).
 *   C_diel = eps0*eps_r*A/(d1+d2);  C_gas = eps0*A/d_gap;  C_cell = послед. соединение.
 * Всё честно из геометрии, без подгонки. Дополнительно даётся пересчёт на 1 см²,
 * чтобы можно было глазами сверить с каноническими 7.97 / 0.885 / 0.797 пФ.
 *
 * @param {{epsR:number, gapMM:number, dielMM:number, radiusMM:number}} g
 *        dielMM — толщина ОДНОГО барьера (их два).
 */
export function analyticCapacitances(g) {
  const EPS0 = 8.8541878128e-12;
  const epsR = g.epsR, dGap = g.gapMM * 1e-3, dDiel = 2 * g.dielMM * 1e-3;
  const R = g.radiusMM * 1e-3;
  const A = Math.PI * R * R;
  const Cd = (EPS0 * epsR * A) / dDiel;
  const Cg = (EPS0 * A) / dGap;
  const Ccell = (Cd * Cg) / (Cd + Cg);
  const k = 1e-4 / A; // пересчёт на 1 см²
  const dEff = dGap + dDiel / epsR;
  return {
    area_m2: A, dEff_m: dEff,
    C_diel: Cd, C_gas: Cg, C_cell: Ccell,
    perCm2: { C_diel: Cd * k, C_gas: Cg * k, C_cell: Ccell * k },
  };
}

// ───────────────────────────────────────────────────── MetricsPanel

const STYLE_ID = 'mp-styles-v1';
function injectStyles() {
  if (typeof document === 'undefined' || document.getElementById(STYLE_ID)) return;
  const el = document.createElement('style');
  el.id = STYLE_ID;
  el.textContent = `
.mp-root{display:flex;flex-direction:column;gap:6px;height:100%;min-height:0;
  font:13px Inter,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;color:${THEME.text};}
.mp-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:6px;flex:1 1 auto;min-height:0;}
.mp-tile{background:${THEME.surface2};border:1px solid ${THEME.border};border-radius:4px;
  padding:6px 8px;display:flex;flex-direction:column;justify-content:space-between;min-width:0;overflow:hidden;}
.mp-tile.mp-na{opacity:.55;}
.mp-label{flex:0 0 auto;font-size:10px;text-transform:uppercase;letter-spacing:.06em;color:${THEME.textMuted};
  white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
.mp-valrow{flex:0 0 auto;display:flex;align-items:baseline;gap:4px;min-width:0;}
.mp-value{font:600 clamp(14px,2.1vh,18px) ui-monospace,"JetBrains Mono",SFMono-Regular,Menlo,monospace;
  font-variant-numeric:tabular-nums;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
.mp-unit{font-size:11px;color:${THEME.textMuted};white-space:nowrap;}
.mp-note{flex:0 0 auto;font-size:10px;color:${THEME.textMuted};white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
.mp-spark{display:block;width:100%;height:16px;flex:1 1 auto;min-height:0;}
.mp-status{display:flex;flex-wrap:wrap;gap:10px;padding:4px 8px;background:${THEME.surface2};
  border:1px solid ${THEME.border};border-radius:4px;
  font:11px ui-monospace,"JetBrains Mono",SFMono-Regular,Menlo,monospace;font-variant-numeric:tabular-nums;
  color:${THEME.textSecondary};flex:0 0 auto;}
.mp-status .mp-k{color:${THEME.textMuted};}
.mp-status .mp-bad{color:${THEME.critical};}
.mp-status .mp-warn{color:${THEME.warning};}
`;
  document.head.appendChild(el);
}

/** Плитки метрик по умолчанию — набор UI_SPEC §2.5, пересчитанный на 2D-ячейку. */
export const DEFAULT_TILES = [
  { key: 'power',      label: 'Power',          fmt: (v) => fmtSI(v, 'W'),      spark: true },
  { key: 'specEnergy', label: 'Spec. energy',   fmt: (v) => fmtSI(v, 'J/L') },
  { key: 'ozone',      label: 'Ozone',          fmt: (v) => fmtFixed(v, 0, 'ppm'), spark: true },
  { key: 'maxEN',      label: 'max E/N',        fmt: (v) => fmtFixed(v, 0, 'Td'),
    thresholds: [1000, 3000, 1e4], spark: true },
  { key: 'Ipeak',      label: 'Peak current',   fmt: (v) => fmtSI(v, 'A') },
  { key: 'breakdowns', label: 'Breakdowns / T', fmt: (v) => fmtFixed(v, 0) },
  { key: 'C_diel',     label: 'C_diel',         fmt: (v) => fmtSI(v, 'F') },
  { key: 'C_cell',     label: 'C_cell',         fmt: (v) => fmtSI(v, 'F') },
  { key: 'energy',     label: 'Energy / period', fmt: (v) => fmtSI(v, 'J') },
];

/**
 * MetricsPanel — плитки + строка статуса. Только DOM, без canvas-осей
 * (спарклайны — маленькие canvas внутри плиток, с ring-буфером фиксированной длины).
 *
 *   const mp = new MetricsPanel(document.querySelector('#metrics'));
 *   mp.update({ power: 4.71, maxEN: 412, _status: { dt: 3.5e-11, cfl: 0.31 } });
 *
 * update() не аллоцирует: строки пишутся только когда изменились, спарклайны —
 * в предвыделенный Float32Array.
 */
export class MetricsPanel {
  constructor(el, opts = {}) {
    if (!el) throw new Error(tr('metrics.noContainer'));
    injectStyles();
    this.el = el;
    this.tiles = opts.tiles || DEFAULT_TILES;
    this.sparkLen = opts.sparkLen || 96;
    this.statusKeys = opts.statusKeys || [
      { key: 'dt', label: 'dt', fmt: (v) => fmtSI(v, 's') },
      { key: 'stepsPerFrame', label: 'steps/frame', fmt: (v) => fmtFixed(v, 0) },
      { key: 'frame', label: 'frame', fmt: (v) => v },
      { key: 'cfl', label: 'CFL', fmt: (v) => fmtFixed(v, 2), thresholds: [0.9, 1.0, 2.0] },
      { key: 'residual', label: 'residual', fmt: (v) => fmtExp(v), thresholds: [1e-8, 1e-6, 1e-4] },
      { key: 'limiter', label: 'limiter', fmt: (v) => v },
    ];
    this._items = new Map();
    this._status = new Map();
    this._build();
  }

  _build() {
    const root = document.createElement('div');
    root.className = 'mp-root';
    const grid = document.createElement('div');
    grid.className = 'mp-grid';
    for (const t of this.tiles) {
      const tile = document.createElement('div');
      tile.className = 'mp-tile';
      const lab = document.createElement('div');
      lab.className = 'mp-label';
      lab.textContent = t.label;
      const row = document.createElement('div');
      row.className = 'mp-valrow';
      const val = document.createElement('div');
      val.className = 'mp-value';
      val.textContent = '—';
      row.appendChild(val);
      const note = document.createElement('div');
      note.className = 'mp-note';
      note.textContent = t.note || '';
      tile.append(lab, row, note);
      let spark = null, ring = null;
      if (t.spark) {
        spark = document.createElement('canvas');
        spark.className = 'mp-spark';
        tile.appendChild(spark);
        ring = { buf: new Float32Array(this.sparkLen), n: 0, head: 0 };
      }
      grid.appendChild(tile);
      this._items.set(t.key, { spec: t, tile, val, note, spark, ring, last: null, lastNote: '' });
    }
    root.appendChild(grid);

    const st = document.createElement('div');
    st.className = 'mp-status';
    for (const s of this.statusKeys) {
      const wrap = document.createElement('span');
      const k = document.createElement('span');
      k.className = 'mp-k';
      k.textContent = s.label + ' ';
      const v = document.createElement('span');
      v.textContent = '—';
      wrap.append(k, v);
      st.appendChild(wrap);
      this._status.set(s.key, { spec: s, node: v, last: null });
    }
    root.appendChild(st);
    this.el.appendChild(root);
    this.root = root;
    this.statusEl = st;
  }

  /**
   * @param {object} data — { <tileKey>: number|null|{v,note,status}, _status:{...} }
   * Ключ, которого нет в data, не трогается вовсе (панель обновляется частями).
   */
  update(data) {
    if (!data) return;
    for (const [key, it] of this._items) {
      if (!(key in data)) continue;
      let raw = data[key];
      let note = null, status = null;
      if (raw && typeof raw === 'object') { note = raw.note ?? null; status = raw.status ?? null; raw = raw.v; }
      const na = raw == null || (typeof raw === 'number' && !Number.isFinite(raw));
      it.tile.classList.toggle('mp-na', na);
      const s = na ? '—' : it.spec.fmt(raw);
      if (s !== it.last) { it.val.textContent = s; it.last = s; }
      const cls = status || (it.spec.thresholds ? classify(raw, it.spec.thresholds) : 'ok');
      const col = STATUS_COLOR[cls] || THEME.text;
      if (it.val.style.color !== col) it.val.style.color = col;
      const nt = note ?? it.spec.note ?? '';
      if (nt !== it.lastNote) { it.note.textContent = nt; it.lastNote = nt; }
      if (it.ring && typeof raw === 'number' && Number.isFinite(raw)) {
        const r = it.ring;
        r.buf[r.head] = raw;
        r.head = (r.head + 1) % r.buf.length;
        if (r.n < r.buf.length) r.n++;
        this._drawSpark(it);
      }
    }
    const st = data._status;
    if (st) {
      for (const [key, it] of this._status) {
        if (!(key in st)) continue;
        const raw = st[key];
        const s = raw == null ? '—' : it.spec.fmt(raw);
        if (s !== it.last) { it.node.textContent = s; it.last = s; }
        const cls = it.spec.thresholds ? classify(raw, it.spec.thresholds) : 'ok';
        it.node.className = cls === 'ok' ? '' : (cls === 'warning' ? 'mp-warn' : 'mp-bad');
      }
    }
  }

  _drawSpark(it) {
    const c = it.spark, r = it.ring;
    const dpr = (typeof devicePixelRatio === 'number' ? devicePixelRatio : 1) || 1;
    const w = c.clientWidth || 60, h = c.clientHeight || 16;
    const W = Math.max(1, Math.round(w * dpr)), H = Math.max(1, Math.round(h * dpr));
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    const g = c.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    if (r.n < 2) return;
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < r.n; i++) {
      const v = r.buf[(r.head - r.n + i + r.buf.length) % r.buf.length];
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    if (!(hi > lo)) { hi = lo + 1; lo -= 1; }
    g.beginPath();
    for (let i = 0; i < r.n; i++) {
      const v = r.buf[(r.head - r.n + i + r.buf.length) % r.buf.length];
      const x = (i / (r.n - 1)) * w;
      const y = h - 1 - ((v - lo) / (hi - lo)) * (h - 2);
      if (i === 0) g.moveTo(x, y); else g.lineTo(x, y);
    }
    g.strokeStyle = THEME.accent;
    g.lineWidth = 1;
    g.stroke();
  }

  destroy() {
    if (this.root && this.root.parentNode) this.root.parentNode.removeChild(this.root);
    this._items.clear();
    this._status.clear();
  }
}

export default MetricsPanel;
