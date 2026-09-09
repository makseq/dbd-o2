// ui.js — all rendering and controls for the 1D DBD-O2 simulator.
// Canvas 2D + plain DOM only. No libraries, no CDN, no build step.
//
// SOLVER SWITCH ─────────────────────────────────────────────────────────────
//   default            → ./solver.js       (the real 1D fluid core)
//   ?mock=1 in the URL → ./mock-solver.js  (analytic stand-in, kept as an option:
//                                           instant, no physics, for UI work)
// Both satisfy the module CONTRACT; nothing below touches anything else.
// Differences the UI absorbs: the real core meshes the whole stack and flags gas
// cells in state.gasMask (see GAS adapter below), carries an extra species O3-,
// and puts condCurrent in the history ring buffer.
const USE_MOCK = new URLSearchParams(location.search).get('mock') === '1';

const solverModule = USE_MOCK
  ? await import('./mock-solver.js')
  : await import('./solver.js');
const DBDSolver = solverModule.DBDSolver || solverModule.default;

// ═══════════════════════════════════════════════════════════ constants ══════
const CSS = getComputedStyle(document.documentElement);
const C = (n, fb) => (CSS.getPropertyValue(n).trim() || fb);
const COL = {
  surface1: C('--surface-1', '#10141a'),
  surface2: C('--surface-2', '#161b23'),
  surface3: C('--surface-3', '#1c2430'),
  border: C('--border', '#242c38'),
  borderStrong: C('--border-strong', '#333f4f'),
  text: C('--text-primary', '#e8eaed'),
  text2: C('--text-secondary', '#a9b0bb'),
  muted: C('--text-muted', '#8b93a1'),
  grid: C('--grid', '#1e2530'),
  axis: C('--axis', '#39424f'),
  accent: C('--accent', '#4da3ff'),
  cursor: C('--cursor', '#c8ced8'),
  warning: C('--warning', '#fab219'),
  critical: C('--critical', '#d03b3b'),
  good: C('--good', '#0ca30c'),
};

// UI_SPEC §3.2 — slot order is validated, do not reshuffle.
// Slots 7-8 (O4p, O2a) are an extension for the 8-species contract; they are
// flagged `ext:true` and are hidden by default in CVD-safe mode.
const SPECIES = [
  { key: 'O2p', label: 'O₂⁺', color: '#e66767', dash: [], w: 2, cvd: true },
  { key: 'e', label: 'e⁻', color: '#3987e5', dash: [], w: 2.5, cvd: true },
  { key: 'Om', label: 'O⁻', color: '#d55181', dash: [6, 3], w: 2, cvd: true },
  { key: 'O2m', label: 'O₂⁻', color: '#c98500', dash: [2, 3], w: 2, cvd: false },
  { key: 'O', label: 'O', color: '#9085e9', dash: [8, 3, 2, 3], w: 2, cvd: false },
  { key: 'O3', label: 'O₃', color: '#199e70', dash: [1, 3], w: 2, cvd: true },
  { key: 'O4p', label: 'O₄⁺', color: '#c06fd8', dash: [5, 2, 1, 2], w: 1.5, cvd: false, ext: true },
  { key: 'O2a', label: 'O₂(a¹Δ)', color: '#5bc8c8', dash: [3, 2], w: 1.5, cvd: false, ext: true },
  // O3- is produced by the real core (O- + 2O2 -> O3- + O2 in ~1.5 ns at 1 atm) and is
  // the dominant afterglow anion; the mock omits it and it degrades to "absent".
  { key: 'O3m', label: 'O₃⁻', color: '#b0723a', dash: [4, 2, 1, 2], w: 1.5, cvd: false, ext: true },
];
const CHARGE = { e: -1, O2p: +1, O4p: +1, Om: -1, O2m: -1, O3m: -1, O: 0, O3: 0, O2a: 0 };

const INFERNO = ['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60', '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4'];
const VIRIDIS = ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'];
const DIVERGING = [
  [-1.0, '#9ec5f4'], [-0.6, '#3987e5'], [-0.25, '#256abf'],
  [0.0, '#383f48'],
  [0.25, '#a72f2f'], [0.6, '#e66767'], [1.0, '#f5a8a8'],
];

const hex2rgb = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];

function makeLUT(stops) {
  const lut = new Uint8ClampedArray(256 * 3);
  const rgb = stops.map(hex2rgb);
  const seg = rgb.length - 1;
  for (let i = 0; i < 256; i++) {
    const f = (i / 255) * seg;
    const k = Math.min(seg - 1, Math.floor(f));
    const u = f - k;
    for (let c = 0; c < 3; c++) lut[i * 3 + c] = rgb[k][c] + (rgb[k + 1][c] - rgb[k][c]) * u;
  }
  return lut;
}
function makeDivLUT(stops) {
  const lut = new Uint8ClampedArray(256 * 3);
  const pos = stops.map((s) => s[0]);
  const rgb = stops.map((s) => hex2rgb(s[1]));
  for (let i = 0; i < 256; i++) {
    const v = (i / 255) * 2 - 1;
    let k = 0;
    while (k < pos.length - 2 && v > pos[k + 1]) k++;
    const u = (v - pos[k]) / (pos[k + 1] - pos[k]);
    for (let c = 0; c < 3; c++) lut[i * 3 + c] = rgb[k][c] + (rgb[k + 1][c] - rgb[k][c]) * Math.max(0, Math.min(1, u));
  }
  return lut;
}
const LUT = { inferno: makeLUT(INFERNO), viridis: makeLUT(VIRIDIS), diverging: makeDivLUT(DIVERGING) };
const lutCss = (lut, i) => `rgb(${lut[i * 3]},${lut[i * 3 + 1]},${lut[i * 3 + 2]})`;

// ══════════════════════════════════════════════════════ gas-window adapter ══
// The real solver meshes the WHOLE stack (dielectric | gas | dielectric) and flags
// gas cells in `state.gasMask`; the mock meshes only the gap (gasMask all 1).
// Every gap-space drawing path must therefore go through this adapter instead of
// assuming `state.x` spans exactly the gap.
//   GAS.i0 / GAS.i1  — first / last gas cell in the global arrays
//   GAS.n            — number of gas cells
//   GAS.x0 / GAS.L   — left face and length of the gas region, metres
const GAS = { i0: 0, i1: 0, n: 1, x0: 0, L: 1, dxMin: 1 };
function updateGasWindow(st) {
  const m = st.gasMask, nx = st.x.length;
  let i0 = -1, i1 = -1;
  for (let i = 0; i < nx; i++) if (!m || m[i]) { if (i0 < 0) i0 = i; i1 = i; }
  if (i0 < 0) { i0 = 0; i1 = nx - 1; }
  GAS.i0 = i0; GAS.i1 = i1; GAS.n = i1 - i0 + 1;
  const xf = st.xFaces;
  GAS.x0 = xf ? xf[i0] : st.x[i0];
  GAS.L = Math.max(1e-12, (xf ? xf[i1 + 1] : st.x[i1]) - GAS.x0);
  let dmin = Infinity;
  if (xf) for (let i = i0; i <= i1; i++) dmin = Math.min(dmin, xf[i + 1] - xf[i]);
  GAS.dxMin = isFinite(dmin) ? dmin : GAS.L / GAS.n;
}
/** gas-local index k (0..GAS.n-1) → global array index */
const gidx = (k) => GAS.i0 + Math.min(GAS.n - 1, Math.max(0, k | 0));
/** global array index → normalised position across the gap, 0..1 */
const gpos = (st, i) => (st.x[i] - GAS.x0) / GAS.L;

// ═══════════════════════════════════════════════════════════ formatting ═════
const SUP = { '-': '⁻', '+': '', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const supStr = (n) => String(n).split('').map((c) => SUP[c] ?? c).join('');
function sci(v, digits = 2) {
  if (!isFinite(v)) return '–';
  if (v === 0) return '0';
  const s = v < 0 ? '−' : '';
  const a = Math.abs(v);
  const e = Math.floor(Math.log10(a));
  const m = a / Math.pow(10, e);
  // Number(...) убирает экспоненциальную форму, которую toPrecision(3) навязывает
  // при порядке >= 3: (4090).toPrecision(3) === '4.09e+3', а UI_SPEC §5.1 запрещает
  // нотацию `e+3` вне полей ввода.
  if (e >= -2 && e <= 3) return s + Number(a.toPrecision(3)).toString();
  return `${s}${m.toFixed(digits)}·10${supStr(e)}`;
}
/** Автомасштаб единиц для плиток: 4090 mA -> «4.09 A». */
function withUnit(v, unit, steps) {
  let u = unit, x = v;
  for (const [mul, name] of steps) {
    if (Math.abs(v) >= mul) { x = v / mul; u = name; }
  }
  return { v: x, unit: u };
}
const sig3 = (v) => (isFinite(v) ? (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0) ? sci(v) : Number(v.toPrecision(3)).toString()) : '–');
const NBSP = ' ';

// ═══════════════════════════════════════════════════════════ canvases ═══════
const dpr = () => Math.min(window.devicePixelRatio || 1, 2);
const $ = (id) => document.getElementById(id);
const canvases = {};
function reg(id) {
  const el = $(id);
  const ctx = el.getContext('2d', { alpha: false });
  const o = { el, ctx, w: 0, h: 0 };
  canvases[id] = o;
  return o;
}
function fit(o) {
  const r = o.el.getBoundingClientRect();
  const d = dpr();
  const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
  if (o.w === w && o.h === h && o.dpr === d) return false;
  o.w = w; o.h = h; o.dpr = d;
  o.el.width = Math.round(w * d); o.el.height = Math.round(h * d);
  o.ctx.setTransform(d, 0, 0, d, 0, 0);
  return true;
}
function clearBg(o, color = COL.surface1) {
  o.ctx.setTransform(o.dpr, 0, 0, o.dpr, 0, 0);
  o.ctx.fillStyle = color;
  o.ctx.fillRect(0, 0, o.w, o.h);
}
function text(ctx, s, x, y, { color = COL.muted, size = 10, align = 'left', base = 'alphabetic', mono = true, weight = '' } = {}) {
  ctx.fillStyle = color;
  ctx.font = `${weight} ${size}px ${mono ? 'JetBrains Mono, ui-monospace, monospace' : 'Inter, system-ui, sans-serif'}`.trim();
  ctx.textAlign = align; ctx.textBaseline = base;
  ctx.fillText(s, x, y);
}

const cvGap = reg('cvGap'), cvWave = reg('cvWave'), cvQV = reg('cvQV'),
  cvProf = reg('cvProf'), cvMacro = reg('cvMacro'), cvMicro = reg('cvMicro');

// ═══════════════════════════════════════════════════════════ presets ════════
const PRESETS = {
  filamentary: {
    label: 'Filamentary',
    p: { U0kV: 10, freqKHz: 10, gapMM: 1.0, dielMM1: 0.5, dielMM2: 0.5, epsR: 9, gamma: 0.02, seedDensity: 1e13, pressureTorr: 760, tempK: 300, areaCM2: 1, nCells: 400, betaGrid: 3.2, ballastOhm: 0, mode: 'default' },
    // lensPeriods = ширина окна осциллографа в ПЕРИОДАХ. Раньше пресеты задавали
    // абсолютный lensLog (2 мкс при периоде 100 мкс), и на первом экране была видна
    // прямая линия без единого импульса тока.
    view: { field: 'ne', log: true, iLog: true, lensPeriods: 1, profMode: 'x' },
  },
  townsend: {
    label: 'Townsend / homogeneous',
    p: { U0kV: 5, freqKHz: 5, gapMM: 2.0, dielMM1: 0.5, dielMM2: 0.5, epsR: 9, gamma: 0.2, seedDensity: 1e15, pressureTorr: 380, tempK: 300, areaCM2: 1, nCells: 300, ballastOhm: 0, mode: 'demo' },
    view: { field: 'sion', log: true, iLog: true, lensPeriods: 1, profMode: 'x' },
  },
  ozonizer: {
    label: 'Ozonizer',
    p: { U0kV: 12, freqKHz: 25, gapMM: 0.6, dielMM1: 0.5, dielMM2: 0.5, epsR: 9, gamma: 0.02, seedDensity: 1e13, pressureTorr: 760, tempK: 300, areaCM2: 10, nCells: 300, ballastOhm: 0, mode: 'demo' },
    view: { field: 'ne', log: true, iLog: true, lensPeriods: 1, profMode: 't' },
  },
};

// ═══════════════════════════════════════════════════════════ view state ═════
const V = {
  preset: 'filamentary',
  field: 'ne',
  log: true,
  auto: true,
  gapMode: 'snapshot',
  iLog: true,
  showItot: false,
  qvWindow: 1,
  profMode: 'x',
  cvdSafe: false,
  running: false,
  speedLog: -3.7,      // sim-seconds per wall-second, log10
  lensLog: -4,         // lens window width = 10^lensLog s (set from lensPeriods)
  lensPeriods: 1,      // ширина окна в периодах — именно это фиксируется при смене f
  cursorT: 0,
  followCursor: true,
  hidden: new Set(),
  solo: null,
  stopAtT: null,
  steady: false,
};

const solver = new DBDSolver(Object.assign({}, PRESETS.filamentary.p));
window.__solver = solver; // debugging handle only

// ═══════════════════════════════════════════════════════════ ui buffers ═════
const MAXW = 4096;
const buf = {
  min: new Float64Array(MAXW), max: new Float64Array(MAXW), last: new Float64Array(MAXW),
  cnt: new Int32Array(MAXW),
  min2: new Float64Array(MAXW), max2: new Float64Array(MAXW), last2: new Float64Array(MAXW),
  cnt2: new Int32Array(MAXW),
  min3: new Float64Array(MAXW), max3: new Float64Array(MAXW), last3: new Float64Array(MAXW),
  ug: new Float64Array(MAXW),
  glow: new Float64Array(4096),          // afterglow buffer G(x)
  evT: new Float64Array(8192), evP: new Int8Array(8192), evN: 0, evHead: 0,
  spark: {},                              // per-metric sparkline ring
  slowT: new Float64Array(2048), slowO3: new Float64Array(2048), slowN: 0,
};
const SPARK_N = 64;
const METRICS = [
  ['power', 'Power', 'W'], ['spec', 'Spec. energy', 'J/L'], ['o3', 'Ozone', 'ppm'],
  ['en', 'max E/N', 'Td'], ['ipk', 'Peak current', 'mA'], ['bd', 'Breakdowns/T', ''],
  ['cdiel', 'C_diel', 'pF'], ['ccell', 'C_cell', 'pF'], ['eper', 'Energy/period', 'mJ'],
];
for (const [k] of METRICS) buf.spark[k] = { a: new Float64Array(SPARK_N), n: 0, h: 0 };
function sparkPush(k, v) {
  const s = buf.spark[k];
  s.a[s.h] = v; s.h = (s.h + 1) % SPARK_N; if (s.n < SPARK_N) s.n++;
}

// x–t streak: rolling offscreen raster
const streak = { cv: document.createElement('canvas'), ctx: null, row: 0, h: 0, w: 0, img: null, t0: 0 };
streak.ctx = streak.cv.getContext('2d', { willReadFrequently: false });

// ═══════════════════════════════════════════════════════════ controls ═══════
const CONTROLS = [
  { k: 'U0kV', label: 'U₀ amplitude', min: 1, max: 20, step: 0.1, unit: 'kV', log: false },
  { k: 'freqKHz', label: 'f frequency', min: 0.5, max: 100, step: 0.01, unit: 'kHz', log: true },
  { k: 'gapMM', label: 'd_gap', min: 0.2, max: 5, step: 0.01, unit: 'mm', log: true, geom: true },
  { k: 'dielMM1', label: 'd_diel ×2', min: 0.1, max: 2, step: 0.05, unit: 'mm', geom: true, twin: 'dielMM2' },
  { k: 'epsR', label: 'ε_r', min: 1, max: 30, step: 0.5, unit: '' },
  { k: 'gamma', label: 'γ sec. emis.', min: 0, max: 0.3, step: 0.005, unit: '' },
  { k: 'seedDensity', label: 'n₀ seed', min: 1e10, max: 1e16, step: 0.01, unit: 'm⁻³', log: true },
  { k: 'pressureTorr', label: 'p pressure', min: 76, max: 1520, step: 1, unit: 'Torr' },
  { k: 'tempK', label: 'T gas', min: 250, max: 600, step: 10, unit: 'K' },
  { k: 'areaCM2', label: 'S area', min: 0.1, max: 100, step: 0.01, unit: 'cm²', log: true },
  { k: 'ballastOhm', label: 'R_series', min: 0, max: 1e7, step: 1, unit: 'Ω' },
  { k: 'nCells', label: 'N_x cells', min: 100, max: 2000, step: 1, unit: '', log: true, hard: true },
];
const ctlEls = {};

function buildControls() {
  const host = $('controls');
  host.textContent = '';
  for (const c of CONTROLS) {
    const row = document.createElement('div');
    row.className = 'ctl';
    const lab = document.createElement('label');
    lab.textContent = c.label;
    const rng = document.createElement('input');
    rng.type = 'range';
    if (c.log) { rng.min = Math.log10(c.min); rng.max = Math.log10(c.max); rng.step = 0.005; }
    else { rng.min = c.min; rng.max = c.max; rng.step = c.step; }
    const num = document.createElement('input');
    num.type = 'number'; num.className = 'num'; num.step = c.step;
    num.min = c.min; num.max = c.max;
    const id = 'ctl_' + c.k;
    rng.id = id; lab.htmlFor = id;
    rng.setAttribute('aria-label', `${c.label}${c.unit ? ', ' + c.unit : ''}`);
    num.setAttribute('aria-label', `${c.label} value${c.unit ? ', ' + c.unit : ''}`);
    num.title = `${c.label}: ${c.min}…${c.max}${c.unit ? ' ' + c.unit : ''}`;
    // единица измерения обязательна (UI_SPEC §6.1): без неё «760» и «1.0e+13»
    // не читаются — непонятно, Торр это или атм, м^-3 или см^-3
    const u = document.createElement('span');
    u.className = 'unit'; u.textContent = c.unit || '';
    row.append(lab, rng, num, u);
    host.append(row);
    ctlEls[c.k] = { rng, num, c };

    const apply = (v, live) => {
      v = Math.max(c.min, Math.min(c.max, v));
      const patch = { [c.k]: c.k === 'nCells' ? Math.round(v) : v };
      if (c.twin) patch[c.twin] = v;
      solver.setParams(patch);
      syncControl(c.k);
      // окно осциллографа привязано к периоду: меняя f, пользователь ожидает
      // увидеть тот же «один период», а не 2 мкс из 2 мс
      if (c.k === 'freqKHz') setLensPeriods(V.lensPeriods);
      onParamsChanged(c.geom || c.hard);
      if (!live) refreshAll();
    };
    rng.addEventListener('input', () => apply(c.log ? Math.pow(10, +rng.value) : +rng.value, true));
    num.addEventListener('change', () => apply(+num.value, false));
  }
  // `params.mode` is the solver's ACCURACY preset (CFL numbers + default cell count),
  // not a waveform shape — the core drives a pure sinusoid only, so offering
  // square/pulsed here would have been a lie. Changing it rebuilds the grid.
  const sec = document.createElement('div');
  sec.className = 'ctl-sec';
  sec.textContent = 'Accuracy / speed';
  host.append(sec);
  const row = document.createElement('div');
  row.className = 'ctl';
  const lab = document.createElement('label'); lab.textContent = 'mode';
  const sel = document.createElement('select');
  const MODE_HINT = { demo: 'demo — fastest, CFL 0.8', default: 'default — CFL 0.4', accurate: 'accurate — CFL 0.2, slow' };
  for (const m of ['demo', 'default', 'accurate']) {
    const o = document.createElement('option'); o.value = m; o.textContent = MODE_HINT[m]; sel.append(o);
  }
  sel.value = solver.params.mode;
  sel.addEventListener('change', () => {
    solver.setParams({ mode: sel.value, nCells: null });
    syncAllControls(); onParamsChanged(true);
  });
  sel.setAttribute('aria-label', 'Solver accuracy preset');
  row.append(lab, sel, document.createElement('span'), document.createElement('span'));
  host.append(row);
  syncAllControls();
}
function syncControl(k) {
  const e = ctlEls[k]; if (!e) return;
  const v = solver.params[k];
  e.rng.value = e.c.log ? Math.log10(Math.max(e.c.min, v)) : v;
  e.num.value = (Math.abs(v) >= 1e5 || (v !== 0 && Math.abs(v) < 1e-3))
    ? v.toExponential(1) : Number(v.toPrecision(4));
}
function syncAllControls() { for (const c of CONTROLS) syncControl(c.k); }

function onParamsChanged(hard) {
  const p = solver.params;
  $('hdrSummary').textContent = `U₀ ${p.U0kV.toFixed(1)} kV · f ${p.freqKHz.toFixed(1)} kHz · gap ${p.gapMM.toFixed(2)} mm · ε_r ${p.epsR}`;
  updateGridInfo();
  // Смена геометрии/сетки ВНУТРИ солвера делает reset() (t=0, history пустая).
  // Значит, весь UI-стейт, привязанный ко времени прогона, обязан обнулиться —
  // иначе evScanT глушит новые события пробоя, cumEnergyJ смешивает два прогона,
  // а панель t-профилей показывает точки из предыдущей симуляции.
  if (hard) { resetUiState(); announce('geometry changed — simulation restarted from t = 0'); }
}
/** Единая точка сброса UI-состояния, привязанного ко времени прогона. */
function resetUiState() {
  buf.evN = 0; buf.evHead = 0; buf.slowN = 0;
  evScanT = 0; evPrevI = 0;
  buf.glow.fill(0); if (buf.glowAbs) buf.glowAbs.fill(0);
  glowMaxHist = 1e-30; sigmaRef = 1e-9;
  streak.row = 0; streak.w = 0; streak.t0 = solver.state.t; streak.tRow = solver.state.t;
  streak.ctx.fillStyle = '#000004';
  tsT.n = 0; tsT.h = 0; tsLastWall = 0;
  for (const s of SPECIES) { tsBuf[s.key].n = 0; tsBuf[s.key].h = 0; }
  for (const [k] of METRICS) { buf.spark[k].n = 0; buf.spark[k].h = 0; }
  cumEnergyJ = 0; lastPeriodIdx = -1;
  gridWorst = Infinity;
  V.cursorT = solver.state.t; V.followCursor = true; V.stopAtT = null;
  steadyHist = []; steadyCount = 0;
}
// худшее (минимальное) значение λ_D/Δx с момента сброса — индикатор обязан
// показывать НЕ мгновенное значение в тихой фазе, а самый плохой случай
let gridWorst = Infinity;
/** Дешёвая (O(N)) выборка худшего λ_D/Δx — каждый кадр, а не раз в период. */
function sampleGridWorst(st) {
  let neMax = 1;
  for (let i = GAS.i0; i <= GAS.i1; i++) if (st.n.e[i] > neMax) neMax = st.n.e[i];
  const lD = Math.sqrt(8.854e-12 * 2 * 1.602e-19 / (neMax * 1.602e-19 * 1.602e-19));
  const ratio = lD / GAS.dxMin;
  if (ratio < gridWorst) gridWorst = ratio;
}
function updateGridInfo() {
  const st = solver.state;
  updateGasWindow(st);
  const dx = GAS.dxMin;
  // Debye length from the peak electron density in the gap (Te ~ 2 eV assumed for the readout)
  let neMax = 1;
  for (let i = GAS.i0; i <= GAS.i1; i++) if (st.n.e[i] > neMax) neMax = st.n.e[i];
  const lD = Math.sqrt(8.854e-12 * 2 * 1.602e-19 / (neMax * 1.602e-19 * 1.602e-19));
  const ratio = lD / dx;
  if (ratio < gridWorst) gridWorst = ratio;
  const el = $('gridInfo');
  // UI_SPEC §5.3: λ_D < 2Δx — уже критично; судим по ХУДШЕМУ за прогон, потому что
  // выборка «сейчас» систематически попадает в тихую фазу (n_e у пола) и показывала
  // зелёную галочку в кадрах, где профиль n_e на 6 порядков выше.
  const bad = gridWorst < 1, warn = gridWorst < 2;
  el.style.color = bad ? COL.critical : warn ? COL.warning : COL.muted;
  el.textContent = `${bad ? '✕ ' : warn ? '⚠ ' : '✓ '}Δx = ${sig3(dx * 1e6)}${NBSP}µm · ` +
    `λ_D(now) = ${sig3(lD * 1e6)}${NBSP}µm · λ_D/Δx now ${sig3(ratio)} / worst ${sig3(gridWorst)}`;
}

// ═══════════════════════════════════════════════════════════ metrics DOM ════
const tileEls = {};
function buildMetrics() {
  const host = $('metrics');
  host.textContent = '';
  for (const [k, label, unit] of METRICS) {
    const d = document.createElement('div');
    d.className = 'tile';
    const kk = document.createElement('div'); kk.className = 'k'; kk.textContent = label;
    const v = document.createElement('div'); v.className = 'v mono';
    v.innerHTML = `<span class="n">–</span>${unit ? `<span class="u">${unit}</span>` : ''}`;
    const sp = document.createElement('canvas'); sp.height = 14;
    d.append(kk, v, sp);
    host.append(d);
    tileEls[k] = { n: v.querySelector('.n'), u: v.querySelector('.u'), spark: sp, ctx: sp.getContext('2d') };
  }
}
function drawSpark(t, key) {
  const s = buf.spark[key];
  const cv = t.spark, ctx = t.ctx;
  const w = Math.max(20, cv.clientWidth || 60), h = 14, d = dpr();
  if (cv.width !== Math.round(w * d)) { cv.width = Math.round(w * d); cv.height = Math.round(h * d); }
  ctx.setTransform(d, 0, 0, d, 0, 0);
  ctx.clearRect(0, 0, w, h);
  if (s.n < 2) return;
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < s.n; i++) { const v = s.a[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
  if (hi - lo < 1e-30) { hi = lo + 1; }
  ctx.beginPath();
  for (let i = 0; i < s.n; i++) {
    const idx = (s.h - s.n + i + SPARK_N) % SPARK_N;
    const x = (i / (s.n - 1)) * (w - 1);
    const y = h - 1 - ((s.a[idx] - lo) / (hi - lo)) * (h - 2);
    i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
  }
  ctx.strokeStyle = COL.accent; ctx.lineWidth = 1; ctx.stroke();
}

// ═══════════════════════════════════════════════════════════ legend ═════════
function buildLegend() {
  const host = $('legend');
  host.textContent = '';
  for (const s of SPECIES) {
    const b = document.createElement('button');
    b.className = 'lg'; b.dataset.k = s.key; b.type = 'button';
    b.title = 'click = solo · alt+click = hide';
    const sw = document.createElement('span');
    sw.className = 'sw';
    sw.style.borderTopColor = s.color;
    sw.style.borderTopStyle = s.dash.length ? 'dashed' : 'solid';
    b.append(sw, document.createTextNode(s.label));
    b.addEventListener('click', (e) => {
      if (e.altKey) {
        V.hidden.has(s.key) ? V.hidden.delete(s.key) : V.hidden.add(s.key);
      } else {
        V.solo = V.solo === s.key ? null : s.key;
      }
      syncLegend(); drawProfiles();
    });
    host.append(b);
  }
  const hint = document.createElement('span');
  hint.className = 'muted';
  hint.style.fontSize = '10px';
  hint.textContent = 'solo: click · hide: alt-click';
  host.append(hint);
  syncLegend();
}
function syncLegend() {
  for (const el of $('legend').querySelectorAll('.lg')) {
    const k = el.dataset.k;
    el.classList.toggle('off', !visible(k));
    el.classList.toggle('solo', V.solo === k);
  }
}
const visible = (k) => {
  const s = SPECIES.find((x) => x.key === k);
  if (V.cvdSafe && !s.cvd) return false;
  if (V.hidden.has(k)) return false;
  if (V.solo && V.solo !== k) return false;
  return true;
};

// ═══════════════════════════════════════════════════════════ field access ═══
const FIELDS = {
  ne: { name: 'n_e', unit: 'm⁻³', log: true, signed: false, get: (st) => st.n.e },
  rho: { name: 'ρ', unit: 'C·m⁻³', log: false, signed: true, get: (st) => st.rho },
  // paint: массив, который РЕАЛЬНО красится (для S_ion это послесвечение в
  // абсолютных единицах); шкала и подписи colorbar строятся по нему же
  sion: {
    name: 'S_ion (afterglow)', unit: 'm⁻³s⁻¹', log: true, signed: false,
    get: (st) => st.ionizRate, paint: () => buf.glowAbs || null,
  },
  Eabs: { name: '|E|', unit: 'V/m', log: false, signed: false, get: (st) => st.E, abs: true },
  EN: { name: 'E/N', unit: 'Td', log: false, signed: false, get: (st) => st.EN },
};
const scaleSmooth = {}; // exponentially smoothed auto-scale per field

function fieldRange(st, key) {
  const f = FIELDS[key];
  const a = (f.paint && f.paint(st)) || f.get(st);
  let lo = Infinity, hi = -Infinity;
  // gas cells only — the dielectric slabs carry no plasma and a very different |E|,
  // and would otherwise dominate the auto-scale
  for (let i = GAS.i0; i <= GAS.i1; i++) {
    const v = f.abs ? Math.abs(a[i]) : a[i];
    if (!isFinite(v)) continue;
    if (v < lo) lo = v; if (v > hi) hi = v;
  }
  if (!isFinite(lo)) { lo = 0; hi = 1; }
  if (f.signed) { const m = Math.max(Math.abs(lo), Math.abs(hi), 1e-30); lo = -m; hi = m; }
  else if (V.log && f.log) { hi = Math.max(hi, 1e-30); lo = Math.max(hi * 1e-7, 1e-30); }
  else { lo = Math.min(0, lo); if (hi - lo < 1e-30) hi = lo + 1; }
  const s = scaleSmooth[key] || (scaleSmooth[key] = { lo, hi });
  if (V.auto) {
    const a1 = 0.12;
    s.lo += (lo - s.lo) * a1; s.hi += (hi - s.hi) * a1;
  }
  return V.auto ? s : { lo, hi };
}
function norm(v, lo, hi, logScale, signed) {
  if (signed) {
    const m = Math.max(Math.abs(lo), Math.abs(hi), 1e-30);
    return Math.max(0, Math.min(1, (v / m + 1) / 2));
  }
  if (logScale) {
    const l = Math.log10(Math.max(v, 1e-30)), a = Math.log10(Math.max(lo, 1e-30)), b = Math.log10(Math.max(hi, 1e-30));
    return b - a < 1e-9 ? 0 : Math.max(0, Math.min(1, (l - a) / (b - a)));
  }
  return hi - lo < 1e-30 ? 0 : Math.max(0, Math.min(1, (v - lo) / (hi - lo)));
}

// ═══════════════════════════════════════════════════════════ B. gap view ════
let gapImg = null, gapImgW = 0;
const gapTmp = document.createElement('canvas');
const gapTmpCtx = gapTmp.getContext('2d');
let lastWallT = performance.now();

// buf.glowAbs — послесвечение в АБСОЛЮТНЫХ единицах (м⁻³с⁻¹). Именно этот массив
// и красится, и задаёт границы шкалы: раньше в полосу шёл glow*glowMaxHist, а шкалу
// считали по мгновенному st.ionizRate, из-за чего между импульсами крашеные значения
// оказывались на 1–3 декады выше верха шкалы, всё насыщалось в жёлтый и филамент
// пропадал, а подписи colorbar относились к другой величине.
function updateGlow(st, dtWall) {
  const nx = st.ionizRate.length;
  if (buf.glow.length < nx) buf.glow = new Float64Array(nx);
  if (!buf.glowAbs || buf.glowAbs.length < nx) buf.glowAbs = new Float64Array(nx);
  let mx = 1e-30;
  for (let i = 0; i < nx; i++) if (st.ionizRate[i] > mx) mx = st.ionizRate[i];
  glowMaxHist = Math.max(glowMaxHist * 0.9995, mx);
  const decay = Math.exp(-dtWall / 0.25); // τ_g = 250 ms wall time (UI_SPEC §4.2)
  let g = 0;
  for (let i = 0; i < nx; i++) {
    const s = st.ionizRate[i] / glowMaxHist;
    const v = Math.max(s, buf.glow[i] * decay);
    buf.glow[i] = v;
    buf.glowAbs[i] = v * glowMaxHist;
    if (v > g) g = v;
  }
  return g;
}
let glowMaxHist = 1e-30;

function drawGap() {
  const o = cvGap; fit(o);
  const ctx = o.ctx; clearBg(o);
  const st = solver.state;
  const p = solver.params;
  const W = o.w, H = o.h;
  const padL = 34, padR = 58, padT = 6, padB = 40;
  const plotW = Math.max(10, W - padL - padR);
  const bandH = Math.max(40, Math.min(150, H - padT - padB - 34));
  const metalW = 12;
  const totalMM = p.dielMM1 + p.gapMM + p.dielMM2;
  const inner = plotW - 2 * metalW;
  const wd1 = (inner * p.dielMM1) / totalMM;
  const wg = (inner * p.gapMM) / totalMM;
  const wd2 = (inner * p.dielMM2) / totalMM;
  const x0 = padL, y0 = padT + 12;
  const gx0 = x0 + metalW + wd1;

  // metal + dielectric slabs with hatching (texture survives forced-colors / print)
  const hatch = (x, w, angle, color, bg) => {
    ctx.save(); ctx.beginPath(); ctx.rect(x, y0, w, bandH); ctx.clip();
    ctx.fillStyle = bg; ctx.fillRect(x, y0, w, bandH);
    ctx.strokeStyle = color; ctx.lineWidth = 1;
    ctx.beginPath();
    if (angle === 90) { for (let i = 0; i < w; i += 4) { ctx.moveTo(x + i, y0); ctx.lineTo(x + i, y0 + bandH); } }
    else { for (let i = -bandH; i < w + bandH; i += 6) { ctx.moveTo(x + i, y0 + bandH); ctx.lineTo(x + i + bandH, y0); } }
    ctx.stroke(); ctx.restore();
  };
  hatch(x0, metalW, 90, '#4b5566', '#3a4250');
  hatch(x0 + metalW, wd1, 45, '#2a3644', '#1c2430');
  hatch(gx0 + wg, wd2, 45, '#2a3644', '#1c2430');
  hatch(gx0 + wg + wd2, metalW, 90, '#4b5566', '#3a4250');

  const f = FIELDS[V.field];
  const arr = (f.paint && f.paint(st)) || f.get(st);
  const nx = GAS.n;                 // cells drawn in the gap band (gas only)
  const px2i = (px) => gidx(Math.floor((px / iw) * nx));
  const { lo, hi } = fieldRange(st, V.field);
  const useLog = V.log && f.log && !f.signed;
  const lut = f.signed ? LUT.diverging : LUT.inferno;
  const iw = Math.max(1, Math.round(wg));

  if (V.gapMode === 'snapshot') {
    if (!gapImg || gapImgW !== iw) { gapImg = ctx.createImageData(iw, 1); gapImgW = iw; }
    const d = gapImg.data;
    for (let px = 0; px < iw; px++) {
      const i = px2i(px);
      const raw = f.abs ? Math.abs(arr[i]) : arr[i];
      const t = norm(raw, lo, hi, useLog, f.signed);
      const c = Math.round(t * 255) * 3;
      d[px * 4] = lut[c]; d[px * 4 + 1] = lut[c + 1]; d[px * 4 + 2] = lut[c + 2]; d[px * 4 + 3] = 255;
    }
    if (gapTmp.width !== iw) { gapTmp.width = iw; gapTmp.height = 1; }
    gapTmpCtx.putImageData(gapImg, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(gapTmp, gx0, y0, wg, bandH);

    // quantitative overlay line (colormap = recognition, line = reading)
    ctx.save(); ctx.beginPath(); ctx.rect(gx0, y0, wg, bandH); ctx.clip();
    ctx.strokeStyle = 'rgba(232,234,237,0.6)'; ctx.lineWidth = 1.5; ctx.beginPath();
    for (let px = 0; px < iw; px++) {
      const i = px2i(px);
      const raw = f.abs ? Math.abs(arr[i]) : arr[i];
      const t = norm(raw, lo, hi, useLog, f.signed);
      const y = y0 + bandH - t * (bandH - 4) - 2;
      px ? ctx.lineTo(gx0 + px, y) : ctx.moveTo(gx0 + px, y);
    }
    ctx.stroke();
    if (f.signed) { // ρ = 0 isoline, dashed white — never colour alone
      ctx.setLineDash([2, 2]); ctx.strokeStyle = 'rgba(255,255,255,.75)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(gx0, y0 + bandH / 2); ctx.lineTo(gx0 + wg, y0 + bandH / 2); ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.restore();

    // filament highlight
    let gmax = 0; for (let i = GAS.i0; i <= GAS.i1; i++) if (buf.glow[i] > gmax) gmax = buf.glow[i];
    if (gmax > 0.05) {
      let a = -1, b = -1;
      for (let i = GAS.i0; i <= GAS.i1; i++) if (buf.glow[i] > 0.05 * gmax) { if (a < 0) a = i; b = i; }
      if (a >= 0) {
        const xa = gx0 + gpos(st, a) * wg;
        const xb = gx0 + Math.min(1, gpos(st, b) + 1 / nx) * wg;
        ctx.strokeStyle = '#ffd166'; ctx.lineWidth = 1;
        ctx.shadowColor = 'rgba(255,209,102,.45)'; ctx.shadowBlur = 12;
        ctx.strokeRect(xa, y0 - 1, Math.max(2, xb - xa), bandH + 2);
        ctx.shadowBlur = 0;
      }
    }
  } else {
    // x–t streak: rolling raster, y = time (down = past)
    const sh = Math.round(bandH);
    if (streak.w !== iw || streak.h !== sh) {
      streak.w = iw; streak.h = sh; streak.cv.width = iw; streak.cv.height = sh;
      streak.ctx.fillStyle = '#000004'; streak.ctx.fillRect(0, 0, iw, sh);
      streak.row = 0; streak.img = streak.ctx.createImageData(iw, 1); streak.tRow = st.t;
    }
    // Растр по РАВНОМЕРНОЙ шкале времени: одна строка = lens/bandH секунд симуляции,
    // а не «один кадр». При адаптивном шаге (1e-12…5e-9 с) и бюджете 12 мс/кадр
    // соседние строки отличались по dt_sim в 1e5 раз — вертикальная ось была
    // нечитаемой, и наклон катодонаправленной волны (главный смысл этой панели)
    // измерить по картинке было нельзя.
    const lens = Math.pow(10, V.lensLog);
    const dtRow = lens / sh;
    if (!(streak.tRow > 0) || st.t < streak.tRow) streak.tRow = st.t;
    let rows = Math.floor((st.t - streak.tRow) / dtRow);
    if (rows > sh) rows = sh;                 // «догоняем» не больше одного экрана
    if (rows > 0) {
      const d = streak.img.data;
      for (let px = 0; px < iw; px++) {
        const i = px2i(px);
        const raw = f.abs ? Math.abs(arr[i]) : arr[i];
        const t = norm(raw, lo, hi, useLog, f.signed);
        const c = Math.round(t * 255) * 3;
        d[px * 4] = lut[c]; d[px * 4 + 1] = lut[c + 1]; d[px * 4 + 2] = lut[c + 2]; d[px * 4 + 3] = 255;
      }
      for (let k = 0; k < rows; k++) {
        streak.ctx.putImageData(streak.img, 0, streak.row);
        streak.row = (streak.row + 1) % sh;
      }
      streak.tRow += rows * dtRow;
    }
    ctx.imageSmoothingEnabled = false;
    const r = streak.row;
    ctx.drawImage(streak.cv, 0, r, iw, sh - r, gx0, y0, wg, ((sh - r) / sh) * bandH);
    if (r > 0) ctx.drawImage(streak.cv, 0, 0, iw, r, gx0, y0 + ((sh - r) / sh) * bandH, wg, (r / sh) * bandH);
    // ось времени слева: 0 = «сейчас» внизу, вверх — прошлое на глубину lens
    ctx.strokeStyle = COL.axis; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(gx0 - 1, y0); ctx.lineTo(gx0 - 1, y0 + bandH); ctx.stroke();
    for (let k = 0; k <= 4; k++) {
      const yy = y0 + (k / 4) * bandH;
      ctx.beginPath(); ctx.moveTo(gx0 - 4, yy); ctx.lineTo(gx0 - 1, yy); ctx.stroke();
      text(ctx, k === 4 ? 'now' : `−${fmtTime(lens * (1 - k / 4), 1)}`,
        gx0 - 6, yy, { size: 8.5, align: 'right', base: 'middle' });
    }
    // граница полупериода — опорная линия для отсчёта наклона фронта
    const Thalf = 0.5 / (solver.params.freqKHz * 1e3);
    if (Thalf < lens) {
      const yh = y0 + bandH * (1 - Thalf / lens);
      ctx.setLineDash([3, 3]); ctx.strokeStyle = 'rgba(232,234,237,.35)';
      ctx.beginPath(); ctx.moveTo(gx0, yh); ctx.lineTo(gx0 + wg, yh); ctx.stroke();
      ctx.setLineDash([]);
      text(ctx, 'T/2', gx0 + wg - 2, yh - 2, { size: 8.5, align: 'right', base: 'bottom', color: COL.text2 });
    }
    text(ctx, `↓ time · 1 row = ${fmtTime(dtRow, 1)} · v_front ≈ ${frontSpeedText()}`,
      gx0 + 4, y0 + 10, { color: '#e8eaed', size: 10 });
  }

  // σ buses at both gas/dielectric interfaces (+/− glyphs = secondary encoding)
  const drawSigma = (xpos, sigma, side) => {
    const mag = Math.min(1, Math.abs(sigma) / (sigmaRef || 1e-9));
    const col = sigma >= 0 ? '#e66767' : '#3987e5';
    const hh = Math.max(2, mag * bandH);
    ctx.globalAlpha = 0.15 + 0.85 * mag;
    ctx.fillStyle = col;
    ctx.fillRect(xpos - 3, y0 + bandH - hh, 6, hh);
    ctx.globalAlpha = 1;
    ctx.strokeStyle = COL.borderStrong; ctx.lineWidth = 1;
    ctx.strokeRect(xpos - 3, y0, 6, bandH);
    const glyphs = Math.max(1, Math.round(mag * 6));
    for (let i = 0; i < glyphs; i++) {
      text(ctx, sigma >= 0 ? '+' : '−', xpos, y0 + bandH - 6 - i * 9, { color: '#fff', size: 8, align: 'center' });
    }
    text(ctx, `${sig3(sigma * 1e5)}${NBSP}nC/cm²`, xpos, y0 + bandH + 11, { color: COL.text2, size: 9.5, align: side === 'L' ? 'left' : 'right' });
  };
  sigmaRef = Math.max(sigmaRef * 0.999, Math.abs(st.sigmaL), Math.abs(st.sigmaR), 1e-9);
  drawSigma(gx0, st.sigmaL, 'L');
  drawSigma(gx0 + wg, st.sigmaR, 'R');

  // x axis
  const ay = y0 + bandH + 24;
  ctx.strokeStyle = COL.axis; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(gx0, ay); ctx.lineTo(gx0 + wg, ay); ctx.stroke();
  const nTicks = 5;
  for (let i = 0; i <= nTicks; i++) {
    const xx = gx0 + (i / nTicks) * wg;
    ctx.beginPath(); ctx.moveTo(xx, ay); ctx.lineTo(xx, ay + 3); ctx.stroke();
    text(ctx, ((i / nTicks) * p.gapMM).toFixed(2), xx, ay + 6, { size: 9, align: 'center', base: 'top' });
  }
  text(ctx, 'x, mm', gx0 + wg + 4, ay + 6, { size: 9, base: 'top' });
  text(ctx, 'metal', x0 + metalW / 2, y0 - 3, { size: 8.5, align: 'center', base: 'bottom' });
  text(ctx, 'dielectric', x0 + metalW + wd1 / 2, y0 - 3, { size: 8.5, align: 'center', base: 'bottom' });
  text(ctx, 'gas gap', gx0 + wg / 2, y0 - 3, { size: 8.5, align: 'center', base: 'bottom', color: COL.text2 });

  // colorbar (a colormap without one is not quantitative)
  const cbx = W - padR + 14, cbw = 12, cbh = bandH;
  for (let i = 0; i < cbh; i++) {
    const t = 1 - i / (cbh - 1);
    ctx.fillStyle = lutCss(lut, Math.round(t * 255));
    ctx.fillRect(cbx, y0 + i, cbw, 1);
  }
  ctx.strokeStyle = COL.border; ctx.strokeRect(cbx, y0, cbw, cbh);
  const labTop = f.signed ? `+${sci(Math.max(Math.abs(lo), Math.abs(hi)))}` : (useLog ? '10' + supStr(Math.round(Math.log10(Math.max(hi, 1e-30)))) : sci(hi));
  const labBot = f.signed ? `−${sci(Math.max(Math.abs(lo), Math.abs(hi)))}` : (useLog ? '10' + supStr(Math.round(Math.log10(Math.max(lo, 1e-30)))) : sci(lo));
  text(ctx, labTop, cbx + cbw + 3, y0 + 4, { size: 9 });
  text(ctx, labBot, cbx + cbw + 3, y0 + cbh, { size: 9 });
  text(ctx, `${f.name}, ${f.unit}`, cbx - 2, y0 - 4, { size: 9, align: 'right', base: 'bottom', color: COL.text2 });

  // interactive readout
  gapGeom = { gx0, wg, y0, bandH, nx };
  if (gapHover >= 0) {
    const i = gidx(gapHover);
    const raw = f.abs ? Math.abs(arr[i]) : arr[i];
    // x measured from the gas-facing surface of the left barrier, not from the metal
    const xmm = (st.x[i] - GAS.x0) * 1e3;
    $('gapReadout').textContent = `${f.name}(x=${xmm.toFixed(3)} mm) = ${sci(raw)} ${f.unit}`;
  } else $('gapReadout').textContent = '';
  $('glowBadge').hidden = !(V.field === 'sion');
}
let sigmaRef = 1e-9;
let gapGeom = null, gapHover = -1;

let frontHist = { t: 0, x: 0, v: 0, age: 0 };
function frontSpeedText() {
  const st = solver.state;
  let mi = GAS.i0, mv = -1;
  for (let i = GAS.i0; i <= GAS.i1; i++) if (st.ionizRate[i] > mv) { mv = st.ionizRate[i]; mi = i; }
  const x = st.x[mi], t = st.t;
  // only meaningful if consecutive frames are close enough in sim time to be on the
  // same ionisation front; otherwise report n/a rather than a fabricated number
  const dtS = t - frontHist.t;
  if (dtS > 0 && dtS < 2e-7 && Math.abs(x - frontHist.x) > 0) {
    const v = Math.abs(x - frontHist.x) / dtS;
    if (isFinite(v) && v > 0) frontHist.v = 0.7 * frontHist.v + 0.3 * v;
    frontHist.age = 0;
  } else if (++frontHist.age > 120) frontHist.v = 0;
  frontHist.t = t; frontHist.x = x;
  return frontHist.v > 0 ? `${sci(frontHist.v)} m/s` : 'n/a (slow down playback)';
}

// ═══════════════════════════════════════════════════════════ history scan ═══
// min/max/last decimation over the visible time window. UI_SPEC §4.4:
// averaging or subsampling destroys the ns current peak — min/max envelope only.
function scanHistory(t0, t1, W, wantUg) {
  const h = solver.history;
  const n = Math.min(W, MAXW);
  buf.cnt.fill(0, 0, n);
  for (let i = 0; i < n; i++) {
    buf.min[i] = Infinity; buf.max[i] = -Infinity;
    buf.min2[i] = Infinity; buf.max2[i] = -Infinity;
    buf.min3[i] = Infinity; buf.max3[i] = -Infinity;
  }
  const span = t1 - t0;
  if (!(span > 0) || h.len === 0) return { samples: 0, bins: n };
  let samples = 0;
  // I_disch = I_total − C_cell·dU_app/dt  (Manley: subtract the cold-cell displacement current)
  const Cc = (solver.periodStats && solver.periodStats.Ccell) || 0;
  const start = (h.head - h.len + h.capacity) % h.capacity;
  let tPrev = 0, uPrev = 0, havePrev = false;
  for (let s = 0; s < h.len; s++) {
    const i = (start + s) % h.capacity;
    const t = h.t[i], u = h.Uapp[i], cur = h.current[i];
    // prefer the solver's own conduction current; fall back to I_total − C_cell·dU/dt
    let disch;
    if (h.condCurrent) disch = h.condCurrent[i];
    else {
      disch = cur;
      if (havePrev && t > tPrev) disch = cur - Cc * ((u - uPrev) / (t - tPrev));
    }
    tPrev = t; uPrev = u; havePrev = true;
    if (t < t0 || t > t1) continue;
    const b = Math.min(n - 1, Math.max(0, Math.floor(((t - t0) / span) * n)));
    if (disch < buf.min[b]) buf.min[b] = disch;
    if (disch > buf.max[b]) buf.max[b] = disch;
    buf.last[b] = disch;
    if (cur < buf.min3[b]) buf.min3[b] = cur;
    if (cur > buf.max3[b]) buf.max3[b] = cur;
    buf.last3[b] = cur;
    if (u < buf.min2[b]) buf.min2[b] = u;
    if (u > buf.max2[b]) buf.max2[b] = u;
    buf.last2[b] = u;
    if (wantUg) buf.ug[b] = h.Ugap[i];
    buf.cnt[b]++;
    samples++;
  }
  return { samples, bins: n };
}

/** Временной охват кольца history: за его пределами данных ФИЗИЧЕСКИ нет. */
function historySpan() {
  const h = solver.history;
  if (!h.len) return { t0: solver.state.t, t1: solver.state.t, span: 0 };
  const start = (h.head - h.len + h.capacity) % h.capacity;
  const last = (h.head - 1 + h.capacity) % h.capacity;
  return { t0: h.t[start], t1: h.t[last], span: h.t[last] - h.t[start] };
}
/** Позиционирование курсора с ограничением по реально доступной истории. */
function setCursor(t) {
  const hs = historySpan();
  V.cursorT = Math.max(hs.t0, Math.min(solver.state.t, t));
  V.followCursor = false;
}

// ═══════════════════════════════════════════════════════════ D. waveforms ═══
function drawWaves() {
  const o = cvWave; fit(o);
  const ctx = o.ctx; clearBg(o);
  const W = o.w, H = o.h;
  const padL = 46, padR = 10, padT = 8, padB = 22;
  const pw = Math.max(10, W - padL - padR);
  const gap = 8;
  const hU = Math.max(30, (H - padT - padB - gap) * 0.45);
  const hI = Math.max(30, H - padT - padB - gap - hU);
  const yU = padT, yI = padT + hU + gap;

  const lens = Math.pow(10, V.lensLog);
  // Пока идёт счёт, курсор приклеен к текущему времени: центрированное окно тогда
  // наполовину пусто. Прижимаем окно к правому краю (5% запаса), при ручном
  // позиционировании курсора оставляем его в центре.
  const t1 = V.followCursor ? V.cursorT + lens * 0.05 : V.cursorT + lens / 2;
  const t0 = t1 - lens;
  const res = scanHistory(t0, t1, Math.round(pw), true);

  const frame = (y, h, label) => {
    ctx.strokeStyle = COL.border; ctx.lineWidth = 1;
    ctx.strokeRect(padL + .5, y + .5, pw, h);
    text(ctx, label, padL + 3, y + 3, { size: 9, base: 'top', color: COL.muted });
  };

  // ── U track ──
  const U0 = solver.params.U0kV * 1e3 * 1.25;
  frame(yU, hU, 'kV');
  ctx.save(); ctx.beginPath(); ctx.rect(padL, yU, pw, hU); ctx.clip();
  const yOfU = (v) => yU + hU / 2 - (v / U0) * (hU / 2 - 2);
  ctx.strokeStyle = COL.grid;
  ctx.beginPath(); ctx.moveTo(padL, yOfU(0)); ctx.lineTo(padL + pw, yOfU(0)); ctx.stroke();
  const drawTrack = (getter, color, width, dash) => {
    ctx.setLineDash(dash); ctx.strokeStyle = color; ctx.lineWidth = width;
    ctx.beginPath(); let started = false;
    for (let i = 0; i < res.bins; i++) {
      if (!buf.cnt[i]) continue;
      const x = padL + i + 0.5, y = yOfU(getter(i));
      started ? ctx.lineTo(x, y) : (ctx.moveTo(x, y), started = true);
    }
    ctx.stroke(); ctx.setLineDash([]);
  };
  drawTrack((i) => buf.last2[i], '#e8eaed', 2, []);
  drawTrack((i) => buf.ug[i], '#a9b0bb', 2, [6, 3]);
  ctx.restore();
  for (const v of [-U0 / 1.25, 0, U0 / 1.25]) {
    text(ctx, (v / 1e3).toFixed(0), padL - 4, yOfU(v), { size: 9, align: 'right', base: 'middle' });
  }
  text(ctx, 'U_app', padL + pw - 4, yU + 3, { size: 9, align: 'right', base: 'top', color: '#e8eaed' });
  text(ctx, 'U_gap', padL + pw - 4, yU + 14, { size: 9, align: 'right', base: 'top', color: '#a9b0bb' });

  // ── I track ──
  frame(yI, hI, V.iLog ? 'A (log)' : 'mA');
  const st = solver.state;
  const iPk = Math.max(1e-9, st.peakCurrent || 1e-6);
  let yOfI;
  if (V.iLog) {
    const hiE = Math.ceil(Math.log10(iPk * 2)), loE = hiE - 6;
    yOfI = (v) => {
      const a = Math.log10(Math.max(Math.abs(v), 1e-30));
      return yI + hI - ((a - loE) / (hiE - loE)) * (hI - 2) - 1;
    };
    for (let e = loE; e <= hiE; e++) {
      const y = yOfI(Math.pow(10, e));
      if (y < yI || y > yI + hI) continue;
      ctx.strokeStyle = COL.grid; ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(padL + pw, y); ctx.stroke();
      text(ctx, '10' + supStr(e), padL - 4, y, { size: 9, align: 'right', base: 'middle' });
    }
  } else {
    const m = iPk * 1.2;
    yOfI = (v) => yI + hI / 2 - (v / m) * (hI / 2 - 2);
    for (const v of [-m / 1.2, 0, m / 1.2]) text(ctx, (v * 1e3).toFixed(1), padL - 4, yOfI(v), { size: 9, align: 'right', base: 'middle' });
    ctx.strokeStyle = COL.grid; ctx.beginPath(); ctx.moveTo(padL, yOfI(0)); ctx.lineTo(padL + pw, yOfI(0)); ctx.stroke();
  }
  ctx.save(); ctx.beginPath(); ctx.rect(padL, yI, pw, hI); ctx.clip();
  // min/max envelope
  ctx.strokeStyle = 'rgba(255,209,102,0.45)'; ctx.lineWidth = 1;
  ctx.beginPath();
  for (let i = 0; i < res.bins; i++) {
    if (!buf.cnt[i] || buf.cnt[i] < 2) continue;
    const x = padL + i + 0.5;
    ctx.moveTo(x, yOfI(buf.min[i])); ctx.lineTo(x, yOfI(buf.max[i]));
  }
  ctx.stroke();
  ctx.strokeStyle = '#ffd166'; ctx.lineWidth = 2;
  ctx.beginPath(); let started = false;
  for (let i = 0; i < res.bins; i++) {
    if (!buf.cnt[i]) continue;
    const x = padL + i + 0.5, y = yOfI(buf.last[i]);
    started ? ctx.lineTo(x, y) : (ctx.moveTo(x, y), started = true);
  }
  ctx.stroke();
  if (V.showItot) { // thin dashed I_total (UI_SPEC §3.2 auxiliary series)
    ctx.setLineDash([3, 3]); ctx.strokeStyle = '#ffd166'; ctx.lineWidth = 1;
    ctx.beginPath(); let st2 = false;
    for (let i = 0; i < res.bins; i++) {
      if (!buf.cnt[i]) continue;
      const x = padL + i + 0.5, y = yOfI(buf.last3[i]);
      st2 ? ctx.lineTo(x, y) : (ctx.moveTo(x, y), st2 = true);
    }
    ctx.stroke(); ctx.setLineDash([]);
  }
  // breakdown event markers
  ctx.strokeStyle = 'rgba(255,209,102,.55)'; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
  for (let k = 0; k < buf.evN; k++) {
    const idx = (buf.evHead - buf.evN + k + buf.evT.length) % buf.evT.length;
    const te = buf.evT[idx];
    if (te < t0 || te > t1) continue;
    const x = padL + ((te - t0) / (t1 - t0)) * pw;
    ctx.beginPath(); ctx.moveTo(x, yI); ctx.lineTo(x, yI + hI); ctx.stroke();
  }
  ctx.setLineDash([]);
  ctx.restore();
  text(ctx, V.showItot ? 'I_disch (bold) · I_total' : 'I_disch', padL + pw - 4, yI + 3, { size: 9, align: 'right', base: 'top', color: '#ffd166' });

  // time axis + cursor
  ctx.strokeStyle = COL.axis; ctx.beginPath();
  ctx.moveTo(padL, yI + hI + .5); ctx.lineTo(padL + pw, yI + hI + .5); ctx.stroke();
  for (let i = 0; i <= 4; i++) {
    const x = padL + (i / 4) * pw, tt = t0 + (i / 4) * (t1 - t0);
    text(ctx, fmtTime(tt), x, yI + hI + 5, { size: 9, align: i === 0 ? 'left' : i === 4 ? 'right' : 'center', base: 'top' });
  }
  const xc = padL + ((V.cursorT - t0) / (t1 - t0)) * pw;
  if (xc >= padL && xc <= padL + pw) {
    ctx.strokeStyle = COL.cursor; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(xc + .5, padT); ctx.lineTo(xc + .5, yI + hI); ctx.stroke();
  }
  if (res.samples === 0) {
    const hs = historySpan();
    text(ctx, hs.span > 0
      ? `no samples in this window — history covers the last ${fmtTime(hs.span)}`
      : 'no samples yet — press ⏵ play',
      padL + pw / 2, yU + hU / 2, { size: 11.5, align: 'center', base: 'middle', color: COL.muted, mono: false });
  }
  $('waveCursor').textContent = `t = ${fmtTime(V.cursorT, 4)}`;
  const dec = res.samples > res.bins ? Math.round(res.samples / res.bins) : 0;
  const badge = $('decBadge');
  badge.hidden = dec < 2;
  if (dec >= 2) badge.textContent = `decimated ×${dec} (min–max envelope)`;
}
function fmtTime(t, digits = 3) {
  const a = Math.abs(t);
  if (a >= 1e-3) return `${(t * 1e3).toFixed(digits)} ms`;
  if (a >= 1e-6) return `${(t * 1e6).toFixed(digits)} µs`;
  return `${(t * 1e9).toFixed(1)} ns`;
}

// ═══════════════════════════════════════════════════════════ E. Lissajous ═══
function drawQV() {
  const o = cvQV; fit(o);
  const ctx = o.ctx; clearBg(o);
  const W = o.w, H = o.h;
  const pad = 34;
  const side = Math.max(20, Math.min(W - pad - 12, H - pad - 14));
  const px0 = pad, py0 = 6, pw = side, ph = side;
  const h = solver.history;
  const T = 1 / (solver.params.freqKHz * 1e3);
  const tEnd = solver.state.t;
  const tStart = tEnd - T * (V.qvWindow === 8 ? 8 : 1);

  let umax = solver.params.U0kV * 1e3 * 1.1, qmax = 1e-12;
  const start = (h.head - h.len + h.capacity) % h.capacity;
  for (let s = 0; s < h.len; s++) {
    const i = (start + s) % h.capacity;
    if (h.t[i] < tStart) continue;
    const q = Math.abs(h.charge[i]); if (q > qmax) qmax = q;
  }
  qmax *= 1.15;
  const X = (u) => px0 + ((u / umax + 1) / 2) * pw;
  const Y = (q) => py0 + ph - ((q / qmax + 1) / 2) * ph;

  // frame + grid
  ctx.strokeStyle = COL.border; ctx.strokeRect(px0 + .5, py0 + .5, pw, ph);
  ctx.strokeStyle = COL.grid;
  ctx.beginPath(); ctx.moveTo(X(0), py0); ctx.lineTo(X(0), py0 + ph);
  ctx.moveTo(px0, Y(0)); ctx.lineTo(px0 + pw, Y(0)); ctx.stroke();
  for (let i = 0; i <= 4; i++) {
    const u = -umax + (i / 4) * 2 * umax;
    text(ctx, (u / 1e3).toFixed(0), X(u), py0 + ph + 4, { size: 9, align: 'center', base: 'top' });
    const q = -qmax + (i / 4) * 2 * qmax;
    text(ctx, (q * 1e9).toFixed(1), px0 - 4, Y(q), { size: 9, align: 'right', base: 'middle' });
  }
  text(ctx, 'U_app, kV', px0 + pw, py0 + ph + 15, { size: 9, align: 'right', base: 'top', color: COL.text2 });
  text(ctx, 'Q, nC', px0 - 4, py0 - 2, { size: 9, align: 'right', base: 'bottom', color: COL.text2 });

  // previous periods (ghosts) then current period
  const drawSeg = (ta, tb, color, width, alpha) => {
    ctx.globalAlpha = alpha; ctx.strokeStyle = color; ctx.lineWidth = width;
    ctx.beginPath(); let started = false;
    for (let s = 0; s < h.len; s++) {
      const i = (start + s) % h.capacity;
      if (h.t[i] < ta || h.t[i] > tb) { started = false; continue; }
      const x = X(h.Uapp[i]), y = Y(h.charge[i]);
      started ? ctx.lineTo(x, y) : (ctx.moveTo(x, y), started = true);
    }
    ctx.stroke(); ctx.globalAlpha = 1;
  };
  // «Призраки» прошлых периодов берутся из solver.qvHistory, а НЕ из history-кольца:
  // кольцо (4096 отсчётов при записи раз в T/2000) вмещает всего ~2 периода, поэтому
  // кнопка «last 8» рисовала 2 периода и создавала ложное впечатление идеальной
  // повторяемости разряда.
  const gh = solver.qvHistory;
  let ghostsShown = 0;
  if (V.qvWindow === 8 && gh && gh.n > 1) {
    const nShow = Math.min(gh.n - 1, 7);
    for (let k = nShow; k >= 1; k--) {
      const idx = (gh.head - 1 - k + gh.cap * 2) % gh.cap;
      const u = gh.u[idx], q = gh.q[idx];
      ctx.globalAlpha = 0.12 + 0.05 * (nShow - k);
      ctx.strokeStyle = '#e8eaed'; ctx.lineWidth = 1;
      ctx.beginPath();
      for (let j = 0; j < u.length; j++) {
        const x = X(u[j]), y = Y(q[j]);
        j ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      }
      ctx.closePath(); ctx.stroke(); ctx.globalAlpha = 1;
      ghostsShown++;
    }
  }
  drawSeg(tEnd - T, tEnd, '#e8eaed', 2, 1);

  // Manley slopes drawn through the centroid
  const ps = solver.periodStats || {};
  const Cd = ps.Cdiel || 0, Cc = ps.Ccell || 0;
  ctx.setLineDash([4, 4]); ctx.lineWidth = 1; ctx.strokeStyle = COL.muted;
  const line = (slope, off) => {
    ctx.beginPath();
    ctx.moveTo(X(-umax), Y(-umax * slope + off));
    ctx.lineTo(X(umax), Y(umax * slope + off));
    ctx.stroke();
  };
  if (Cd > 0) { line(Cd, qmax * 0.35); line(Cd, -qmax * 0.35); }
  if (Cc > 0) { line(Cc, 0); }
  ctx.setLineDash([]);

  // marker for the live point
  const stq = solver.state;
  ctx.fillStyle = COL.cursor;
  ctx.beginPath(); ctx.arc(X(stq.Uapp), Y(stq.charge), 3, 0, 7); ctx.fill();

  const E = ps.energyPerPeriodJ || 0;
  $('qvFoot').innerHTML =
    `C<sub>diel</sub> ${sig3(Cd * 1e12)}${NBSP}pF <span class="muted">(geom ${sig3((ps.CdielGeom || 0) * 1e12)})</span> · ` +
    `C<sub>cell</sub> ${sig3(Cc * 1e12)}${NBSP}pF · ` +
    `E ${sig3(E * 1e3)}${NBSP}mJ · P ${sig3(ps.powerW || 0)}${NBSP}W · U_burn ${sig3(ps.UburnkV || 0)}${NBSP}kV` +
    (V.qvWindow === 8 ? ` · ghosts ${ghostsShown}` : '');
  // Прежняя проверка (C_diel > 1.05*C_cell) не могла сработать НИКОГДА: при вырождении
  // фита обе величины падают на геометрические константы с фиксированным отношением.
  // Судим по качеству фита сторон: R² и число найденных дуг (UI_SPEC §2.3).
  const warn = $('qvWarn');
  const bad = !(ps.qvOk === true);
  warn.hidden = !bad;
  if (bad) {
    warn.textContent = (ps.qvFallback !== false || !(ps.qvArcsOn > 0))
      ? `⚠ sides not resolved (${ps.qvArcsOn || 0}/${ps.qvArcsOff || 0}) · geometric C shown`
      : `⚠ not a parallelogram · R² ${sig3(ps.qvR2on || 0)}/${sig3(ps.qvR2off || 0)} · arcs ${ps.qvArcsOn}/${ps.qvArcsOff}`;
    warn.title = 'Manley slopes are taken from least-squares fits of the burning / dark arcs. '
      + 'R² < 0.98 means the side is not straight, i.e. the figure is not a parallelogram '
      + 'and C_diel is a fitted slope, not a parallelogram slope.';
  }
}

// ═══════════════════════════════════════════════════════════ F. profiles ════
const tsBuf = {};
for (const s of SPECIES) tsBuf[s.key] = { a: new Float64Array(1024), n: 0, h: 0 };
const tsT = { a: new Float64Array(1024), n: 0, h: 0 };
// Раньше временной ряд писался ТОЛЬКО на границе периода. Реальная скорость счёта
// ~1e-6 сим-с за секунду wall => одна точка раз в 40 с, и панель «time-series»
// (режим по умолчанию у пресета ozonizer) была пуста первые полминуты и деградирована
// минутами. Теперь — по стенным часам, с троттлом 100 мс.
let tsLastWall = 0;
function maybePushTimeSeries(nowMs) {
  if (nowMs - tsLastWall < 100) return;
  tsLastWall = nowMs;
  pushTimeSeries();
}
function pushTimeSeries() {
  const st = solver.state;
  for (const s of SPECIES) {
    const arr = st.n[s.key];
    let m = 0;
    // average over GAS cells only — dielectric cells hold no species
    if (arr) { for (let i = GAS.i0; i <= GAS.i1; i++) m += arr[i]; m /= GAS.n; }
    const b = tsBuf[s.key];
    b.a[b.h] = m; b.h = (b.h + 1) % 1024; if (b.n < 1024) b.n++;
  }
  tsT.a[tsT.h] = st.t; tsT.h = (tsT.h + 1) % 1024; if (tsT.n < 1024) tsT.n++;
}

let profGeom = null, profHover = -1;
function drawProfiles() {
  const o = cvProf; fit(o);
  const ctx = o.ctx; clearBg(o);
  const st = solver.state;
  const W = o.w, H = o.h;
  const padL = 44, padR = 62, padT = 8, padB = 20;
  const pw = Math.max(10, W - padL - padR), ph = Math.max(10, H - padT - padB);

  // y range (log10 m^-3, never narrower than 6 decades)
  let hi = 12;
  const vis = SPECIES.filter((s) => visible(s.key) && st.n[s.key]);
  for (const s of vis) {
    const a = st.n[s.key];
    for (let i = GAS.i0; i <= GAS.i1; i++) { const l = Math.log10(Math.max(a[i], 1e-30)); if (l > hi) hi = l; }
  }
  hi = Math.ceil(hi);
  const lo = Math.min(12, hi - 6);
  const Y = (v) => padT + ph - ((Math.log10(Math.max(v, 1e-30)) - lo) / (hi - lo)) * ph;

  ctx.strokeStyle = COL.border; ctx.strokeRect(padL + .5, padT + .5, pw, ph);
  for (let e = Math.ceil(lo); e <= hi; e++) {
    const y = Y(Math.pow(10, e));
    ctx.strokeStyle = COL.grid; ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(padL + pw, y); ctx.stroke();
    text(ctx, '10' + supStr(e), padL - 4, y, { size: 9, align: 'right', base: 'middle' });
  }

  const isX = V.profMode === 'x';
  const gapMM = solver.params.gapMM;
  const nx = GAS.n;
  let X, count;
  // in x-mode the curve index is GAS-LOCAL (0..GAS.n-1); gpos() maps a global
  // cell index onto 0..1 across the gap regardless of where the gas starts
  if (isX) { count = nx; X = (k) => padL + gpos(st, gidx(k)) * pw; }
  else {
    count = tsT.n;
    if (count < 2) {
      text(ctx, 'collecting data…', padL + pw / 2, padT + ph / 2,
        { size: 11, align: 'center', base: 'middle', color: COL.muted, mono: false });
    }
    const t0 = tsT.n ? tsT.a[(tsT.h - tsT.n + 1024) % 1024] : 0;
    const t1 = st.t;
    X = (i) => padL + ((tsT.a[(tsT.h - tsT.n + i + 1024) % 1024] - t0) / Math.max(1e-12, t1 - t0)) * pw;
  }

  const labels = [];
  for (const s of vis) {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.w; ctx.setLineDash(s.dash);
    ctx.beginPath();
    let ylast = 0;
    for (let i = 0; i < count; i++) {
      const v = isX ? st.n[s.key][gidx(i)] : tsBuf[s.key].a[(tsBuf[s.key].h - tsBuf[s.key].n + i + 1024) % 1024];
      const x = X(i), y = Y(v);
      i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      ylast = y;
    }
    ctx.stroke(); ctx.setLineDash([]);
    if (count) labels.push({ y: ylast, s });
  }
  // direct end-of-curve labels (never colour alone), de-overlapped
  labels.sort((a, b) => a.y - b.y);
  let prev = -Infinity;
  for (const l of labels) {
    const y = Math.max(prev + 11, Math.min(padT + ph, l.y));
    prev = y;
    text(ctx, l.s.label, padL + pw + 4, y, { size: 10, base: 'middle', color: l.s.color });
  }

  ctx.strokeStyle = COL.axis; ctx.beginPath();
  ctx.moveTo(padL, padT + ph + .5); ctx.lineTo(padL + pw, padT + ph + .5); ctx.stroke();
  for (let i = 0; i <= 4; i++) {
    const x = padL + (i / 4) * pw;
    let lab;
    if (isX) lab = ((i / 4) * gapMM).toFixed(2);
    else if (tsT.n < 2) lab = '–';
    else {
      const tt = tsT.a[(tsT.h - tsT.n + Math.floor((i / 4) * (tsT.n - 1)) + 1024) % 1024];
      const spanT = st.t - tsT.a[(tsT.h - tsT.n + 1024) % 1024];
      lab = fmtTime(tt, spanT > 0 && spanT < 1e-5 ? 4 : 3);
    }
    text(ctx, lab, x, padT + ph + 4, { size: 9, align: 'center', base: 'top' });
  }
  text(ctx, isX ? 'x, mm' : 't', padL + pw + 4, padT + ph + 4, { size: 9, base: 'top', color: COL.text2 });

  profGeom = { padL, padT, pw, ph, lo, hi, isX, count };
  if (profHover >= 0 && isX) {
    const k = Math.min(nx - 1, Math.max(0, profHover));
    const i = gidx(k);
    const x = X(k);
    ctx.strokeStyle = COL.cursor; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + ph); ctx.stroke(); ctx.setLineDash([]);
    const rows = SPECIES.filter((s) => st.n[s.key]).map((s) => ({ s, v: st.n[s.key][i] })).sort((a, b) => b.v - a.v);
    let qsum = 0;
    for (const r of rows) qsum += CHARGE[r.s.key] * r.v;
    const tip = $('profTip');
    tip.hidden = false;
    tip.textContent = `x = ${((st.x[i] - GAS.x0) * 1e3).toFixed(3)} mm\n` +
      rows.map((r) => `${r.s.label.padEnd(7)} ${sci(r.v)}`).join('\n') +
      `\nΣq/e = ${sci(qsum)}`;
    const bx = Math.min(o.w - 150, x + 10);
    tip.style.left = bx + 'px'; tip.style.top = padT + 4 + 'px';
    $('profCursor').textContent = `x = ${((st.x[i] - GAS.x0) * 1e3).toFixed(3)} mm`;
  } else { $('profTip').hidden = true; $('profCursor').textContent = 'x = –'; }
}

// ═══════════════════════════════════════════════════════════ H. timelines ═══
function drawMacro() {
  const o = cvMacro; fit(o);
  const ctx = o.ctx; clearBg(o, COL.surface2);
  const W = o.w, H = o.h;
  const t1 = Math.max(solver.state.t, 1e-9), t0 = 0;
  const T = 1 / (solver.params.freqKHz * 1e3);
  const X = (t) => ((t - t0) / (t1 - t0)) * W;
  // period ticks
  ctx.strokeStyle = COL.grid; ctx.lineWidth = 1;
  const nP = t1 / T;
  const stepP = Math.max(1, Math.pow(10, Math.floor(Math.log10(Math.max(1, nP / 10)))));
  ctx.beginPath();
  for (let k = 0; k * stepP * T <= t1; k++) { const x = X(k * stepP * T); ctx.moveTo(x, H - 6); ctx.lineTo(x, H); }
  ctx.stroke();
  // breakdown events: ▲ up for positive half-period, ▼ down for negative
  for (let k = 0; k < buf.evN; k++) {
    const idx = (buf.evHead - buf.evN + k + buf.evT.length) % buf.evT.length;
    const x = X(buf.evT[idx]);
    const up = buf.evP[idx] > 0;
    ctx.fillStyle = up ? '#e66767' : '#3987e5';
    ctx.beginPath();
    if (up) { ctx.moveTo(x, 3); ctx.lineTo(x - 3, 9); ctx.lineTo(x + 3, 9); }
    else { ctx.moveTo(x, H - 3); ctx.lineTo(x - 3, H - 9); ctx.lineTo(x + 3, H - 9); }
    ctx.closePath(); ctx.fill();
  }
  // область без данных (за пределами кольца history) — заштрихована, туда курсор
  // не ставится: раньше туда можно было перетащить курсор и получить пустые панели
  const hs = historySpan();
  if (hs.span > 0 && hs.t0 > t0) {
    ctx.fillStyle = 'rgba(0,0,0,.45)';
    ctx.fillRect(0, 0, X(hs.t0), H);
    ctx.strokeStyle = COL.border; ctx.beginPath();
    ctx.moveTo(X(hs.t0) + .5, 0); ctx.lineTo(X(hs.t0) + .5, H); ctx.stroke();
    if (X(hs.t0) > 60) text(ctx, 'no history', 4, 9, { size: 8.5, base: 'top' });
  }
  // lens window (min 4 px — otherwise invisible, the classic brush bug)
  const lens = Math.pow(10, V.lensLog);
  const lx0 = X(V.cursorT - lens / 2), lw = Math.max(4, X(V.cursorT + lens / 2) - lx0);
  ctx.fillStyle = 'rgba(77,163,255,.18)'; ctx.fillRect(lx0, 0, lw, H);
  ctx.strokeStyle = COL.accent; ctx.strokeRect(lx0 + .5, .5, lw, H - 1);
  ctx.strokeStyle = COL.cursor; ctx.beginPath();
  ctx.moveTo(X(V.cursorT) + .5, 0); ctx.lineTo(X(V.cursorT) + .5, H); ctx.stroke();
  text(ctx, `0`, 3, H - 3, { size: 9 });
  text(ctx, `${fmtTime(t1)} (${Math.floor(nP)} T)`, W - 3, H - 3, { size: 9, align: 'right' });
  macroGeom = { t0, t1, W };
}
let macroGeom = null;

function drawMicro() {
  const o = cvMicro; fit(o);
  const ctx = o.ctx; clearBg(o, COL.surface2);
  const W = o.w, H = o.h;
  const lens = Math.pow(10, V.lensLog);
  const t0 = V.cursorT - lens / 2, t1 = V.cursorT + lens / 2;
  const res = scanHistory(t0, t1, W, false);
  const iPk = Math.max(1e-9, solver.state.peakCurrent || 1e-6);
  ctx.strokeStyle = '#ffd166'; ctx.lineWidth = 1;
  ctx.beginPath();
  for (let i = 0; i < res.bins; i++) {
    if (!buf.cnt[i]) continue;
    const a = Math.min(1, Math.abs(buf.max[i]) / iPk);
    const b = Math.min(1, Math.abs(buf.min[i]) / iPk);
    const m = Math.max(a, b);
    ctx.moveTo(i + .5, H); ctx.lineTo(i + .5, H - m * (H - 4));
  }
  ctx.stroke();
  ctx.strokeStyle = COL.cursor; ctx.beginPath();
  ctx.moveTo(W / 2 + .5, 0); ctx.lineTo(W / 2 + .5, H); ctx.stroke();
  text(ctx, `−${fmtTime(lens / 2)}`, 3, H - 3, { size: 9 });
  text(ctx, `+${fmtTime(lens / 2)}`, W - 3, H - 3, { size: 9, align: 'right' });
  microGeom = { t0, t1, W };
}
let microGeom = null;

// ═══════════════════════════════════════════════════════ event detection ════
// Breakdown = rising crossing of 0.1·I_peak with a 5 ns dead time (UI_SPEC §8.2).
let evScanT = 0, evPrevI = 0;
function scanEvents() {
  const h = solver.history;
  if (!h.len) return;
  const thr = 0.1 * Math.max(1e-12, solver.state.peakCurrent || 0);
  const Cc = (solver.periodStats && solver.periodStats.Ccell) || 0;
  const start = (h.head - h.len + h.capacity) % h.capacity;
  let tP = 0, uP = 0, has = false;
  for (let s = 0; s < h.len; s++) {
    const i = (start + s) % h.capacity;
    const t = h.t[i];
    // discharge current only — the displacement current must not trigger events.
    // Prefer the solver's own conduction current when the ring buffer carries it;
    // otherwise fall back to subtracting C_cell*dU/dt numerically.
    let disch;
    if (h.condCurrent) disch = h.condCurrent[i];
    else {
      disch = h.current[i];
      if (has && t > tP) disch -= Cc * ((h.Uapp[i] - uP) / (t - tP));
    }
    tP = t; uP = h.Uapp[i]; has = true;
    if (t <= evScanT) continue;
    const a = Math.abs(disch);
    if (a > thr && evPrevI <= thr) {
      const last = buf.evN ? buf.evT[(buf.evHead - 1 + buf.evT.length) % buf.evT.length] : -1;
      if (t - last > 5e-9) {
        buf.evT[buf.evHead] = t;
        buf.evP[buf.evHead] = h.Uapp[i] >= 0 ? 1 : -1;
        buf.evHead = (buf.evHead + 1) % buf.evT.length;
        if (buf.evN < buf.evT.length) buf.evN++;
        announce(`breakdown at ${fmtTime(t)}`);
      }
    }
    evPrevI = a;
    evScanT = t;
  }
}
let liveTimer = 0;
function announce(msg) {
  const now = performance.now();
  if (now - liveTimer < 1500) return;
  liveTimer = now;
  $('live').textContent = msg;
}

// ═══════════════════════════════════════════════════════════ metrics ════════
let lastMetricT = 0, lastPeriodIdx = -1, cumEnergyJ = 0;
function updateMetrics(force) {
  const now = performance.now();
  if (!force && now - lastMetricT < 250) return;   // 4 Hz — faster is unreadable
  lastMetricT = now;
  const st = solver.state, ps = solver.periodStats || {};
  const p = solver.params;
  const volumeL = (p.areaCM2 * p.gapMM * 0.1) / 1000; // cm^3 → L
  const set = (k, v, unit) => {
    if (!tileEls[k]) return;
    tileEls[k].n.textContent = v;
    if (unit !== undefined && tileEls[k].u) tileEls[k].u.textContent = unit;
  };
  set('power', sig3(ps.powerW ?? st.power ?? 0));
  // Specific energy = cumulative deposited energy per litre of gas in the gap
  set('spec', sig3(cumEnergyJ / Math.max(1e-12, volumeL)));
  const o3 = withUnit(st.o3ppm || 0, 'ppm', [[1e4, '%']]);
  set('o3', o3.unit === '%' ? sig3(o3.v / 100) : sig3(o3.v), o3.unit);
  set('en', sig3(st.maxEN || 0));
  const ipk = withUnit((st.peakCurrent || 0) * 1e3, 'mA', [[1e3, 'A']]);
  set('ipk', sig3(ipk.v), ipk.unit);
  set('bd', String(st.breakdownsPerPeriod ?? 0));
  set('cdiel', sig3((ps.Cdiel || 0) * 1e12));
  set('ccell', sig3((ps.Ccell || 0) * 1e12));
  set('eper', sig3((ps.energyPerPeriodJ || 0) * 1e3));

  const T = 1 / (p.freqKHz * 1e3);
  const pIdx = Math.floor(st.t / T);
  if (pIdx !== lastPeriodIdx) {
    lastPeriodIdx = pIdx;
    sparkPush('power', ps.powerW || 0); sparkPush('o3', st.o3ppm || 0);
    sparkPush('en', st.maxEN || 0); sparkPush('ipk', st.peakCurrent || 0);
    sparkPush('bd', st.breakdownsPerPeriod || 0); sparkPush('cdiel', ps.Cdiel || 0);
    sparkPush('ccell', ps.Ccell || 0); sparkPush('eper', ps.energyPerPeriodJ || 0);
    cumEnergyJ += ps.energyPerPeriodJ || 0;
    sparkPush('spec', cumEnergyJ);
    for (const [k] of METRICS) drawSpark(tileEls[k], k);
    $('panelMetrics').setAttribute('aria-label',
      `period ${pIdx}: ${st.breakdownsPerPeriod} breakdowns, peak current ${sig3((st.peakCurrent || 0) * 1e3)} mA, power ${sig3(ps.powerW || 0)} W`);
    if (V.steady) checkSteady(ps);
  }
  updateGridInfo();   // 4 Гц, по ТЕКУЩЕМУ n_e (а не раз в период, в тихой фазе)
  $('periodText').textContent = `period #${pIdx} · t = ${fmtTime(st.t)}`;
}

let steadyHist = [], steadyCount = 0;
function checkSteady(ps) {
  const e = ps.energyPerPeriodJ || 0;
  steadyHist.push(e);
  if (steadyHist.length > 2) steadyHist.shift();
  if (steadyHist.length === 2) {
    const rel = Math.abs(steadyHist[1] - steadyHist[0]) / Math.max(1e-15, Math.abs(steadyHist[1]));
    steadyCount = rel < 0.005 ? steadyCount + 1 : 0;
    if (steadyCount >= 5) {
      V.steady = false; setRunning(false);
      $('btnSteady').classList.remove('on');
      announce(`steady state reached (Δ ${(rel * 100).toFixed(2)} %)`);
    }
  }
}

function updateStatus(steps, simDt, wallMs) {
  const st = solver.state;
  const el = $('solverStatus');
  const res = st.residual ?? null;
  // доля отбракованных шагов — главный индикатор «схема ушла в нефизичный режим»
  const rejFrac = st.steps ? (st.rejects || 0) / st.steps : 0;
  let cls = COL.muted;
  if (rejFrac > 0.01 || (st.wallCfl || 0) > 1) cls = COL.warning;
  if (res !== null && res > 1e-8) cls = COL.critical;
  el.style.color = cls;
  el.textContent =
    `dt ${sci(st.dt)}${NBSP}s · steps/frame ${steps} · wall ${wallMs.toFixed(1)}${NBSP}ms → sim ${fmtTime(simDt)}` +
    ` · rejects ${sig3(rejFrac * 100)}%` +
    (res !== null ? ` · resid ${sci(res)}` : '') +
    ` · clip ${sci(st.clipCharge || 0)}${NBSP}C`;
}

// ═══════════════════════════════════════════════════════════ main loop ══════
let fpsT = performance.now(), fpsN = 0, fps = 0;
let lastFrameWall = performance.now();

function frame() {
  requestAnimationFrame(frame);
  const now = performance.now();
  const dtWall = Math.min(0.25, (now - lastFrameWall) / 1000);
  lastFrameWall = now;

  let steps = 0, simDt = 0, wallMs = 0;
  if (V.running) {
    const speed = Math.pow(10, V.speedLog);   // sim-seconds per wall-second
    const want = solver.state.t + speed * dtWall;
    solver.targetSimTime = V.stopAtT !== null ? Math.min(want, V.stopAtT) : want;
    const w0 = performance.now();
    const r = solver.advance(12);
    wallMs = performance.now() - w0;
    steps = r.steps; simDt = r.simTime;
    if (V.stopAtT !== null && solver.state.t >= V.stopAtT) { V.stopAtT = null; setRunning(false); }
    if (V.followCursor) V.cursorT = solver.state.t;
  }

  const st = solver.state;
  updateGasWindow(st);   // must precede every gap-space draw call
  updateGlow(st, dtWall);
  maybePushTimeSeries(now);
  sampleGridWorst(st);
  scanEvents();
  drawGap();
  drawWaves();
  drawQV();
  drawProfiles();
  drawMacro();
  drawMicro();
  updateMetrics(false);
  if (V.running) updateStatus(steps, simDt, wallMs);

  fpsN++;
  if (now - fpsT > 500) {
    fps = (fpsN * 1000) / (now - fpsT); fpsT = now; fpsN = 0;
    $('fpsText').textContent = `${fps.toFixed(0)} fps`;
  }
}

function refreshAll() { /* redraw happens next rAF; kept for explicit call sites */ }

// ═══════════════════════════════════════════════════════════ transport ══════
function setRunning(on) {
  V.running = on;
  $('btnPlay').textContent = on ? '⏸ pause' : '⏵ play';
  $('runBadge').classList.toggle('running', on);
  $('runText').textContent = on ? 'RUNNING' : 'PAUSED';
}
function updateSpeedLabel() {
  const s = Math.pow(10, V.speedLog);
  $('speedText').textContent = `⟳ ${sci(s)} s/s`;
  $('speedLabel').textContent = `${sci(s)} s/s  (1 s → ${fmtTime(s)})`;
}
function updateLensLabel() {
  const w = Math.pow(10, V.lensLog);
  const T = 1 / (solver.params.freqKHz * 1e3);
  $('lensLabel').textContent = `lens ${fmtTime(w)} (${sig3(w / T)} T)`;
}
/** Поставить окно осциллографа равным k периодам текущей частоты. */
function setLensPeriods(k) {
  V.lensPeriods = k;
  const T = 1 / (solver.params.freqKHz * 1e3);
  const el = $('lens');
  const lo = +el.min, hi = +el.max;
  V.lensLog = Math.max(lo, Math.min(hi, Math.log10(k * T)));
  el.value = V.lensLog;
  updateLensLabel();
}

$('btnPlay').addEventListener('click', () => setRunning(!V.running));
$('btnReset').addEventListener('click', () => doReset());
$('btnResetSoft').addEventListener('click', () => doReset());
$('btnResetPreset').addEventListener('click', () => applyPreset(V.preset, true));
$('btnStep').addEventListener('click', () => { solver.targetSimTime = Infinity; solver.step(); V.cursorT = solver.state.t; updateMetrics(true); });
$('btnPeriod').addEventListener('click', () => {
  V.stopAtT = solver.state.t + 1 / (solver.params.freqKHz * 1e3);
  V.followCursor = true; setRunning(true);
});
$('btnSteady').addEventListener('click', (e) => {
  V.steady = !V.steady; steadyHist = []; steadyCount = 0;
  e.currentTarget.classList.toggle('on', V.steady);
  if (V.steady) { V.followCursor = true; setRunning(true); }
});
$('btnPrevEv').addEventListener('click', () => snapEvent(-1));
$('btnNextEv').addEventListener('click', () => snapEvent(+1));
function snapEvent(dir) {
  let best = null;
  for (let k = 0; k < buf.evN; k++) {
    const t = buf.evT[(buf.evHead - buf.evN + k + buf.evT.length) % buf.evT.length];
    if (dir < 0 && t < V.cursorT - 1e-12 && (best === null || t > best)) best = t;
    if (dir > 0 && t > V.cursorT + 1e-12 && (best === null || t < best)) best = t;
  }
  if (best !== null) { V.cursorT = best; V.followCursor = false; setRunning(false); }
  else announce('no further breakdown event in buffer');
}

$('speed').addEventListener('input', (e) => { V.speedLog = +e.target.value; updateSpeedLabel(); });
$('lens').addEventListener('input', (e) => {
  V.lensLog = +e.target.value;
  V.lensPeriods = Math.pow(10, V.lensLog) * solver.params.freqKHz * 1e3;
  updateLensLabel();
});
// быстрые зумы: один период / восемь периодов / масштаб филамента
for (const [id, k] of [['lens1T', 1], ['lens8T', 8], ['lensFil', 2e-8 * solver.params.freqKHz * 1e3]]) {
  const b = $(id);
  if (!b) continue;
  b.addEventListener('click', () => {
    if (id === 'lensFil') { V.lensLog = -7.3; V.lensPeriods = Math.pow(10, -7.3) * solver.params.freqKHz * 1e3; $('lens').value = V.lensLog; updateLensLabel(); }
    else setLensPeriods(k);
    for (const o of ['lens1T', 'lens8T', 'lensFil']) $(o).classList.toggle('on', o === id);
  });
}

function doReset() {
  solver.reset();
  resetUiState();
  V.steady = false;
  $('btnSteady').classList.remove('on');
  updateMetrics(true);
}

// ═══════════════════════════════════════════════════════════ presets ════════
function applyPreset(name, hard) {
  const pr = PRESETS[name]; if (!pr) return;
  V.preset = name;
  $('preset').value = name;
  solver.setParams(Object.assign({}, pr.p));
  V.field = pr.view.field; V.log = pr.view.log; V.iLog = pr.view.iLog;
  V.profMode = pr.view.profMode;
  setLensPeriods(pr.view.lensPeriods || 1);
  for (const o of ['lens1T', 'lens8T', 'lensFil']) {
    const b = $(o); if (b) b.classList.toggle('on', o === (pr.view.lensPeriods === 8 ? 'lens8T' : 'lens1T'));
  }
  $('gapLog').checked = V.log;
  for (const b of $('fieldSeg').children) b.classList.toggle('on', b.dataset.field === V.field);
  $('iLog').classList.toggle('on', V.iLog); $('iLin').classList.toggle('on', !V.iLog);
  $('profX').classList.toggle('on', V.profMode === 'x'); $('profT').classList.toggle('on', V.profMode !== 'x');
  syncAllControls();
  onParamsChanged(true);
  if (hard !== false) doReset();
}
$('preset').addEventListener('change', (e) => applyPreset(e.target.value, true));

// ═══════════════════════════════════════════════════════════ panel tools ════
for (const b of $('fieldSeg').children) {
  b.addEventListener('click', () => {
    V.field = b.dataset.field;
    for (const o of $('fieldSeg').children) o.classList.toggle('on', o === b);
    streak.w = 0;
  });
}
$('gapLog').addEventListener('change', (e) => { V.log = e.target.checked; });
$('gapAuto').addEventListener('change', (e) => { V.auto = e.target.checked; });
$('viewSnap').addEventListener('click', () => { V.gapMode = 'snapshot'; $('viewSnap').classList.add('on'); $('viewStreak').classList.remove('on'); });
$('viewStreak').addEventListener('click', () => { V.gapMode = 'streak'; streak.w = 0; $('viewStreak').classList.add('on'); $('viewSnap').classList.remove('on'); });
$('iLog').addEventListener('click', () => { V.iLog = true; $('iLog').classList.add('on'); $('iLin').classList.remove('on'); });
$('iLin').addEventListener('click', () => { V.iLog = false; $('iLin').classList.add('on'); $('iLog').classList.remove('on'); });
$('showItot').addEventListener('change', (e) => { V.showItot = e.target.checked; });
$('qv1').addEventListener('click', () => { V.qvWindow = 1; $('qv1').classList.add('on'); $('qv8').classList.remove('on'); });
$('qv8').addEventListener('click', () => { V.qvWindow = 8; $('qv8').classList.add('on'); $('qv1').classList.remove('on'); });
$('profX').addEventListener('click', () => { V.profMode = 'x'; $('profX').classList.add('on'); $('profT').classList.remove('on'); });
$('profT').addEventListener('click', () => { V.profMode = 't'; $('profT').classList.add('on'); $('profX').classList.remove('on'); });
$('cvdSafe').addEventListener('change', (e) => { V.cvdSafe = e.target.checked; syncLegend(); });

// hover interactions
cvGap.el.addEventListener('mousemove', (e) => {
  if (!gapGeom) return;
  const r = cvGap.el.getBoundingClientRect();
  const x = e.clientX - r.left;
  gapHover = (x >= gapGeom.gx0 && x <= gapGeom.gx0 + gapGeom.wg)
    ? Math.floor(((x - gapGeom.gx0) / gapGeom.wg) * gapGeom.nx) : -1;
});
cvGap.el.addEventListener('mouseleave', () => { gapHover = -1; });
cvProf.el.addEventListener('mousemove', (e) => {
  if (!profGeom) return;
  const r = cvProf.el.getBoundingClientRect();
  const x = e.clientX - r.left;
  profHover = (x >= profGeom.padL && x <= profGeom.padL + profGeom.pw)
    ? Math.floor(((x - profGeom.padL) / profGeom.pw) * GAS.n) : -1;
});
cvProf.el.addEventListener('mouseleave', () => { profHover = -1; });

const seekFrom = (o, geom) => (e) => {
  if (!geom) return;
  const r = o.el.getBoundingClientRect();
  const f = (e.clientX - r.left) / r.width;
  setCursor(geom.t0 + f * (geom.t1 - geom.t0));
};
cvMacro.el.addEventListener('mousedown', (e) => { seekFrom(cvMacro, macroGeom)(e); macroDrag = true; });
window.addEventListener('mousemove', (e) => { if (macroDrag) seekFrom(cvMacro, macroGeom)(e); });
window.addEventListener('mouseup', () => { macroDrag = false; });
let macroDrag = false;
cvWave.el.addEventListener('click', (e) => {
  const r = cvWave.el.getBoundingClientRect();
  const padL = 46, padR = 10;
  const f = (e.clientX - r.left - padL) / Math.max(1, r.width - padL - padR);
  if (f < 0 || f > 1) return;
  const lens = Math.pow(10, V.lensLog);
  setCursor(V.cursorT - lens / 2 + f * lens);
});

// ═══════════════════════════════════════════════════════════ table / help ═══
function openModal(title, html) {
  $('modalTitle').textContent = title;
  $('modalBody').innerHTML = html;
  $('modal').hidden = false;
}
$('modalClose').addEventListener('click', () => { $('modal').hidden = true; });
$('modal').addEventListener('click', (e) => { if (e.target === $('modal')) $('modal').hidden = true; });

for (const b of document.querySelectorAll('[data-table]')) {
  b.addEventListener('click', () => showTable(b.dataset.table));
}
function showTable(which) {
  const st = solver.state;
  if (which === 'prof' || which === 'gap') {
    const keys = SPECIES.filter((s) => st.n[s.key]).map((s) => s.key);
    let rows = `<table><tr><th>x, mm</th>${keys.map((k) => `<th>${k}</th>`).join('')}<th>E, V/m</th><th>E/N, Td</th><th>ρ</th></tr>`;
    const stride = Math.max(1, Math.floor(GAS.n / 60));
    for (let i = GAS.i0; i <= GAS.i1; i += stride) {
      rows += `<tr><td>${((st.x[i] - GAS.x0) * 1e3).toFixed(4)}</td>${keys.map((k) => `<td>${sci(st.n[k][i])}</td>`).join('')}` +
        `<td>${sci(st.E[i])}</td><td>${sig3(st.EN[i])}</td><td>${sci(st.rho[i])}</td></tr>`;
    }
    openModal('Profiles @ t = ' + fmtTime(st.t), rows + '</table>');
  } else if (which === 'wave') {
    const h = solver.history;
    const start = (h.head - h.len + h.capacity) % h.capacity;
    const stride = Math.max(1, Math.floor(h.len / 200));
    let rows = '<table><tr><th>t, s</th><th>U_app, V</th><th>U_gap, V</th><th>I, A</th><th>Q, C</th></tr>';
    for (let s = 0; s < h.len; s += stride) {
      const i = (start + s) % h.capacity;
      rows += `<tr><td>${h.t[i].toExponential(6)}</td><td>${h.Uapp[i].toFixed(1)}</td><td>${h.Ugap[i].toFixed(1)}</td><td>${h.current[i].toExponential(3)}</td><td>${h.charge[i].toExponential(3)}</td></tr>`;
    }
    openModal('Waveform samples', rows + '</table>');
  } else {
    const ps = solver.periodStats || {};
    openModal('Manley report', `# dbd-o2 Q–V (Manley) report
C_diel     = ${sig3((ps.Cdiel || 0) * 1e12)} pF
C_cell     = ${sig3((ps.Ccell || 0) * 1e12)} pF
U_burn     = ${sig3(ps.UburnkV || 0)} kV
E/period   = ${sig3((ps.energyPerPeriodJ || 0) * 1e3)} mJ
P          = ${sig3(ps.powerW || 0)} W
params     = ${JSON.stringify(solver.params)}`);
  }
}
$('btnHelp').addEventListener('click', () => openModal('Keyboard shortcuts',
  `Space      play / pause
←  →       step frame / seek
Shift+← →  previous / next breakdown event
+  −       speed
1 … 5      colormap field (n_e, ρ, S_ion, |E|, E/N)
L          lin / log current
T          table view of profiles
E          export CSV
?          this help`));

// ═══════════════════════════════════════════════════════════ export ═════════
function download(name, mime, data) {
  const blob = new Blob([data], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
const paramHeader = () => {
  const hs = historySpan();
  return `# dbd-o2 v0.1 ${USE_MOCK ? '(MOCK SOLVER — not physical)' : ''}\n`
    + `# ${new Date().toISOString()}\n`
    + `# params ${JSON.stringify(solver.params)}\n`
    + `# ВНИМАНИЕ: выгружается только кольцевой буфер истории — последние `
    + `${(hs.span * 1e6).toFixed(2)} мкс (~${(hs.span * solver.params.freqKHz * 1e3).toFixed(2)} периодов) `
    + `из ${(solver.state.t * 1e6).toFixed(2)} мкс прогона\n`;
};

$('btnCSV').addEventListener('click', () => {
  const h = solver.history;
  const start = (h.head - h.len + h.capacity) % h.capacity;
  let s = paramHeader() + 't_s,U_app_V,U_gap_V,I_A,Q_C\n';
  for (let k = 0; k < h.len; k++) {
    const i = (start + k) % h.capacity;
    s += `${h.t[i].toExponential(9)},${h.Uapp[i].toExponential(6)},${h.Ugap[i].toExponential(6)},${h.current[i].toExponential(6)},${h.charge[i].toExponential(6)}\n`;
  }
  download('dbd-o2-timeseries.csv', 'text/csv', s);
});
$('btnJSON').addEventListener('click', () => {
  const ps = solver.periodStats || {};
  download('dbd-o2-run.json', 'application/json', JSON.stringify({
    version: '0.1', mock: USE_MOCK, timestamp: new Date().toISOString(),
    params: solver.params, view: { field: V.field, log: V.log, iLog: V.iLog, lensLog: V.lensLog },
    periodStats: ps, t: solver.state.t,
  }, null, 2));
});
$('btnPNG').addEventListener('click', () => {
  const list = [['Discharge gap', cvGap], ['Waveforms', cvWave], ['Q–V Lissajous', cvQV], ['Density profiles', cvProf]];
  const scale = 2, padx = 12, hdr = 26, foot = 30;
  const cw = Math.max(...list.map(([, c]) => c.w));
  const ch = list.reduce((a, [, c]) => a + c.h + hdr, 0);
  const out = document.createElement('canvas');
  out.width = (cw + padx * 2) * scale; out.height = (ch + foot + padx) * scale;
  const g = out.getContext('2d');
  g.scale(scale, scale);
  g.fillStyle = COL.surface1; g.fillRect(0, 0, cw + padx * 2, ch + foot + padx);
  let y = padx;
  for (const [name, c] of list) {
    text(g, name.toUpperCase(), padx, y + 12, { color: COL.muted, size: 11 });
    y += hdr;
    g.drawImage(c.el, padx, y, c.w, c.h);
    y += c.h;
  }
  const p = solver.params;
  text(g, `dbd-o2${USE_MOCK ? ' [MOCK]' : ''} · U₀ ${p.U0kV} kV · f ${p.freqKHz} kHz · gap ${p.gapMM} mm · ε_r ${p.epsR} · γ ${p.gamma} · t = ${fmtTime(solver.state.t)} · ${new Date().toISOString()}`,
    padx, y + 14, { color: COL.text2, size: 10 });
  out.toBlob((b) => {
    const url = URL.createObjectURL(b);
    const a = document.createElement('a'); a.href = url; a.download = 'dbd-o2-dashboard.png'; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
});

// ═══════════════════════════════════════════════════════════ keyboard ═══════
window.addEventListener('keydown', (e) => {
  if (e.target.matches('input,select,textarea')) return;
  const k = e.key;
  if (k === ' ') { e.preventDefault(); setRunning(!V.running); }
  else if (k === 'ArrowRight' && e.shiftKey) snapEvent(+1);
  else if (k === 'ArrowLeft' && e.shiftKey) snapEvent(-1);
  else if (k === 'ArrowRight') setCursor(V.cursorT + Math.pow(10, V.lensLog) / 10);
  else if (k === 'ArrowLeft') setCursor(V.cursorT - Math.pow(10, V.lensLog) / 10);
  else if (k === '+' || k === '=') { V.speedLog = Math.min(-2, V.speedLog + 0.25); $('speed').value = V.speedLog; updateSpeedLabel(); }
  else if (k === '-') { V.speedLog = Math.max(-9, V.speedLog - 0.25); $('speed').value = V.speedLog; updateSpeedLabel(); }
  else if (k >= '1' && k <= '5') { const b = $('fieldSeg').children[+k - 1]; if (b) b.click(); }
  else if (k === 'l' || k === 'L') { (V.iLog ? $('iLin') : $('iLog')).click(); }
  else if (k === 't' || k === 'T') showTable('prof');
  else if (k === 'e' || k === 'E') $('btnCSV').click();
  else if (k === '?') $('btnHelp').click();
  else if (k === 'Escape') $('modal').hidden = true;
});

// mobile controls drawer
$('btnDrawer').addEventListener('click', () => {
  $('panelCtl').scrollIntoView({ behavior: 'smooth', block: 'center' });
});

// ═══════════════════════════════════════════════════════════ boot ══════════
buildControls();
buildMetrics();
buildLegend();
applyPreset('filamentary', false);
$('speed').value = V.speedLog;
updateSpeedLabel(); updateLensLabel();
setRunning(true);
window.addEventListener('resize', () => { streak.w = 0; });
requestAnimationFrame(frame);
