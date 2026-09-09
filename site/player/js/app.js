// app.js — управляющий слой плеера: загрузка прогонов, состояние, транспорт,
// двухмасштабный таймлайн, режим сравнения «фото вкл / фото выкл».

import Playback from './loader.js';
import Scene from './scene.js';
import { COLORMAP_NAMES, colormapLabel, defaultColormap } from './colormaps.js';
import {
  WaveformPlot, LissajousPlot, RadialProfilePlot, AxialProfilePlot, PhotoPanel,
  PhotoSeriesBuilder, cellVolumes, updateAxialFromPlayback,
} from './plots.js';
import { MetricsPanel } from './metrics.js';
import { t as tr, fieldLabel, fieldUnit, dataText, applyStatic, initLangSwitch } from './i18n.js';


const DATA_ROOT = '../data/';

const RUNS = [
  { id: 'run-default', key: 'run.default' },
  { id: 'run-nophoto', key: 'run.nophoto' },
  { id: 'run-low', key: 'run.low' },
  { id: 'run-high', key: 'run.high' },
  { id: 'run-fast', key: 'run.fast' },
  { id: 'run-synth', key: 'run.synth' },
];

const BASE_FPS = 24;

const $ = (id) => document.getElementById(id);

const el = {
  runSel: $('runSel'), detailSel: $('detailSel'), cmpBtn: $('cmpBtn'), status: $('status'),
  slotA: $('slotA'), slotB: $('slotB'), canvasA: $('canvasA'), canvasB: $('canvasB'),
  overlay: $('overlay'), ovTitle: $('ovTitle'), ovText: $('ovText'),
  fieldList: $('fieldList'), segLog: $('segLog'), segScope: $('segScope'),
  cmapSel: $('cmapSel'), glowChk: $('glowChk'), glowRange: $('glowRange'),
  sigChk: $('sigChk'), probeChk: $('probeChk'),
  readout: $('readout').querySelector('tbody'), runInfo: $('runInfo'), limitNote: $('limitNote'),
  btnFirst: $('btnFirst'), btnPrev: $('btnPrev'), btnPlay: $('btnPlay'),
  btnNext: $('btnNext'), btnLast: $('btnLast'), speedSel: $('speedSel'),
  btnA: $('btnA'), btnB: $('btnB'), btnLoopClr: $('btnLoopClr'), loopInfo: $('loopInfo'),
  clock: $('clock'), timeline: $('timeline'),
  plots: $('plots'), plotsBtn: $('plotsBtn'), plotsSrc: $('plotsSrc'), photoStat: $('photoStat'),
  metrics: $('metrics'),
  cvWave: $('cvWave'), cvLiss: $('cvLiss'), cvAxial: $('cvAxial'),
  cvRadial: $('cvRadial'), cvPhoto: $('cvPhoto'),
  btnLogI: $('btnLogI'), btnWinZoom: $('btnWinZoom'), btnLissMode: $('btnLissMode'),
  btnAxMode: $('btnAxMode'), btnSigAsinh: $('btnSigAsinh'),
  expPng: $('expPng'), expWebm: $('expWebm'), expCsv: $('expCsv'),
};

const state = {
  runId: 'run-default',
  detail: 'compact',
  compare: false,
  field: 'n_e',
  logScale: true,
  scaleMode: 'frame',
  cmap: defaultColormap('n_e'),
  afterglow: false,
  decay: 0.86,
  showSigma: true,
  showProbe: true,
  playing: false,
  speed: 1,
  t: 0,
  zoomSpan: 4e-7,
  loop: null,
  loading: false,
  showPlots: true,
  logI: true,
  winFollowZoom: false,
  lissMode: 'period',
  axMode: 'axis',
  sigAsinh: false,
  recording: false,
};

// графики (plots.js / metrics.js) — создаются заново на каждый прогон
const P = {
  wf: null, lj: null, rp: null, ap: null, pp: null, mp: null,
  vol: null, Vgas: 0, wcum: null, ipk: 0, breakdowns: 0, T: 1e-4,
  token: 0, photo: null,
};

const sceneA = new Scene(el.canvasA);
const sceneB = new Scene(el.canvasB);
let pbA = null, pbB = null;
const cache = new Map();          // key -> Playback

// ────────────────────────────── утилиты ──────────────────────────────────
function fmtSci(v, d = 2) {
  if (!Number.isFinite(v)) return '—';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return (a >= 100 ? v.toFixed(1) : a >= 1 ? v.toFixed(3) : v.toFixed(4)).replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
  const e = Math.floor(Math.log10(a));
  return `${(v / Math.pow(10, e)).toFixed(d)}e${e}`;
}
function timeUnit(span) {
  if (span >= 1e-3) return { k: 1e3, u: tr('unit.ms') };
  if (span >= 1e-6) return { k: 1e6, u: tr('unit.us') };
  if (span >= 1e-9) return { k: 1e9, u: tr('unit.ns') };
  return { k: 1e12, u: tr('unit.ps') };
}
function fmtTime(t) {
  if (t >= 1e-3) return `${(t * 1e3).toFixed(4)} ${tr('unit.ms')}`;
  if (t >= 1e-7) return `${(t * 1e6).toFixed(4)} ${tr('unit.us')}`;
  return `${(t * 1e9).toFixed(2)} ${tr('unit.ns')}`;
}
function lowerBound(arr, x) {
  let lo = 0, hi = arr.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < x) lo = m + 1; else hi = m; }
  return lo;
}
function clamp(v, a, b) { return v < a ? a : v > b ? b : v; }

function showOverlay(title, text) {
  el.ovTitle.textContent = title;
  el.ovText.textContent = text || '';
  el.overlay.classList.remove('hidden');
}
function hideOverlay() { el.overlay.classList.add('hidden'); }
function setStatus(s, isErr) {
  el.status.textContent = s;
  el.status.classList.toggle('err', !!isErr);
}

// ────────────────────────────── загрузка ─────────────────────────────────
function dirFor(id, detail) {
  return `${DATA_ROOT}${id}${detail === 'full' ? '' : '-compact'}/`;
}

async function getRun(id, detail) {
  const key = `${id}|${detail}`;
  if (cache.has(key)) return cache.get(key);
  const dir = dirFor(id, detail);
  const pb = await Playback.load(dir, {
    cacheFrames: 96,
    onProgress: (got, total) => {
      const mb = (got / 1048576).toFixed(1);
      setStatus(total
        ? tr('status.loadingRunOf', { id, mb, total: (total / 1048576).toFixed(1) })
        : tr('status.loadingRun', { id, mb }));
    },
  });
  if (!pb.frameCount) throw new Error(tr('err.zeroFrames', { id }));
  if (cache.size > 4) { const k = cache.keys().next().value; cache.delete(k); }
  cache.set(key, pb);
  return pb;
}

async function loadAll() {
  if (state.loading) return;
  state.loading = true;
  state.playing = false; syncPlayBtn();
  hideOverlay();
  try {
    setStatus(tr('status.loading'));
    if (state.compare) {
      pbA = await getRun('run-default', state.detail);
      pbB = await getRun('run-nophoto', state.detail);
    } else {
      pbA = await getRun(state.runId, state.detail);
      pbB = null;
    }
    sceneA.setPlayback(pbA);
    sceneA.title = state.compare ? tr('scene.photoOn') : (RUNS.find((r) => r.id === state.runId)?.id || '');
    sceneA.subtitle = state.compare ? 'run-default' : '';
    if (pbB) { sceneB.setPlayback(pbB); sceneB.title = tr('scene.photoOff'); sceneB.subtitle = 'run-nophoto'; }
    else sceneB.setPlayback(null);
    el.slotB.classList.toggle('hidden', !state.compare);

    buildFieldChips();
    applyOptions();
    state.loop = null;
    state.t = pbA.frameTimes[0];
    state.zoomSpan = Math.min(4e-7, duration() / 4 || 4e-7);
    renderRunInfo();
    buildPlots(pbA);
    setStatus(tr('status.loaded', {
      runId: pbA.manifest.runId, frames: pbA.frameCount,
      nr: pbA.grid.nrOut, nz: pbA.grid.nzOut,
      tLast: (pbA.manifest.stats.tLast * 1e6).toFixed(3),
    }));
  } catch (err) {
    pbA = null; pbB = null;
    destroyPlots();
    sceneA.setPlayback(null); sceneB.setPlayback(null);
    sceneA.error = tr('err.noData');
    setStatus(String(err.message || err), true);
    showOverlay(tr('overlay.loadFailed'), tr('overlay.loadFailedText', {
      err: String(err.message || err),
      dir: dirFor(state.compare ? 'run-default' : state.runId, state.detail),
    }));
  } finally {
    state.loading = false;
  }
}

// ────────────────────────────── UI: поля ─────────────────────────────────
function buildFieldChips() {
  el.fieldList.innerHTML = '';
  const specs = pbA ? pbA.manifest.fields : [];
  const has = new Set(specs.map((f) => f.name));
  if (!has.has(state.field)) state.field = specs[0] ? specs[0].name : 'n_e';
  for (const f of specs) {
    const b = document.createElement('button');
    b.className = 'chip' + (f.name === state.field ? ' on' : '');
    b.textContent = shortLabel(f);
    b.title = `${fieldLabel(f)} [${fieldUnit(f)}]`;
    b.onclick = () => setField(f.name);
    el.fieldList.appendChild(b);
  }
}
function shortLabel(f) {
  switch (f.name) {
    case 'n_e': return 'n_e';
    case 'ionizRate': return tr('field.chip.ionizRate');
    case 'rho': return tr('field.chip.rho');
    case 'Emag': return '|E|';
    case 'EN': return 'E/N';
    case 'photoIonRate': return tr('field.chip.photoIonRate');
    case 'photoDetachRate': return tr('field.chip.photoDetachRate');
    case 'n_O3m': return 'O₃⁻';
    case 'n_O3': return 'O₃';
    default: return f.name;
  }
}

function syncCmapSel() {
  el.cmapSel.value = state.cmap;
}

/** Единственная точка смены поля: держит чипы, палитру и обе сцены в согласии. */
function setField(name) {
  state.field = name;
  state.cmap = defaultColormap(name);
  buildFieldChips(); syncCmapSel(); applyOptions();
}
/** Частичная смена опций отображения (используется UI и headless-проверками). */
function setOpts(o) {
  Object.assign(state, o);
  if (o.cmap) syncCmapSel();
  syncScaleButtons();
  applyOptions();
}
function syncScaleButtons() {
  for (const [root, v] of [[el.segLog, state.logScale ? 'log' : 'lin'], [el.segScope, state.scaleMode]]) {
    root.querySelectorAll('button').forEach((b) => b.classList.toggle('on', b.dataset.v === v));
  }
  el.glowChk.checked = state.afterglow;
  el.sigChk.checked = state.showSigma;
  el.probeChk.checked = state.showProbe;
}

function applyOptions() {
  const o = {
    scaleMode: state.scaleMode, logScale: state.logScale, cmap: state.cmap,
    afterglow: state.afterglow, decay: state.decay,
    showSigma: state.showSigma, showProbe: state.showProbe,
  };
  for (const s of [sceneA, sceneB]) { s.setField(state.field); s.setOptions(o); }
}

// ────────────────────────────── время/кадры ──────────────────────────────
function tMin() { return pbA ? pbA.frameTimes[0] : 0; }
function tMax() { return pbA ? pbA.frameTimes[pbA.frameCount - 1] : 1; }
function duration() { return Math.max(1e-12, tMax() - tMin()); }
function curIndex() { return pbA ? pbA.findFrame(state.t) : 0; }

function seek(t) {
  state.t = clamp(t, tMin(), tMax());
  if (state.afterglow) { /* след сохраняем: он и есть смысл режима */ }
}
function stepFrames(d) {
  if (!pbA) return;
  const i = clamp(curIndex() + d, 0, pbA.frameCount - 1);
  state.t = pbA.frameTimes[i];
}

function syncPlayBtn() { el.btnPlay.textContent = state.playing ? '❚❚' : '▶'; }

// ────────────────────────────── таймлайн ─────────────────────────────────
const TL = { padL: 46, padR: 40 };

function timelineBoxes(w, h) {
  const x0 = TL.padL, x1 = w - TL.padR;
  return {
    A: { x0, x1, y0: 14, h: 50, dens: { y: 65, h: 8 }, labY: 75 },
    B: { x0, x1, y0: 108, h: 50, labY: 162 },
  };
}

/** Амплитуды U и I по ВСЕМУ прогону (полный проход, кэш на прогон). */
function seriesNorm() {
  if (TL._norm && TL._normKey === pbA.manifest.runId + pbA.frameCount) return TL._norm;
  const s = pbA.series;
  let uMax = 0, iMax = 0;
  if (s && s.t) {
    for (let k = 0; k < s.t.length; k++) {
      const u = Math.abs(s.Uapp[k]); if (u > uMax) uMax = u;
      if (s.Itot) { const c = Math.abs(s.Itot[k]); if (c > iMax) iMax = c; }
    }
  }
  TL._norm = { uMax: uMax || 1, iMax: iMax || 1 };
  TL._normKey = pbA.manifest.runId + pbA.frameCount;
  return TL._norm;
}

function drawTimeline() {
  const cv = el.timeline;
  const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
  const w = cv.clientWidth, h = cv.clientHeight || 150;
  if (!w) return;
  const pw = Math.round(w * dpr), ph = Math.round(h * dpr);
  if (cv.width !== pw || cv.height !== ph) { cv.width = pw; cv.height = ph; }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  if (!pbA) return;

  const box = timelineBoxes(w, h);
  TL.box = box; TL.w = w;

  // ── верхняя дорожка: весь прогон
  const a0 = tMin(), a1 = tMax();
  drawTrack(ctx, box.A, a0, a1, { title: tr('tl.whole'), density: true, frames: false });
  // ── нижняя дорожка: лупа
  let span = clamp(state.zoomSpan, 2e-9, duration());
  let b0 = state.t - span / 2, b1 = state.t + span / 2;
  if (b0 < a0) { b0 = a0; b1 = a0 + span; }
  if (b1 > a1) { b1 = a1; b0 = Math.max(a0, a1 - span); }
  TL.zoom = [b0, b1];
  drawTrack(ctx, box.B, b0, b1, { title: tr('tl.zoom'), density: false, frames: true, local: true });

  // связка окон: показать окно лупы на верхней дорожке
  const xa = box.A.x0 + ((b0 - a0) / (a1 - a0)) * (box.A.x1 - box.A.x0);
  const xb = box.A.x0 + ((b1 - a0) / (a1 - a0)) * (box.A.x1 - box.A.x0);
  ctx.fillStyle = 'rgba(77,163,255,0.14)';
  ctx.fillRect(xa, box.A.y0, Math.max(1.5, xb - xa), box.A.h);
  ctx.strokeStyle = 'rgba(77,163,255,0.6)';
  ctx.lineWidth = 1;
  ctx.strokeRect(xa - 0.5, box.A.y0 + 0.5, Math.max(1.5, xb - xa) + 1, box.A.h - 1);
}

function drawTrack(ctx, b, t0, t1, opt) {
  const w = b.x1 - b.x0, h = b.h, y0 = b.y0;
  const span = t1 - t0 || 1e-12;
  const s = pbA.series;

  ctx.fillStyle = '#0c1017';
  ctx.fillRect(b.x0, y0, w, h);
  ctx.strokeStyle = '#1f2733'; ctx.lineWidth = 1;
  ctx.strokeRect(b.x0 + 0.5, y0 + 0.5, w - 1, h - 1);

  // сетка + подписи времени
  const un = timeUnit(span);
  const nTicks = Math.max(2, Math.min(10, Math.round(w / 110)));
  const rawStep = span / nTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((v) => v >= rawStep) || mag * 10;
  ctx.font = '9.5px "JetBrains Mono", ui-monospace, monospace';
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let t = Math.ceil(t0 / step) * step; t <= t1; t += step) {
    const x = b.x0 + ((t - t0) / span) * w;
    ctx.strokeStyle = 'rgba(255,255,255,0.05)';
    ctx.beginPath(); ctx.moveTo(x + 0.5, y0); ctx.lineTo(x + 0.5, y0 + h); ctx.stroke();
    ctx.fillStyle = '#6b7482';
    ctx.fillText((t * un.k).toFixed(step * un.k < 1 ? 3 : step * un.k < 10 ? 2 : 1), x, b.labY);
  }
  ctx.textAlign = 'left';
  ctx.fillStyle = '#8b93a1';
  ctx.fillText(un.u, b.x1 + 5, b.labY);

  // заголовок
  ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  ctx.font = '9px Inter, system-ui, sans-serif';
  ctx.fillStyle = '#59616e';
  ctx.fillText(opt.title, b.x0, y0 - 4);

  // цикл
  if (state.loop) {
    const [la, lb] = state.loop;
    const xa = clamp(b.x0 + ((la - t0) / span) * w, b.x0, b.x1);
    const xb = clamp(b.x0 + ((lb - t0) / span) * w, b.x0, b.x1);
    if (xb > xa) {
      ctx.fillStyle = 'rgba(250,178,25,0.13)';
      ctx.fillRect(xa, y0, xb - xa, h);
      ctx.strokeStyle = 'rgba(250,178,25,0.75)';
      ctx.beginPath();
      ctx.moveTo(xa + 0.5, y0); ctx.lineTo(xa + 0.5, y0 + h);
      ctx.moveTo(xb - 0.5, y0); ctx.lineTo(xb - 0.5, y0 + h);
      ctx.stroke();
    }
  }

  // кривые
  if (s && s.t && s.t.length > 1) {
    const i0 = Math.max(0, lowerBound(s.t, t0) - 1);
    const i1 = Math.min(s.t.length - 1, lowerBound(s.t, t1) + 1);
    const nPts = i1 - i0 + 1;
    const stride = Math.max(1, Math.floor(nPts / (w * 6)));

    // Верхняя дорожка — амплитуды всего прогона; лупа — амплитуды окна,
    // иначе форма импульса тока в окне 100 нс вырождается в прямую.
    let uMax, iMax;
    if (opt.local) {
      uMax = 0; iMax = 0;
      for (let k = i0; k <= i1; k++) {
        const u = Math.abs(s.Uapp[k]); if (u > uMax) uMax = u;
        if (s.Itot) { const c = Math.abs(s.Itot[k]); if (c > iMax) iMax = c; }
      }
      const g = seriesNorm();
      uMax = uMax || g.uMax; iMax = iMax || g.iMax;
    } else {
      ({ uMax, iMax } = seriesNorm());
    }
    const yMid = y0 + h / 2;

    ctx.strokeStyle = 'rgba(255,255,255,0.10)';
    ctx.beginPath(); ctx.moveTo(b.x0, yMid + 0.5); ctx.lineTo(b.x1, yMid + 0.5); ctx.stroke();

    const path = (get, scaleFn, color, wid) => {
      ctx.strokeStyle = color; ctx.lineWidth = wid;
      ctx.beginPath();
      let started = false;
      for (let k = i0; k <= i1; k += stride) {
        const x = b.x0 + ((s.t[k] - t0) / span) * w;
        const y = yMid - scaleFn(get(k)) * (h / 2 - 3);
        if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
      }
      ctx.stroke();
    };
    const lin = (v) => clamp(v / uMax, -1, 1);
    const den = Math.asinh(iMax / (iMax * 1e-4)) || 1;
    const sym = (v) => clamp(Math.asinh(v / (iMax * 1e-4)) / den, -1, 1);
    path((k) => s.Uapp[k], lin, 'rgba(77,163,255,0.85)', 1.4);
    if (s.Itot) path((k) => s.Itot[k], sym, 'rgba(240,165,60,0.9)', 1.2);

    ctx.font = '9px "JetBrains Mono", ui-monospace, monospace';
    ctx.textAlign = 'right'; ctx.textBaseline = 'top';
    ctx.fillStyle = 'rgba(77,163,255,0.85)';
    ctx.fillText(tr('tl.uAmp', { v: fmtSci(uMax) }), b.x1 - 4, y0 + 3);
    ctx.fillStyle = 'rgba(240,165,60,0.9)';
    ctx.fillText(tr('tl.iAmp', { v: fmtSci(iMax) }), b.x1 - 4, y0 + 14);
    if (opt.local) {
      ctx.fillStyle = '#59616e';
      ctx.fillText(tr('tl.localNorm'), b.x1 - 4, y0 + 25);
    }
  }

  // плотность кадров
  if (opt.density) {
    const bins = Math.max(1, Math.round(w));
    const counts = new Uint16Array(bins);
    const ft = pbA.frameTimes;
    for (let k = 0; k < ft.length; k++) {
      const u = (ft[k] - t0) / span;
      if (u < 0 || u > 1) continue;
      counts[Math.min(bins - 1, (u * bins) | 0)]++;
    }
    let mx = 1;
    for (let k = 0; k < bins; k++) if (counts[k] > mx) mx = counts[k];
    ctx.fillStyle = '#0c1017';
    ctx.fillRect(b.x0, b.dens.y, w, b.dens.h);
    for (let k = 0; k < bins; k++) {
      if (!counts[k]) continue;
      const q = counts[k] / mx;
      ctx.fillStyle = q > 0.5 ? `rgba(250,178,25,${0.35 + 0.65 * q})`
        : `rgba(77,163,255,${0.25 + 0.6 * q})`;
      ctx.fillRect(b.x0 + k, b.dens.y + 1, 1, b.dens.h - 2);
    }
    ctx.strokeStyle = '#1f2733';
    ctx.strokeRect(b.x0 + 0.5, b.dens.y + 0.5, w - 1, b.dens.h - 1);
    ctx.font = '8.5px Inter, system-ui, sans-serif';
    ctx.fillStyle = '#59616e'; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillText(tr('tl.frames'), b.x0 - 4, b.dens.y + b.dens.h / 2);
  }

  // отметки кадров в лупе
  if (opt.frames) {
    const ft = pbA.frameTimes;
    let k0 = lowerBound(ft, t0), k1 = lowerBound(ft, t1);
    if (k1 - k0 < 800) {
      ctx.strokeStyle = 'rgba(200,215,235,0.35)';
      ctx.beginPath();
      for (let k = k0; k < k1; k++) {
        const x = b.x0 + ((ft[k] - t0) / span) * w;
        ctx.moveTo(x + 0.5, y0 + h - 7); ctx.lineTo(x + 0.5, y0 + h - 1);
      }
      ctx.stroke();
    }
    const cur = pbA.findFrame(state.t);
    const xc = b.x0 + ((ft[cur] - t0) / span) * w;
    if (xc >= b.x0 - 2 && xc <= b.x1 + 2) {
      ctx.fillStyle = '#fab219';
      ctx.beginPath(); ctx.arc(xc, y0 + h - 4, 2.6, 0, 6.2832); ctx.fill();
    }
  }

  // курсор
  const xc = b.x0 + ((state.t - t0) / span) * w;
  if (xc >= b.x0 - 1 && xc <= b.x1 + 1) {
    ctx.strokeStyle = '#e8eaed'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(xc + 0.5, y0); ctx.lineTo(xc + 0.5, y0 + h); ctx.stroke();
    ctx.fillStyle = '#e8eaed';
    ctx.beginPath();
    ctx.moveTo(xc, y0); ctx.lineTo(xc - 4, y0 - 5); ctx.lineTo(xc + 4, y0 - 5); ctx.closePath();
    ctx.fill();
  }
}

// клики/драг по таймлайну
function timelineSeek(ev, capture) {
  if (!pbA || !TL.box) return;
  const rect = el.timeline.getBoundingClientRect();
  const x = ev.clientX - rect.left, y = ev.clientY - rect.top;
  const A = TL.box.A, B = TL.box.B;
  let box = null, t0 = 0, t1 = 0;
  if (y < (A.y0 + A.h + 22)) { box = A; t0 = tMin(); t1 = tMax(); }
  else { box = B; [t0, t1] = TL.zoom || [tMin(), tMax()]; }
  if (capture) TL.drag = box === A ? 'A' : 'B';
  const u = clamp((x - box.x0) / (box.x1 - box.x0), 0, 1);
  seek(t0 + u * (t1 - t0));
}

el.timeline.addEventListener('pointerdown', (e) => {
  el.timeline.setPointerCapture(e.pointerId);
  state.playing = false; syncPlayBtn();
  timelineSeek(e, true);
});
el.timeline.addEventListener('pointermove', (e) => { if (TL.drag) timelineSeek(e, false); });
el.timeline.addEventListener('pointerup', () => { TL.drag = null; });
el.timeline.addEventListener('pointercancel', () => { TL.drag = null; });
el.timeline.addEventListener('wheel', (e) => {
  if (!TL.box) return;
  const rect = el.timeline.getBoundingClientRect();
  const y = e.clientY - rect.top;
  if (y < TL.box.B.y0 - 10) return;
  e.preventDefault();
  state.zoomSpan = clamp(state.zoomSpan * Math.exp(e.deltaY * 0.0015), 2e-9, duration());
}, { passive: false });

// ────────────────────────────── боковые данные ───────────────────────────
const READOUT_ROWS = [
  ['ro.t', (s, i) => fmtTime(s.t[i])],
  ['ro.Uapp', (s, i) => fmtSci(s.Uapp[i])],
  ['ro.Ugap', (s, i) => fmtSci(s.Ugap[i])],
  ['ro.Icond', (s, i) => fmtSci(s.Icond[i])],
  ['ro.Idisp', (s, i) => fmtSci(s.Idisp[i])],
  ['ro.Itot', (s, i) => fmtSci(s.Itot[i])],
  ['ro.maxEN', (s, i) => fmtSci(s.maxEN[i])],
  ['ro.sigmaMax', (s, i) => fmtSci(s.sigmaMax[i])],
  ['ro.o3ppm', (s, i) => fmtSci(s.o3ppm[i])],
  ['ro.photoL', (s, i) => fmtSci(s.photoEmitTotalL ? s.photoEmitTotalL[i] : NaN)],
  ['ro.photoR', (s, i) => fmtSci(s.photoEmitTotalR ? s.photoEmitTotalR[i] : NaN)],
];

function renderReadout() {
  if (!pbA || !pbA.series || !pbA.series.t) { el.readout.innerHTML = ''; return; }
  const s = pbA.series;
  let i = lowerBound(s.t, state.t);
  i = clamp(i, 0, s.t.length - 1);
  if (i > 0 && Math.abs(s.t[i - 1] - state.t) < Math.abs(s.t[i] - state.t)) i--;
  let html = '';
  for (const [k, f] of READOUT_ROWS) {
    let v;
    try { v = f(s, i); } catch (e) { v = '—'; }
    html += `<tr><td>${tr(k)}</td><td>${v}</td></tr>`;
  }
  el.readout.innerHTML = html;
}

function renderRunInfo() {
  if (!pbA) { el.runInfo.textContent = '—'; el.limitNote.textContent = ''; return; }
  const m = pbA.manifest, p = m.params || {};
  const crash = m.validation && m.validation.run && m.validation.run.crashed;
  let html = `<b>${m.runId}</b> · ${m.level}<br>`;
  html += `${tr('info.params', { U0: p.U0kV, f: p.freqKHz })}<br>`;
  html += `${tr('info.photo', { on: !!p.photoModule })}<br>`;
  html += `${tr('info.grid', { nr: p.nr, nz: p.nz, nrOut: m.grid.nrOut, nzOut: m.grid.nzOut })}<br>`;
  html += `${tr('info.counts', { frames: m.frames.length, steps: m.stats.steps })}<br>`;
  html += tr('info.tspan', { tLast: (m.stats.tLast * 1e6).toFixed(4) });
  if (crash) {
    html += `<div class="crash">${tr('info.crash', { t: (crash.t * 1e6).toFixed(4) })}`
      + `<br>${dataText(String(crash.message)).slice(0, 420)}</div>`;
  }
  el.runInfo.innerHTML = html;
  el.limitNote.textContent = (m.limits && m.limits.note)
    ? dataText(m.limits.note) : tr('limit.default');
}

// ────────────────────────────── графики ──────────────────────────────────
function destroyPlots() {
  P.token++;
  for (const k of ['wf', 'lj', 'rp', 'ap', 'pp']) { if (P[k]) P[k].destroy(); P[k] = null; }
  if (P.mp) { P.mp.destroy && P.mp.destroy(); P.mp = null; }
  el.metrics.innerHTML = '';
  P.photo = null;
}

function buildPlots(pb) {
  destroyPlots();
  const my = P.token;
  const par = pb.manifest.params || {};
  P.T = 1 / ((par.freqKHz || 10) * 1e3);
  P.vol = cellVolumes(pb.grid);
  P.Vgas = P.vol.reduce((a, b) => a + b, 0);

  P.wf = new WaveformPlot(el.cvWave, pb.series, {
    logCurrent: state.logI,
    onSeek: (t) => { state.playing = false; syncPlayBtn(); seek(t); },
  });
  P.lj = new LissajousPlot(el.cvLiss, pb.series, {
    freqHz: (par.freqKHz || 10) * 1e3, mode: state.lissMode,
    geometry: { epsR: par.epsR, gapMM: par.gapMM, dielMM: par.dielMM, radiusMM: par.radiusMM },
  });
  P.rp = new RadialProfilePlot(el.cvRadial, { scale: state.sigAsinh ? 'asinh' : 'lin' });
  P.ap = new AxialProfilePlot(el.cvAxial);
  P.pp = new PhotoPanel(el.cvPhoto);
  P.mp = new MetricsPanel(el.metrics);

  const s = pb.series;
  P.wf.setWindow(pb.frameTimes[0], s ? s.t[s.t.length - 1] : tMax());
  P.pp.setWindow(pb.frameTimes[0], pb.frameTimes[pb.frameCount - 1]);

  // накопленная энергия ∫U·I dt, пик тока, число пробоев — один проход по ряду
  P.ipk = 0; P.wcum = null; P.breakdowns = 0;
  if (s && s.t) {
    for (let k = 0; k < s.t.length; k++) P.ipk = Math.max(P.ipk, Math.abs(s.Icond[k]));
    P.wcum = new Float64Array(s.t.length);
    for (let k = 1; k < s.t.length; k++) {
      P.wcum[k] = P.wcum[k - 1]
        + 0.5 * (s.Uapp[k] * s.Itot[k] + s.Uapp[k - 1] * s.Itot[k - 1]) * (s.t[k] - s.t[k - 1]);
    }
    let above = false;
    for (let k = 0; k < s.t.length; k++) {
      const a = P.ipk > 0 && Math.abs(s.Icond[k]) > 0.05 * P.ipk;
      if (a && !above) P.breakdowns++;
      above = a;
    }
  }

  el.plotsSrc.textContent = state.compare ? tr('metrics.srcLeft') : pb.manifest.runId;
  el.photoStat.textContent = tr('metrics.photoIntegrating');
  P.pp.setData(null);
  PhotoSeriesBuilder.build(pb, {}).then((d) => {
    if (my !== P.token || !P.pp) return;          // прогон уже сменили
    P.photo = d; P.pp.setData(d);
    el.photoStat.textContent = d.missing && d.missing.length
      ? tr('metrics.photoMissing', { list: d.missing.join(', ') }) : tr('metrics.photoRate');
  }).catch((e) => { if (my === P.token) el.photoStat.textContent = tr('metrics.error', { msg: e.message }); });
}

let lastMetricsTs = 0;
function updatePlots(t, i) {
  if (!pbA || !P.wf || !state.showPlots) return;
  if (state.winFollowZoom && TL.zoom) P.wf.setWindow(TL.zoom[0], TL.zoom[1]);
  P.wf.setCursor(t); P.lj.setCursor(t); P.pp.setCursor(t);
  const su = pbA.getSurface(i);
  P.rp.setData(pbA.grid.rCenters, su.sigmaL, su.sigmaR);
  updateAxialFromPlayback(P.ap, pbA, i, { mode: state.axMode, vol: P.vol });

  const now = performance.now();
  if (now - lastMetricsTs < 250) return;          // метрики — 4 Гц (UI_SPEC §4.1)
  lastMetricsTs = now;
  const s = pbA.series;
  if (!s || !s.t) return;
  const k = clamp(lowerBound(s.t, t), 0, s.t.length - 1);
  const fit = P.lj.fit;
  P.mp.update({
    // На оборванном прогоне петля Лиссажу не замкнута: C_cell из неё ещё
    // извлекается (ёмкостная ветвь чистая), а C_diel, мощность и «энергия за
    // период» — нет. Показывать их числом здесь значило бы обойти собственный
    // бейдж графика, поэтому при fit.warn они гасятся (UI_SPEC §2.3).
    power: fit && fit.ok && !fit.warn ? fit.P : null,
    specEnergy: { v: P.wcum ? P.wcum[k] / P.Vgas / 1000 : null, note: tr('metrics.note.cumulative') },
    ozone: s.o3ppm[k],
    maxEN: s.maxEN[k],
    Ipeak: P.ipk,
    breakdowns: P.breakdowns / Math.max(1, s.t[s.t.length - 1] / P.T),
    C_diel: fit && fit.ok && !fit.warn ? fit.C_diel : null,
    C_cell: fit && fit.ok ? { v: fit.C_cell, note: tr('metrics.note.medianDQ') } : null,
    energy: fit && fit.ok
      ? { v: fit.W, note: fit.warn ? tr('metrics.note.loopOpen') : tr('metrics.note.perPeriod') } : null,
    _status: {
      frame: `${i + 1}/${pbA.frameCount}`,
      dt: s.t[Math.min(k + 1, s.t.length - 1)] - s.t[k],
      stepsPerFrame: pbA.manifest.stats.steps / pbA.frameCount,
      limiter: fit && fit.warn ? tr('metrics.limiterLiss') : 'ok',
    },
  });
}

// ────────────────────────────── экспорт ──────────────────────────────────
function download(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 4000);
}
function stamp() {
  return `${(state.t * 1e9).toFixed(1)}ns`;
}
function runTag() {
  return state.compare ? 'cmp-default-nophoto' : state.runId;
}

/** Композит сцены (одной или двух) + подпись — общий для PNG и WebM. */
function composite(cv) {
  const list = state.compare && pbB ? [el.canvasA, el.canvasB] : [el.canvasA];
  const gap = 4, cap = 26;
  const w = list.reduce((a, c) => a + c.width, 0) + gap * (list.length - 1);
  const h = Math.max(...list.map((c) => c.height)) + cap;
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  const g = cv.getContext('2d');
  g.fillStyle = '#080b10'; g.fillRect(0, 0, w, h);
  let x = 0;
  for (const c of list) { g.drawImage(c, x, cap); x += c.width + gap; }
  g.fillStyle = '#e8eaed';
  g.font = `${Math.round(cap * 0.52)}px ui-monospace, "JetBrains Mono", monospace`;
  g.textBaseline = 'middle'; g.textAlign = 'left';
  const f = pbA && pbA.fieldByName.get(state.field);
  g.fillText(tr('export.caption', {
    tag: runTag(), field: f ? fieldLabel(f) : state.field, t: fmtTime(state.t),
    scaleKind: state.logScale ? 'log' : 'lin',
    scope: state.scaleMode === 'global' ? tr('scale.global') : tr('scale.frame'),
  }), 8, cap / 2);
  return cv;
}

function exportPNG() {
  if (!pbA) return;
  const cv = composite(document.createElement('canvas'));
  cv.toBlob((b) => download(b, `dbd2d-${runTag()}-${state.field}-${stamp()}.png`), 'image/png');
  setStatus(tr('status.pngSaved', { field: state.field, t: fmtTime(state.t) }));
}

function exportCSV() {
  if (!pbA || !pbA.series) return;
  const s = pbA.series;
  const names = pbA.manifest.series.names;
  const n = s.t.length;
  const parts = [names.join(',') + '\n'];
  const chunk = [];
  for (let i = 0; i < n; i++) {
    chunk.length = 0;
    for (const nm of names) chunk.push(s[nm][i].toExponential(9));
    parts.push(chunk.join(',') + '\n');
  }
  download(new Blob(parts, { type: 'text/csv' }), `dbd2d-${pbA.manifest.runId}-series.csv`);
  setStatus(tr('status.csv', { rows: n, cols: names.length }));
}

const REC = { mr: null, cv: null, chunks: null, stopAt: 0 };
function toggleWebM() {
  if (state.recording) { stopWebM(); return; }
  if (!pbA) return;
  if (typeof MediaRecorder === 'undefined') { setStatus(tr('status.webmUnsupported'), true); return; }
  REC.cv = document.createElement('canvas');
  composite(REC.cv);
  const stream = REC.cv.captureStream(30);
  const mime = ['video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm']
    .find((m) => MediaRecorder.isTypeSupported(m)) || '';
  REC.chunks = [];
  try { REC.mr = new MediaRecorder(stream, mime ? { mimeType: mime, videoBitsPerSecond: 8e6 } : undefined); }
  catch (e) { setStatus('WebM: ' + e.message, true); return; }
  REC.mr.ondataavailable = (e) => { if (e.data && e.data.size) REC.chunks.push(e.data); };
  REC.mr.onstop = () => {
    download(new Blob(REC.chunks, { type: 'video/webm' }), `dbd2d-${runTag()}-${state.field}.webm`);
    setStatus(tr('status.webmDone', { mb: (REC.chunks.reduce((a, b) => a + b.size, 0) / 1048576).toFixed(1) }));
  };
  REC.mr.start(200);
  // пишем от текущего кадра до конца прогона (или до конца участка A–B)
  REC.stopAt = state.loop ? state.loop[1] : tMax();
  state.recording = true;
  state.playing = true; syncPlayBtn();
  el.expWebm.classList.add('rec');
  el.expWebm.textContent = tr('btn.webm.stop');
  setStatus(tr('status.webmRecording'));
}
function stopWebM() {
  state.recording = false;
  el.expWebm.classList.remove('rec');
  el.expWebm.textContent = tr('btn.webm');
  if (REC.mr && REC.mr.state !== 'inactive') REC.mr.stop();
  REC.mr = null;
}

// ────────────────────────────── главный цикл ─────────────────────────────
let lastTs = 0, acc = 0;
function frameLoop(ts) {
  const dt = lastTs ? Math.min(0.1, (ts - lastTs) / 1000) : 0;
  lastTs = ts;

  if (state.playing && pbA) {
    acc += dt * BASE_FPS * state.speed;
    while (acc >= 1) {
      acc -= 1;
      let i = curIndex() + 1;
      const lo = state.loop ? pbA.findFrame(state.loop[0]) : 0;
      const hi = state.loop ? pbA.findFrame(state.loop[1]) : pbA.frameCount - 1;
      if (i > hi) { i = lo; sceneA.resetAfterglow(); sceneB.resetAfterglow(); }
      state.t = pbA.frameTimes[clamp(i, 0, pbA.frameCount - 1)];
    }
  }

  if (pbA) {
    const iA = curIndex();
    sceneA.setFrameIndex(iA);
    if (state.compare && pbB) {
      // Сравнение честно только при ОБЩЕЙ шкале: иначе покадровая автошкала
      // выравнивает яркость и разница фото вкл/выкл визуально исчезает.
      const iB = pbB.findFrame(state.t);
      sceneB.setFrameIndex(iB);
      const amp = (pb, i) => {
        const c = pb.getCodec(i, state.field) || {};
        return Math.max(Math.abs(c.max || 0), Math.abs(c.min || 0));
      };
      const m = state.scaleMode === 'global'
        ? Math.max(sceneA.globalMax(state.field), sceneB.globalMax(state.field))
        : Math.max(amp(pbA, iA), amp(pbB, iB));
      sceneA.setForceMax(m); sceneB.setForceMax(m);
      sceneA.render(); sceneB.render();
    } else {
      sceneA.setForceMax(0); sceneB.setForceMax(0);
      sceneA.render();
    }
    drawTimeline();
    renderReadout();
    updatePlots(state.t, iA);
    if (state.recording) {
      composite(REC.cv);
      if (state.t >= REC.stopAt - 1e-15 || !state.playing) stopWebM();
    }
    const i = curIndex();
    el.clock.textContent = tr('tl.clock', { t: fmtTime(state.t), i: i + 1, n: pbA.frameCount });
    el.loopInfo.textContent = state.loop
      ? tr('tl.loop', { a: fmtTime(state.loop[0]), b: fmtTime(state.loop[1]) }) : '';
  } else {
    sceneA.render();
  }
  requestAnimationFrame(frameLoop);
}

// ────────────────────────────── события ──────────────────────────────────
function initControls() {
  for (const r of RUNS) {
    const o = document.createElement('option');
    o.value = r.id; o.textContent = `${r.id} — ${tr(r.key)}`;
    el.runSel.appendChild(o);
  }
  el.runSel.value = state.runId;
  for (const c of COLORMAP_NAMES) {
    const o = document.createElement('option');
    o.value = c; o.textContent = colormapLabel(c);
    el.cmapSel.appendChild(o);
  }
  syncCmapSel();

  el.runSel.onchange = () => { state.runId = el.runSel.value; state.compare = false; el.cmpBtn.classList.remove('on'); loadAll(); };
  el.detailSel.onchange = () => { state.detail = el.detailSel.value; loadAll(); };
  el.cmpBtn.onclick = () => {
    state.compare = !state.compare;
    el.cmpBtn.classList.toggle('on', state.compare);
    el.runSel.disabled = state.compare;
    loadAll();
  };

  const seg = (root, key, cb) => {
    root.querySelectorAll('button').forEach((b) => {
      b.onclick = () => {
        root.querySelectorAll('button').forEach((x) => x.classList.remove('on'));
        b.classList.add('on');
        cb(b.dataset.v);
      };
    });
  };
  seg(el.segLog, 'log', (v) => { state.logScale = v === 'log'; applyOptions(); });
  seg(el.segScope, 'scope', (v) => { state.scaleMode = v; applyOptions(); });

  el.cmapSel.onchange = () => { state.cmap = el.cmapSel.value; applyOptions(); };
  el.glowChk.onchange = () => { state.afterglow = el.glowChk.checked; applyOptions(); };
  el.glowRange.oninput = () => { state.decay = Number(el.glowRange.value) / 100; applyOptions(); };
  el.sigChk.onchange = () => { state.showSigma = el.sigChk.checked; applyOptions(); };
  el.probeChk.onchange = () => { state.showProbe = el.probeChk.checked; applyOptions(); };

  el.btnPlay.onclick = () => { state.playing = !state.playing; syncPlayBtn(); };
  el.btnPrev.onclick = () => { state.playing = false; syncPlayBtn(); stepFrames(-1); };
  el.btnNext.onclick = () => { state.playing = false; syncPlayBtn(); stepFrames(1); };
  el.btnFirst.onclick = () => { state.t = tMin(); sceneA.resetAfterglow(); sceneB.resetAfterglow(); };
  el.btnLast.onclick = () => { state.t = tMax(); };
  el.speedSel.onchange = () => { state.speed = Number(el.speedSel.value); };

  el.btnA.onclick = () => {
    const b = state.loop ? state.loop[1] : tMax();
    state.loop = [Math.min(state.t, b - 1e-12), b];
  };
  el.btnB.onclick = () => {
    const a = state.loop ? state.loop[0] : tMin();
    state.loop = [a, Math.max(state.t, a + 1e-12)];
  };
  el.btnLoopClr.onclick = () => { state.loop = null; };

  // ── графики и экспорт
  el.plotsBtn.onclick = () => {
    state.showPlots = !state.showPlots;
    el.plotsBtn.classList.toggle('on', state.showPlots);
    el.plots.classList.toggle('hidden', !state.showPlots);
  };
  el.btnLogI.onclick = () => {
    state.logI = !state.logI;
    el.btnLogI.classList.toggle('on', state.logI);
    if (P.wf) P.wf.setLogCurrent(state.logI);
  };
  el.btnWinZoom.onclick = () => {
    state.winFollowZoom = !state.winFollowZoom;
    el.btnWinZoom.classList.toggle('on', state.winFollowZoom);
    if (!state.winFollowZoom && P.wf && pbA) {
      const s = pbA.series;
      P.wf.setWindow(pbA.frameTimes[0], s ? s.t[s.t.length - 1] : tMax());
    }
  };
  const LISS = ['period', 'last8', 'all'];
  const LISS_LAB = { period: 'liss.period', last8: 'liss.last8', all: 'liss.all' };
  el.btnLissMode.onclick = () => {
    state.lissMode = LISS[(LISS.indexOf(state.lissMode) + 1) % LISS.length];
    el.btnLissMode.textContent = tr(LISS_LAB[state.lissMode]);
    if (P.lj) P.lj.setMode(state.lissMode);
  };
  el.btnAxMode.onclick = () => {
    state.axMode = state.axMode === 'axis' ? 'mean' : 'axis';
    el.btnAxMode.textContent = tr(state.axMode === 'axis' ? 'ax.axis' : 'ax.mean');
  };
  el.btnSigAsinh.onclick = () => {
    state.sigAsinh = !state.sigAsinh;
    el.btnSigAsinh.classList.toggle('on', state.sigAsinh);
    if (P.rp) P.rp.setScaleMode(state.sigAsinh ? 'asinh' : 'lin');
  };
  el.expPng.onclick = exportPNG;
  el.expCsv.onclick = exportCSV;
  el.expWebm.onclick = toggleWebM;

  // зонд
  for (const [cv, sc] of [[el.canvasA, sceneA], [el.canvasB, sceneB]]) {
    cv.addEventListener('pointermove', (e) => {
      const r = cv.getBoundingClientRect();
      sc.setPointer({ x: e.clientX - r.left, y: e.clientY - r.top });
    });
    cv.addEventListener('pointerleave', () => sc.setPointer(null));
  }

  window.addEventListener('keydown', (e) => {
    if (e.target && /^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) return;
    const big = e.shiftKey ? 10 : 1;
    switch (e.key) {
      case ' ': e.preventDefault(); state.playing = !state.playing; syncPlayBtn(); break;
      case 'ArrowLeft': e.preventDefault(); state.playing = false; syncPlayBtn(); stepFrames(-big); break;
      case 'ArrowRight': e.preventDefault(); state.playing = false; syncPlayBtn(); stepFrames(big); break;
      case 'Home': state.t = tMin(); break;
      case 'End': state.t = tMax(); break;
      case '[': el.btnA.onclick(); break;
      case ']': el.btnB.onclick(); break;
      case '\\': state.loop = null; break;
      default: break;
    }
  });
}

applyStatic();
initLangSwitch();
initControls();
loadAll();
requestAnimationFrame(frameLoop);

// диагностика для headless-проверки
window.__player = {
  state, get pbA() { return pbA; }, get pbB() { return pbB; }, sceneA, sceneB, loadAll,
  plots: P, seek, exportPNG, exportCSV, composite, setField, setOpts,
};
