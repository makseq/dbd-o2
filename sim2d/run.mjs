#!/usr/bin/env node
// run.mjs — CLI-раннер 2D осесимметричного ДБР в чистом O2.
//
// Опорные документы: docs/ERRATA.md (приоритетен), docs/NUMERICS_2D.md,
// docs/PHOTO_PROCESSES.md, docs/REFERENCE_TARGETS.md (ДЕЙСТВУЮЩИЕ эталоны).
//
// Использование:
//   node sim2d/run.mjs --preset validation --periods 2 --out data/val
//   node sim2d/run.mjs --preset filament --periods 1 --nr 48 --nz 168 --photo off
//
// Прогресс идёт в stderr (stdout остаётся чистым: туда пишется ТОЛЬКО итоговый JSON,
// чтобы раннер можно было включать в пайплайн `node run.mjs ... | jq`).

import fs from 'node:fs';
import path from 'node:path';
import { DBD2D, DEFAULTS } from './solver2d.mjs';
import { createRecorders } from './recorder.mjs';
import { EPS0 } from './physics2d.mjs';

// ──────────────────────────────────────────────────────────────── разбор CLI

const HELP = `
node sim2d/run.mjs [опции]

  --preset <name>     validation | filament | quick   (по умолчанию validation)
  --periods <x>       сколько периодов считать (дробное можно), по умолчанию 2
  --out <dir>         каталог прогона (по умолчанию data/<preset>)
  --nr <n> --nz <n>   сетка
  --photo on|off      все три фотопроцесса разом
  --progress <n>      печатать прогресс каждые n шагов (0 = выключить), по умолчанию 2000
  --record on|off     писать frames.bin/series.bin (по умолчанию on)
  --level full,compact  какие уровни писать
  --maxSeconds <s>    жёсткий бюджет wall-clock; по достижении прогон закрывается штатно
  --maxSteps <n>      то же по числу шагов
  --resumeless        не писать чекпойнт при аварии

  Любой параметр DBD2D можно переопределить напрямую:
  --U0kV --freqKHz --gammaIon --seedSpotAmp --seedBackground --d6MassConvention
  --chemSafety --chemActiveFrac --dtMax --dtMin --strictGauss --betaG --nzDiel ...
`;

function parseArgs(argv) {
  const out = { _: [] };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (!a.startsWith('--')) { out._.push(a); continue; }
    const key = a.slice(2);
    const eq = key.indexOf('=');
    if (eq >= 0) { out[key.slice(0, eq)] = key.slice(eq + 1); continue; }
    const nxt = argv[i + 1];
    if (nxt === undefined || nxt.startsWith('--')) { out[key] = 'true'; continue; }
    out[key] = nxt; i++;
  }
  return out;
}

const BOOL = { on: true, off: false, true: true, false: false, yes: true, no: false, 1: true, 0: false };
function asBool(v, d) { if (v === undefined) return d; const b = BOOL[String(v).toLowerCase()]; return b === undefined ? d : b; }
function asNum(v, d) { if (v === undefined) return d; const x = Number(v); return Number.isFinite(x) ? x : d; }

// ──────────────────────────────────────────────────────────────── пресеты

// ВАЖНО (ERRATA B6): при радиально ОДНОРОДНОЙ затравке решение остаётся радиально
// однородным навсегда — это не баг, а свойство осесимметричной постановки. Поэтому
// валидационный режим (seedSpotAmp = 0) имеет право работать на узкой r-сетке:
// он проверяет схемные/интегральные критерии V4..V13, которым радиальная структура
// не нужна. Филамент требует ОТДЕЛЬНОГО пресета с затравочным пятном.
const PRESETS = {
  validation: {
    seedSpotAmp: 0,          // радиально однородно
    nr: 4, nz: 140, nzDiel: 14,
    betaR: 0, betaG: 3.0,
    photoModule: true,
    dtMax: 2e-9,
  },
  filament: {
    seedSpotAmp: 1e18, seedSpotSigmaUM: 40,
    nr: 48, nz: 168, nzDiel: 14,
    betaR: 1.6, betaG: 3.0,
    photoModule: true,
  },
  quick: {
    seedSpotAmp: 0,
    nr: 2, nz: 84, nzDiel: 8,
    betaR: 0, betaG: 2.5,
    photoModule: true,
    dtMax: 5e-9,
  },
};

// параметры DBD2D, которые можно переопределить с командной строки «как есть»
const NUM_KEYS = new Set(Object.keys(DEFAULTS).filter((k) => typeof DEFAULTS[k] === 'number'));
const BOOL_KEYS = new Set(Object.keys(DEFAULTS).filter((k) => typeof DEFAULTS[k] === 'boolean'));

export function buildParams(args) {
  const preset = args.preset || 'validation';
  if (!PRESETS[preset]) throw new Error(`неизвестный пресет '${preset}', доступны: ${Object.keys(PRESETS)}`);
  const p = { ...PRESETS[preset] };
  if (args.nr !== undefined) p.nr = asNum(args.nr, p.nr);
  if (args.nz !== undefined) p.nz = asNum(args.nz, p.nz);
  if (args.photo !== undefined) {
    const on = asBool(args.photo, true);
    p.photoIonization = on; p.photoEmission = on; p.photoDetachment = on;
    p.photoModule = on;
  }
  for (const k of NUM_KEYS) if (args[k] !== undefined) p[k] = asNum(args[k], DEFAULTS[k]);
  for (const k of BOOL_KEYS) if (args[k] !== undefined) p[k] = asBool(args[k], DEFAULTS[k]);
  if (args.d6MassConvention !== undefined) p.d6MassConvention = args.d6MassConvention;
  if (args.wallBC !== undefined) p.wallBC = args.wallBC;
  if (args.nzDiel !== undefined) p.nzDiel = asNum(args.nzDiel, p.nzDiel);
  return { preset, params: p };
}

// ──────────────────────────────────────────────────── накопитель диагностики

/**
 * Онлайновый анализ для VALIDATION.md. Считает всё, что нужно для V4..V12,
 * прямо по ходу прогона (по КАЖДОМУ шагу, а не по прореженным кадрам).
 */
class Analyzer {
  constructor(sim, periods) {
    this.sim = sim;
    this.period = sim.period;
    this.periods = periods;
    this.reset();
    // геометрия ячейки — для сверки с эталонами на 1 см^2
    const S = sim.poisson;
    this.area = Math.PI * S.R * S.R;                          // м^2
    this.Cd = (EPS0 * sim.p.epsR * this.area) / (sim.p.dielMM * 1e-3) / 2; // два барьера последовательно
    this.Cg = (EPS0 * this.area) / (sim.p.gapMM * 1e-3);
    this.Ccell = 1 / (1 / this.Cd + 1 / this.Cg);
  }

  reset() {
    this.lyss = [];        // [Uapp, Q] прорежённая петля последнего периода
    this.pulses = [];      // импульсы тока
    this._inPulse = false;
    this.energy = 0;       // ∫ U dQ = ∫ U I dt
    this.energyLast = 0;   // за последний полный период
    this.tLastPeriodStart = null;
    this.sigmaMaxRun = 0;
    this.maxENrun = 0;
    this.qCheckStart = null;
    this.o3max = 0;
    this.neMaxRun = 0;
    this.samples = [];     // прорежённые ряды для отчёта
    this._prev = null;
    this.QatHalf = [];     // заряд в узлах смены полярности
    this.chargeErr = 0;
    this.negOverE = [];    // моменты «n(-) > n_e» после импульса
    this.halves = [];      // помуприодная бухгалтерия заряда (V6)
    this._hIdx = -1;
  }

  step(st, dt) {
    const s = this.sim;
    const I = st.Itot, U = st.Uapp;
    this.energy += U * I * dt;
    if (st.sigmaMax > this.sigmaMaxRun) this.sigmaMaxRun = st.sigmaMax;
    if (st.maxEN > this.maxENrun) this.maxENrun = st.maxEN;
    if (st.o3ppm > this.o3max) this.o3max = st.o3ppm;
    if (s.neMax > this.neMaxRun) this.neMaxRun = s.neMax;

    // --- снимок филамента в момент ГЛОБАЛЬНОГО пика тока (для развёртки по фото) ---
    if (Math.abs(st.Icond) > (this._IcondMax || 0)) {
      this._IcondMax = Math.abs(st.Icond);
      this.peakSnap = filamentProfile(s, st);
    }

    // --- детектор импульсов тока по |Icond| ---
    const Ic = Math.abs(st.Icond);
    const thr = 1e-5;                    // 10 мкА — заведомо выше тёмного тока
    if (!this._inPulse && Ic > thr) {
      this._inPulse = true;
      this._cur = { t0: st.t, tPeak: st.t, Ipeak: Ic, q: 0, UgapPeak: st.Ugap,
                    UgapAtStart: st.Ugap, UappAtStart: st.Uapp, ENpeak: st.maxEN, half: [] };
    }
    if (this._inPulse) {
      const c = this._cur;
      c.q += Ic * dt;
      if (Ic > c.Ipeak) { c.Ipeak = Ic; c.tPeak = st.t; c.UgapPeak = st.Ugap; c.ENpeak = st.maxEN; }
      c.half.push([st.t, Ic]);
      if (c.half.length > 200000) c.half.shift();
      if (Ic < thr) {
        this._inPulse = false;
        c.t1 = st.t;
        c.fwhm = fwhm(c.half, c.Ipeak);
        c.Jpeak = c.Ipeak / this.area;
        delete c.half;
        // напряжение горения = Ugap на пике импульса
        this.pulses.push(c);
        // отношение отрицательных ионов к электронам сразу после импульса
        this.negOverE.push({ t: st.t, ratio: negRatio(s) });
      }
    }

    // --- бухгалтерия по полупериодам (V6: перенесённый заряд) ---
    const hi = Math.floor(st.t / (0.5 * this.period));
    if (hi !== this._hIdx) {
      this._hIdx = hi;
      this.halves.push({ idx: hi, t0: st.t, Qstart: st.Q, absQcond: 0, Qcond: 0, pulses: 0,
                         sigmaMax: 0, UgapMax: 0 });
    }
    // V12: полный заряд системы (объём газа + оба sigma) обязан сохраняться ТОЖДЕСТВЕННО:
    // реакции идут через общий extent (A5), фото- и фоновый источник рождают ПАРЫ,
    // стенка лишь перекладывает заряд из газа в sigma. Единственный сток — floor-клип.
    if (this._qTot0 === undefined) { this._qTot0 = this.sim.charge().total; this._qScale = 0; }
    const qt = this.sim.charge();
    this._qTotLast = qt.total;
    const sc = Math.abs(qt.volume) + Math.abs(qt.surface);
    if (sc > this._qScale) this._qScale = sc;

    const H = this.halves[this.halves.length - 1];
    H.t1 = st.t; H.Qend = st.Q;
    H.absQcond += Math.abs(st.Icond) * dt;
    H.Qcond += st.Icond * dt;
    if (st.sigmaMax > H.sigmaMax) H.sigmaMax = st.sigmaMax;
    if (Math.abs(st.Ugap) > H.UgapMax) H.UgapMax = Math.abs(st.Ugap);

    if (this.samples.length === 0 || st.t - this.samples[this.samples.length - 1][0] > this.period / 4000) {
      this.samples.push([st.t, U, st.Ugap, st.Icond, st.Itot, st.Q, st.maxEN, st.sigmaMax, s.neMax, st.o3ppm]);
    }
  }

  /** Ёмкости из петли Лиссажу за указанный интервал времени. */
  lissajous(t0, t1) {
    const pts = this.samples.filter((p) => p[0] >= t0 && p[0] <= t1);
    if (pts.length < 20) return null;
    // Наклон ветвей: «разряд горит» (|Ugap| велик, |Icond| велик) => C_diel,
    // «разряд погашен» => C_cell. Разделяем по |Icond|.
    const Ithr = 1e-5;
    const off = pts.filter((p) => Math.abs(p[3]) < Ithr);
    const on = pts.filter((p) => Math.abs(p[3]) >= Ithr);
    return { Ccell: slope(off), Cdiel: slope(on), nOff: off.length, nOn: on.length,
             Umin: Math.min(...pts.map((p) => p[1])), Umax: Math.max(...pts.map((p) => p[1])),
             Qmin: Math.min(...pts.map((p) => p[5])), Qmax: Math.max(...pts.map((p) => p[5])) };
  }
}

/** Линейная регрессия Q(U) -> наклон (Ф). */
function slope(pts) {
  const n = pts.length;
  if (n < 5) return NaN;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (const p of pts) { const x = p[1], y = p[5]; sx += x; sy += y; sxx += x * x; sxy += x * y; }
  return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

function fwhm(pairs, Ipeak) {
  if (!pairs || pairs.length < 3) return NaN;
  const h = 0.5 * Ipeak;
  let a = null, b = null;
  for (const [t, v] of pairs) { if (v >= h) { if (a === null) a = t; b = t; } }
  return a === null ? NaN : b - a;
}

/**
 * Радиальный профиль филамента: n_e(r) = max по z в газе, и радиус на полувысоте
 * (FWHM/2) + «эффективный токовый радиус» sqrt(<r^2>) по весу n_e.
 */
function filamentProfile(s, st) {
  const S = s.poisson, nz = s.nz;
  const prof = new Float64Array(s.nr);
  for (let i = 0; i < s.nr; i++) {
    let m = 0;
    for (let j = s.JG0; j <= s.JG1; j++) { const v = s.n.e[i * nz + j]; if (v > m) m = v; }
    prof[i] = m;
  }
  const p0 = prof[0];
  let rHalf = NaN;
  for (let i = 1; i < s.nr; i++) {
    if (prof[i] < 0.5 * p0) {
      const f = (0.5 * p0 - prof[i - 1]) / (prof[i] - prof[i - 1]);
      rHalf = S.rc[i - 1] + f * (S.rc[i] - S.rc[i - 1]);
      break;
    }
  }
  // средний квадрат радиуса по объёмному весу n_e
  let w = 0, wr2 = 0;
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) {
      const q = s.n.e[i * nz + j] * S.Acell[i] * S.dz[j];
      w += q; wr2 += q * S.rc[i] * S.rc[i];
    }
  }
  return {
    t: st.t, Icond: st.Icond, Ugap: st.Ugap, maxEN: st.maxEN, neAxisMax: p0,
    rHalf_um: rHalf * 1e6, rRms_um: Math.sqrt(wr2 / Math.max(w, 1e-300)) * 1e6,
    profile_ne: Array.from(prof), r_um: Array.from(S.rc, (x) => x * 1e6),
  };
}

/** Отношение полного числа отрицательных ионов к числу электронов (V9). */
function negRatio(s) {
  const S = s.poisson, nz = s.nz;
  let ne = 0, nneg = 0;
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) {
      const k = i * nz + j, w = S.Acell[i] * S.dz[j];
      ne += s.n.e[k] * w;
      nneg += (s.n.Om[k] + s.n.O2m[k] + s.n.O3m[k]) * w;
    }
  }
  return { ne, nneg, ratio: nneg / Math.max(ne, 1e-300) };
}

// ──────────────────────────────────────────────────────────────── прогресс

function fmt(x, d = 3) {
  if (!Number.isFinite(x)) return String(x);
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return x.toExponential(d);
  return x.toFixed(d);
}

function progressLine(sim, st, k, t0wall, tEnd) {
  const el = (Date.now() - t0wall) / 1000;
  const frac = st.t / tEnd;
  const eta = frac > 1e-6 ? el * (1 / frac - 1) : NaN;
  return `[${String(k).padStart(9)}] t=${st.t.toExponential(4)} dt=${st.dt.toExponential(2)}`
    + ` (${sim.limiter}) U=${fmt(st.Uapp / 1e3, 2)}кВ Ug=${fmt(st.Ugap / 1e3, 2)}кВ`
    + ` I=${fmt(st.Itot * 1e3, 3)}мА E/N=${fmt(st.maxEN, 1)}Тд ne=${st.neMax.toExponential(2)}`
    + ` | ${fmt(el, 0)}с, ETA ${Number.isFinite(eta) ? fmt(eta, 0) : '?'}с, rej=${sim.rejects}`;
}

// ──────────────────────────────────────────────────────────────── main

export async function main(argv = process.argv.slice(2)) {
  const args = parseArgs(argv);
  if (args.help || args.h) { process.stderr.write(HELP); return 0; }

  const { preset, params } = buildParams(args);
  const periods = asNum(args.periods, 2);
  const outDir = args.out || path.join('data', preset);
  const doRecord = asBool(args.record, true);
  const progEvery = asNum(args.progress, 2000);
  const levels = (args.level || 'full,compact').split(',').map((s) => s.trim()).filter(Boolean);
  const maxSeconds = asNum(args.maxSeconds, Infinity);
  const maxSteps = asNum(args.maxSteps, Infinity);

  const sim = new DBD2D(params);
  const tEnd = periods * sim.period;

  process.stderr.write(
    `# preset=${preset} nr=${sim.nr} nz=${sim.nz} (gas ${sim.JG1 - sim.JG0 + 1})`
    + ` U0=${sim.p.U0kV}кВ f=${sim.p.freqKHz}кГц gamma=${sim.p.gammaIon}`
    + ` photo=${sim.p.photoModule ? `${+sim.p.photoIonization}${+sim.p.photoEmission}${+sim.p.photoDetachment}` : 'off'}`
    + ` d6=${sim.p.d6MassConvention} seedSpot=${sim.p.seedSpotAmp}\n`
    + `# tEnd=${tEnd.toExponential(3)} с, dz_wall=${(sim.poisson.dz[sim.JG0] * 1e6).toFixed(2)} мкм,`
    + ` dr0=${(sim.poisson.dr[0] * 1e6).toFixed(2)} мкм, C_cell=${(sim.Ccell * 1e15).toFixed(3)} фФ\n`,
  );

  let recs = null;
  if (doRecord) {
    fs.mkdirSync(path.dirname(outDir) || '.', { recursive: true });
    recs = createRecorders({
      dataDir: path.dirname(outDir) || '.',
      runId: path.basename(outDir),
      levels,
      seriesMaxRecords: asNum(args.seriesMax, undefined),
      maxBytes: asNum(args.maxBytes, undefined),
      // бюджет по уровням: полный — научный, компактный — для плеера (< 30 МБ)
      perLevel: {
        full: { maxBytes: asNum(args.maxBytesFull, 350e6) },
        compact: { maxBytes: asNum(args.maxBytesCompact, 24e6) },
      },
      dtFrameMin: asNum(args.dtFrameMin, undefined),
      dtFrameMax: asNum(args.dtFrameMax, undefined),
      params: { ...sim.p, preset, periods },
      meta: { generator: 'sim2d/run.mjs', node: process.version, startedAt: new Date().toISOString() },
      tEnd,
    });
  }

  const an = new Analyzer(sim, periods);
  const t0wall = Date.now();
  let k = 0, crashed = null;
  let lastProg = 0;

  // Мягкая остановка по сигналу: без неё убитый прогон теряет manifest.json
  // (индекс кадров живёт в памяти рекордера) и все записанные байты становятся мусором.
  let stopSignal = null;
  const ckptEvery = asNum(args.checkpointSec, 300);   // как часто фиксировать manifest
  let lastCkpt = Date.now();
  for (const sig of ['SIGINT', 'SIGTERM', 'SIGHUP']) {
    process.on(sig, () => {
      if (stopSignal) process.exit(130);      // второй сигнал — жёстко
      stopSignal = sig;
      process.stderr.write(`\n# получен ${sig}: закрываю прогон штатно (manifest будет записан)\n`);
    });
  }

  while (sim.t < tEnd) {
    let dt;
    try {
      dt = sim.step();
    } catch (e) {
      crashed = { step: k, t: sim.t, message: e.message };
      process.stderr.write(`\n!! АВАРИЯ на шаге ${k}, t=${sim.t.toExponential(6)} с:\n   ${e.message}\n`);
      break;
    }
    k++;
    const st = sim.state;
    an.step(st, dt);
    if (recs) recs.record(st, k === 1);
    if (progEvery > 0 && (k - lastProg >= progEvery)) {
      lastProg = k;
      process.stderr.write(progressLine(sim, st, k, t0wall, tEnd) + '\n');
      // промежуточная фиксация manifest (см. Recorder.checkpoint): длинный прогон
      // должен переживать kill без потери индекса кадров
      if (recs && Date.now() - lastCkpt > ckptEvery * 1e3) {
        lastCkpt = Date.now();
        try { recs.checkpoint({ validation: null, partial: true }); } catch (e) { /* не фатально */ }
      }
    }
    if (stopSignal) { process.stderr.write(`\n# остановка по сигналу ${stopSignal}\n`); break; }
    if (k >= maxSteps) { process.stderr.write(`\n# остановка: maxSteps=${maxSteps}\n`); break; }
    if ((Date.now() - t0wall) / 1000 > maxSeconds) {
      process.stderr.write(`\n# остановка: maxSeconds=${maxSeconds}\n`); break;
    }
  }

  const wall = (Date.now() - t0wall) / 1000;
  process.stderr.write(progressLine(sim, sim.state, k, t0wall, tEnd) + '  [конец]\n');

  const summary = buildSummary(sim, an, { preset, periods, tEnd, steps: k, wall, crashed, outDir });
  if (recs) { try { recs.close({ validation: summary, crashed }); } catch (e) { /* пустой прогон */ } }
  fs.mkdirSync(outDir, { recursive: true });
  fs.writeFileSync(path.join(outDir, 'summary.json'), JSON.stringify(summary, null, 2));
  process.stdout.write(JSON.stringify(summary, null, 2) + '\n');
  return crashed ? 2 : 0;
}

function buildSummary(sim, an, info) {
  const period = sim.period;
  const tEnd = Math.min(info.tEnd, sim.t);
  // последний ПОЛНЫЙ период (для Лиссажу и мощности)
  const tA = Math.max(0, tEnd - period), tB = tEnd;
  const inLast = an.samples.filter((p) => p[0] >= tA);
  let Wlast = 0;
  for (let i = 1; i < inLast.length; i++) {
    const dt = inLast[i][0] - inLast[i - 1][0];
    Wlast += 0.5 * (inLast[i][1] * inLast[i][4] + inLast[i - 1][1] * inLast[i - 1][4]) * dt;
  }
  const Plast = period > 0 ? Wlast / period : 0;
  const lis = an.lissajous(tA, tB);
  const areaCM2 = an.area * 1e4;

  const pulses = an.pulses.map((p) => ({
    t: p.t0, tPeak: p.tPeak, dur: p.t1 - p.t0, fwhm: p.fwhm,
    Ipeak: p.Ipeak, Jpeak: p.Jpeak, q: p.q, Ugap: p.UgapPeak, EN: p.ENpeak,
    UgapAtStart: p.UgapAtStart, UappAtStart: p.UappAtStart,
  }));

  return {
    run: info,
    geometry: {
      nr: sim.nr, nz: sim.nz, JG0: sim.JG0, JG1: sim.JG1,
      R_m: sim.poisson.R, area_m2: an.area, area_cm2: areaCM2,
      dzWall_um: sim.poisson.dz[sim.JG0] * 1e6,
      dEff_m: sim.dEff,
      C_cell_analytic_F: an.Ccell, C_diel_analytic_F: an.Cd, C_gas_analytic_F: an.Cg,
      C_cell_solver_F: sim.Ccell,
    },
    perf: { steps: info.steps, wall_s: info.wall, ms_per_step: (info.wall * 1e3) / Math.max(1, info.steps),
            rejects: sim.rejects, clipCount: sim.clipCount, qClip: sim.qClip },
    final: { t: sim.t, dt: sim.dt, limiter: sim.limiter, Uapp: sim.Uapp, Ugap: sim.Ugap,
             maxEN: an.maxENrun, neMax: an.neMaxRun, sigmaMax: an.sigmaMaxRun, o3ppm: sim.o3ppm },
    power: { W_lastPeriod_J: Wlast, P_lastPeriod_W: Plast,
             P_per_cm2_W: Plast / areaCM2, W_total_J: an.energy },
    lissajous: lis,
    halfPeriods: an.halves.map((h) => ({ ...h, dQ: h.Qend - h.Qstart })),
    pulses,
    filamentAtPeak: an.peakSnap || null,
    negRatioAfterPulse: an.negOverE.slice(-6),
    o3: { ppm_final: sim.o3ppm, n_O3_mean: (sim.o3ppm * 1e-6) * sim.N },
    charge: {
      Q_final: sim.Q, qClip: sim.qClip,
      qTotalStart: an._qTot0, qTotalEnd: an._qTotLast, qScaleMax: an._qScale,
      relDrift: an._qScale > 0 ? Math.abs(an._qTotLast - an._qTot0) / an._qScale : null,
    },
    samplesCount: an.samples.length,
  };
}

const isMain = process.argv[1] && path.resolve(process.argv[1]) === path.resolve(new URL(import.meta.url).pathname);
if (isMain) {
  main().then((c) => process.exit(c)).catch((e) => { console.error(e); process.exit(1); });
}

export default main;
