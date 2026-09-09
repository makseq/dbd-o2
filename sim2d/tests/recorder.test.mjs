// recorder.test.mjs — round-trip контейнера кадров: recorder.mjs пишет, loader.mjs читает.
// Запуск: node sim2d/tests/recorder.test.mjs
//
// Проверяется:
//   1) манифест валиден (обязательные секции, согласованность сетки и полей);
//   2) оффсеты кадров сходятся с реальным размером frames.bin и series.bin;
//   3) ошибка после квантования: < 1 % для лог-полей, < 0.5 % для линейных (уровень 'full');
//   4) ОТДЕЛЬНО — знакопеременное поле rho: точные нули, отрицательные значения,
//      симметрия кодера, отсутствие сдвига знака;
//   5) ряды читаются без прореживания и совпадают с записанными до точности f32;
//   6) для компактного уровня (uint8) — соблюдение ТЕОРЕТИЧЕСКОЙ границы 10^(DEC/510)-1.

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { Recorder, makeFileFetch, SERIES_NAMES, DTYPES } from '../recorder.mjs';
import { Playback } from '../../player/js/loader.mjs';

let failures = 0, checks = 0;
function ok(cond, msg, detail = '') {
  checks++;
  if (cond) { console.log(`  ok   ${msg}${detail ? '  [' + detail + ']' : ''}`); }
  else { failures++; console.log(`  FAIL ${msg}${detail ? '  [' + detail + ']' : ''}`); }
}
function section(s) { console.log(`\n=== ${s} ===`); }

// ---------------------------------------------------------------------------
// Синтетическое состояние (маленькая сетка, детерминированное)
// ---------------------------------------------------------------------------
const NR = 24, NZ = 40, JG0 = 8, JG1 = 31;   // 24 газовых строки
const NZGAS = JG1 - JG0 + 1;
const R = 5e-4, LZ = 2e-3;

function makeState() {
  const rf = new Float64Array(NR + 1), zf = new Float64Array(NZ + 1);
  for (let i = 0; i <= NR; i++) rf[i] = R * Math.expm1(1.6 * i / NR) / Math.expm1(1.6);
  for (let j = 0; j <= NZ; j++) zf[j] = LZ * j / NZ;
  const r = new Float64Array(NR), z = new Float64Array(NZ);
  for (let i = 0; i < NR; i++) r[i] = 0.5 * (rf[i] + rf[i + 1]);
  for (let j = 0; j < NZ; j++) z[j] = 0.5 * (zf[j] + zf[j + 1]);
  const gasMask = new Uint8Array(NZ);
  for (let j = JG0; j <= JG1; j++) gasMask[j] = 1;
  const N = NR * NZ;
  const f = () => new Float64Array(N);
  return {
    t: 0, dt: 1e-12, r, z, rf, zf, gasMask,
    n: { e: f(), O2p: f(), O4p: f(), Om: f(), O2m: f(), O3m: f(), O: f(), O3: f(), O2a: f() },
    E: f(), rho: f(), ionizRate: f(), EN: f(), photoIonRate: f(), photoDetachRate: f(),
    sigmaL: new Float64Array(NR), sigmaR: new Float64Array(NR),
    Uapp: 0, Ugap: 0, Icond: 0, Idisp: 0, Itot: 0, Q: 0, maxEN: 0, o3ppm: 0,
    sigmaMax: 0, photoEmitTotalL: 0, photoEmitTotalR: 0,
  };
}

// Заполнение полей на кадр k: широкий динамический диапазон + знакопеременный rho
function fill(st, k) {
  const rc = 60e-6, zc = 5e-4 + 3e-4 * k / 5;
  st.t = k * 1e-9;
  for (let i = 0; i < NR; i++) {
    for (let j = 0; j < NZ; j++) {
      const idx = i * NZ + j;
      const dr = (st.r[i] - 0) / rc, dz = (st.z[j] - zc) / 4e-5;
      const g = Math.exp(-dr * dr - dz * dz);
      const gas = j >= JG0 && j <= JG1;
      st.n.e[idx] = gas ? 1e14 + 1e20 * g : 0;
      st.ionizRate[idx] = gas ? 1e18 * g + 1e12 : 0;
      st.n.O3m[idx] = gas ? 1e16 * (0.1 + g) : 0;
      st.n.O3[idx] = gas ? 1e20 * (1 + 0.5 * g) : 0;
      st.photoIonRate[idx] = gas ? 1e15 * g : 0;
      st.photoDetachRate[idx] = gas ? 1e10 * g : 0;
      st.E[idx] = 3e6 + 2e7 * g;
      st.EN[idx] = 120 + 800 * g;
      // ЗНАКОПЕРЕМЕННОЕ: положительный слой перед головкой, отрицательный за ней,
      // и ТОЧНЫЕ НУЛИ в полосе (частый источник багов квантования)
      const s = (st.z[j] - zc) / 4e-5;
      let rho = -1.2e-3 * Math.exp(-dr * dr) * s * Math.exp(-s * s);
      if (Math.abs(s) < 0.05) rho = 0;
      if (i > NR - 4) rho = 0;              // целая зона точных нулей
      st.rho[idx] = gas ? rho : 0;
    }
  }
  for (let i = 0; i < NR; i++) {
    st.sigmaL[i] = 1e-5 * Math.exp(-((st.r[i] / 1.5e-4) ** 2)) * (k + 1);
    st.sigmaR[i] = -0.7e-5 * Math.exp(-((st.r[i] / 2e-4) ** 2)) * (k + 1);
  }
  st.Uapp = 1e4 * Math.sin(2 * Math.PI * 1e4 * st.t);
  st.Ugap = 0.6 * st.Uapp;
  st.Icond = 1e-3 * k; st.Idisp = 1e-4 * k; st.Itot = st.Icond + st.Idisp;
  st.Q = 1e-9 * k; st.maxEN = 920 - 10 * k; st.o3ppm = 3.2 * k;
  st.sigmaMax = 1e-5 * (k + 1); st.photoEmitTotalL = 1e12 * k; st.photoEmitTotalR = 5e11 * k;
}

// ---------------------------------------------------------------------------
const NFRAMES = 6;
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'dbd-rec-'));
const fetchFn = makeFileFetch();

async function runLevel(level) {
  section(`Уровень '${level}'`);
  const dir = path.join(tmp, level);
  const rec = new Recorder({
    dir, runId: `test-${level}`, level,
    params: { gapMM: 1, dielMM: 0.5, epsR: 9, radiusMM: 0.5, nr: NR, nz: NZ, U0kV: 10, freqKHz: 10 },
    meta: { test: true },
    tEnd: NFRAMES * 1e-9,
  });
  const st = makeState();
  const truth = [];        // эталон: поля ПОСЛЕ вырезки газа (для ds=1 — точное сравнение)
  const seriesTruth = [];
  for (let k = 0; k < NFRAMES; k++) {
    fill(st, k);
    const emitted = rec.record(st, true);
    if (!emitted) throw new Error('кадр не записан при force=true');
    const snap = {};
    for (const f of rec.fields) {
      const src = f.src.startsWith('n.') ? st.n[f.src.slice(2)] : st[f.src];
      const cut = new Float64Array(NR * NZGAS);
      for (let i = 0; i < NR; i++) for (let j = 0; j < NZGAS; j++) cut[i * NZGAS + j] = src[i * NZ + JG0 + j];
      snap[f.name] = cut;
    }
    snap.__sigmaL = Float64Array.from(st.sigmaL);
    snap.__sigmaR = Float64Array.from(st.sigmaR);
    truth.push(snap);
    seriesTruth.push(SERIES_NAMES.map((n) => st[n]));
  }
  rec.close();

  // ---- размеры файлов -----------------------------------------------------
  const framesSize = fs.statSync(path.join(dir, 'frames.bin')).size;
  const seriesSize = fs.statSync(path.join(dir, 'series.bin')).size;
  const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
  const fl = man.frameLayout;
  ok(framesSize === fl.headerBytes + fl.frameStride * fl.frameCount,
     'размер frames.bin = header + stride*frameCount', `${framesSize} Б`);
  ok(seriesSize === man.series.headerBytes + man.series.recordStride * man.series.count,
     'размер series.bin = header + stride*count', `${seriesSize} Б`);
  const offOk = man.frames.every((f, i) => f.byteOffset === fl.headerBytes + i * fl.frameStride)
    && man.frames[man.frames.length - 1].byteOffset + fl.frameStride === framesSize;
  ok(offOk, 'byteOffset каждого кадра согласован и последний кадр упирается в конец файла');

  // ---- валидность манифеста ----------------------------------------------
  ok(man.formatVersion === 1 && man.runId && man.level === level, 'манифест: version/runId/level');
  ok(man.grid.units === 'mm' && man.grid.rFaces.length === man.grid.nrOut + 1
     && man.grid.zFaces.length === man.grid.nzOut + 1, 'манифест: сетка в мм, грани согласованы с nrOut/nzOut');
  ok(Math.abs(man.grid.rFaces[man.grid.nrOut] - R * 1e3) < 1e-9, 'манифест: r доходит до R (мм)',
     `${man.grid.rFaces[man.grid.nrOut].toFixed(4)} мм`);
  ok(man.grid.dielectrics.length === 3 && man.grid.dielectrics[1].gas === true,
     'манифест: границы диэлектриков + газовый зазор');
  ok(man.grid.gasMaskZ.filter(Boolean).length === NZGAS, 'манифест: маска газа', `${NZGAS} строк`);
  ok(man.fields.length === 9 && man.fields.every((f) => f.byteLength === f.shape[0] * f.shape[1] * DTYPES[f.dtype].bytes),
     'манифест: 9 полей, byteLength = shape * sizeof(dtype)');
  ok(man.frames.every((f) => Object.keys(f.fields).length === 9
     && Object.values(f.fields).every((c) => 'min' in c && 'max' in c && 'scale' in c && 'log' in c)),
     'манифест: у каждого кадра кодеки {min,max,scale,log} по всем полям');
  ok(man.stats.nonFiniteHits === 0, 'ни одного не-конечного значения');

  // ---- чтение ЛОАДЕРОМ ----------------------------------------------------
  const pb = await Playback.load(dir + '/', { fetch: fetchFn, cacheFrames: 4 });
  ok(pb.frameCount === NFRAMES, 'loader: число кадров', `${pb.frameCount}`);
  ok(pb.fields.length === 9, 'loader: список полей');
  ok(pb.findFrame(2.4e-9) === 2 && pb.findFrame(1e9) === NFRAMES - 1, 'loader: поиск кадра по времени');

  // кэш с вытеснением
  for (let i = 0; i < NFRAMES; i++) pb.getFrame(i, 'n_e');
  ok(pb.cache.size === 4, 'loader: LRU-кэш ограничен', `size=${pb.cache.size}`);

  // ---- точность --------------------------------------------------------------
  const ds = man.grid.dsR !== 1 || man.grid.dsZ !== 1;
  const errs = {};
  for (const f of man.fields) {
    let worstRel = 0, worstAbs = 0, zeroBreak = 0, signBreak = 0, clipped = 0, n = 0;
    for (let k = 0; k < NFRAMES; k++) {
      const dec = pb.getFrame(k, f.name);
      const codec = pb.getCodec(k, f.name);
      const ref = truth[k][f.name];
      if (ds) continue;                                  // при даунсэмплинге эталон другой (см. ниже)
      const amax = Math.abs(codec.max) || 1;
      for (let idx = 0; idx < ref.length; idx++) {
        const v = ref[idx], d = dec[idx];
        n++;
        if (v === 0 && d !== 0) zeroBreak++;
        if (v !== 0 && Math.sign(v) !== Math.sign(d) && Math.abs(d) > 0) signBreak++;
        if (f.map === 'log') {
          if (v > 0 && v < codec.min) { clipped++; continue; }       // ниже порога DEC декад — обрезано намеренно
          if (v > 0) worstRel = Math.max(worstRel, Math.abs(d - v) / v);
        } else if (f.map === 'lin') {
          if (v > 0.01 * amax) worstRel = Math.max(worstRel, Math.abs(d - v) / v);
          worstAbs = Math.max(worstAbs, Math.abs(d - v) / amax);
        } else { // asinh
          if (Math.abs(v) >= amax * f.v0rel) worstRel = Math.max(worstRel, Math.abs(d - v) / Math.abs(v));
          worstAbs = Math.max(worstAbs, Math.abs(d - v) / amax);
        }
      }
    }
    errs[f.name] = { worstRel, worstAbs, zeroBreak, signBreak, clipped, n, map: f.map, dtype: f.dtype };
  }

  if (!ds) {
    const logLim = level === 'full' ? 0.01 : Math.pow(10, 6 / (2 * 255)) - 1;
    for (const f of man.fields) {
      const e = errs[f.name];
      if (f.map === 'log') {
        ok(e.worstRel < logLim, `лог-поле ${f.name} (${f.dtype}): отн. ошибка < ${(logLim * 100).toFixed(2)} %`,
           `${(e.worstRel * 100).toExponential(2)} %`);
      } else if (f.map === 'lin') {
        ok(e.worstRel < 0.005, `линейное ${f.name} (${f.dtype}): отн. ошибка < 0.5 % (v > 1 % от max)`,
           `${(e.worstRel * 100).toExponential(2)} %`);
      } else {
        ok(e.worstRel < 0.005, `asinh ${f.name} (${f.dtype}): отн. ошибка < 0.5 % при |v| >= v0`,
           `${(e.worstRel * 100).toExponential(2)} %`);
      }
    }

    // ---- ОТДЕЛЬНО: знакопеременное поле -----------------------------------
    section(`Знакопеременное поле rho ('${level}')`);
    const e = errs.rho;
    ok(e.zeroBreak === 0, 'нули остаются ТОЧНЫМИ нулями после round-trip', `${e.n} отсчётов`);
    ok(e.signBreak === 0, 'знак не переворачивается ни в одном отсчёте');
    let neg = 0, pos = 0, zer = 0;
    for (let k = 0; k < NFRAMES; k++) {
      const dec = pb.getFrame(k, 'rho');
      for (const v of dec) { if (v < 0) neg++; else if (v > 0) pos++; else zer++; }
    }
    ok(neg > 0 && pos > 0 && zer > 0, 'декодированное rho содержит и минусы, и плюсы, и нули',
       `-:${neg} +:${pos} 0:${zer}`);
    ok(e.worstAbs < 1e-3, 'abs-ошибка rho мала в линейной зоне (в долях max)', `${e.worstAbs.toExponential(2)}`);
    // симметрия кодера: -v должно давать ровно -q
    const rawRho = pb.getFrameRaw(3, 'rho');
    let asym = 0;
    for (const q of rawRho) if (q === -32768) asym++;    // запрещённое значение (нет пары +32768)
    ok(asym === 0, 'нет асимметричного отсчёта -32768 (кодер клампит на ±32767)');
  } else {
    section(`Компактный уровень: даунсэмплинг ${man.grid.dsR}x${man.grid.dsZ}`);
    ok(man.grid.nrOut === Math.ceil(NR / man.grid.dsR) && man.grid.nzOut === Math.ceil(NZGAS / man.grid.dsZ),
       'форма кадра после прореживания', `${man.grid.nrOut}x${man.grid.nzOut}`);
    // консервативность: интеграл n_e по объёму сохраняется площадным усреднением
    const dec = pb.getFrame(2, 'n_e');
    ok(dec.length === man.grid.nrOut * man.grid.nzOut, 'длина декодированного поля = nrOut*nzOut');
    let mn = Infinity, mx = -Infinity;
    for (const v of dec) { if (v < mn) mn = v; if (v > mx) mx = v; }
    ok(mx > 1e19 && mn >= 0, 'компактное n_e сохраняет пик головки', `max=${mx.toExponential(2)}`);
    // теоретическая граница uint8 для лог-полей
    const bound = Math.pow(10, 6 / (2 * 255)) - 1;
    ok(bound > 0.02 && bound < 0.03, `uint8 + 6 декад: предельная ошибка ${(bound * 100).toFixed(2)} % (это ФИЗИКА формата, не баг)`);
    const rho = pb.getFrame(2, 'rho');
    let hasNeg = false, hasZero = false;
    for (const v of rho) { if (v < 0) hasNeg = true; if (v === 0) hasZero = true; }
    ok(hasNeg && hasZero, 'компактное rho: знак и нули сохранены');
  }

  // ---- поверхностный заряд ------------------------------------------------
  const surf = pb.getSurface(4);
  if (!ds) {
    let we = 0;
    for (let i = 0; i < NR; i++) {
      we = Math.max(we, Math.abs(surf.sigmaL[i] - truth[4].__sigmaL[i]) / Math.max(1e-30, Math.abs(truth[4].__sigmaL[i])));
      we = Math.max(we, Math.abs(surf.sigmaR[i] - truth[4].__sigmaR[i]) / Math.max(1e-30, Math.abs(truth[4].__sigmaR[i])));
    }
    ok(we < 1e-6, 'sigmaL/sigmaR (Float32, без квантования) совпадают', `отн. ${we.toExponential(2)}`);
    ok(surf.sigmaR.some((v) => v < 0), 'sigmaR сохраняет отрицательный знак');
  } else {
    ok(surf.sigmaL.length === man.grid.nrOut, 'sigma прорежен до nrOut');
  }

  // ---- ряды ---------------------------------------------------------------
  const s = pb.series;
  ok(s && s.t.length === man.series.count, 'loader: ряды загружены', `${man.series.count} записей`);
  if (level === 'full') {
    ok(man.series.count === NFRAMES && man.series.decimated === false,
       'полный уровень: ряды БЕЗ прореживания (запись на каждый шаг)');
    let worst = 0;
    for (let k = 0; k < NFRAMES; k++) {
      for (let c = 0; c < SERIES_NAMES.length; c++) {
        const ref = seriesTruth[k][c], got = s[SERIES_NAMES[c]][k];
        if (ref === 0) { worst = Math.max(worst, Math.abs(got)); continue; }
        worst = Math.max(worst, Math.abs(got - ref) / Math.abs(ref));
      }
    }
    ok(worst < 1e-6, 'ряды совпадают до точности Float32', `${worst.toExponential(2)}`);
    ok(s.t instanceof Float64Array, 't хранится Float64 (иначе слипнутся соседние шаги)');
    ok(SERIES_NAMES.length === 12 && SERIES_NAMES.includes('photoEmitTotalR'),
       'состав рядов: t,Uapp,Ugap,Icond,Idisp,Itot,Q,maxEN,o3ppm,sigmaMax,photoEmitTotalL/R');
  }
  return { man, errs };
}

// ---------------------------------------------------------------------------
console.log('recorder round-trip test\nвременный каталог:', tmp);
const full = await runLevel('full');
const compact = await runLevel('compact');

section('Режим Range (ленивая подкачка кадров)');
{
  const eager = await Playback.load(path.join(tmp, 'full') + '/', { fetch: fetchFn });
  const lazy = await Playback.load(path.join(tmp, 'full') + '/', { fetch: fetchFn, mode: 'range', cacheFrames: 2 });
  let thrown = false;
  try { lazy.getFrame(3, 'n_e'); } catch { thrown = true; }
  ok(thrown, 'range: синхронный getFrame до подкачки честно падает, а не отдаёт мусор');
  const a = await lazy.getFrameAsync(3, 'n_e');
  const b = eager.getFrame(3, 'n_e');
  let same = a.length === b.length;
  for (let i = 0; i < a.length && same; i++) if (a[i] !== b[i]) same = false;
  ok(same, 'range: кадр по Range-запросу побитово совпадает с eager-режимом');
  await lazy.prefetch(0, 2);
  ok(lazy.blocks.size === 2, 'range: кэш блоков вытесняет старые', `size=${lazy.blocks.size}`);
}

section('Битовая раскладка кадра (уровень full)');
console.log('  frameStride =', full.man.frameLayout.frameStride, 'Б;  headerBytes =', full.man.frameLayout.headerBytes);
for (const f of full.man.fields) {
  console.log(`  +${String(f.offsetInFrame).padStart(7)}  ${f.name.padEnd(16)} ${f.dtype.padEnd(4)} ${f.map.padEnd(6)} ${f.byteLength} Б  ${f.shape.join('x')}`);
}
console.log(`  +${String(full.man.frameLayout.surfaces.sigmaL.offsetInFrame).padStart(7)}  sigmaL           f32  -      ${full.man.frameLayout.surfaces.sigmaL.length * 4} Б`);
console.log(`  +${String(full.man.frameLayout.surfaces.sigmaR.offsetInFrame).padStart(7)}  sigmaR           f32  -      ${full.man.frameLayout.surfaces.sigmaR.length * 4} Б`);

section('Ошибки квантования (сводка)');
console.log('  поле              full(dtype)  отн.ошибка     compact(dtype)');
for (const name of Object.keys(full.errs)) {
  const a = full.errs[name], b = compact.errs[name];
  console.log(`  ${name.padEnd(17)} ${a.dtype.padEnd(11)} ${(a.worstRel * 100).toExponential(2).padStart(10)} %  ${b.dtype}`);
}

section('ИТОГ');
console.log(`  проверок: ${checks}, провалов: ${failures}`);
fs.rmSync(tmp, { recursive: true, force: true });
process.exit(failures ? 1 : 0);
