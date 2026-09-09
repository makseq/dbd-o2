#!/usr/bin/env node
// summarize.mjs — сводка по всем прогонам в data/: читает <run>/manifest.json,
// frames.bin и series.bin и печатает таблицу приёмки.
//
// Ключевой критерий задания: ЕСТЬ ЛИ ФИЛАМЕНТ — отношение n_e на оси к n_e на краю
// в момент пика тока > 3. Считается по кадру, ближайшему к пику |Icond| из series.bin.
//
//   node sim2d/summarize.mjs [dataDir]        (по умолчанию data/)

import fs from 'node:fs';
import path from 'node:path';

const dataDir = process.argv[2] || 'data';

/** Обратное отображение кодека кадра (см. player/js/loader.mjs). */
function decode(codec, q, map, decades) {
  if (!codec || !(codec.scale > 0)) return 0;
  if (map === 'log') return q === 0 ? 0 : Math.pow(10, Math.log10(codec.min) + q / codec.scale);
  if (map === 'asinh') return codec.v0 * Math.sinh(q / codec.scale);
  return q / codec.scale;
}

function readSeries(dir, man) {
  const s = man.series;
  if (!s || !s.count) return null;
  const buf = fs.readFileSync(path.join(dir, s.file));
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const idx = Object.fromEntries(s.names.map((n, i) => [n, i]));
  const rec = (k, name) => {
    const i = idx[name];
    const off = s.headerBytes + k * s.recordStride + s.offsets[i];
    return s.dtypes[i] === 'f64' ? dv.getFloat64(off, true) : dv.getFloat32(off, true);
  };
  return { count: s.count, get: rec, has: (n) => n in idx };
}

/** Профиль n_e(r) = max по z для кадра с индексом fi. */
function frameNeProfile(dir, man, fi) {
  const fl = man.frameLayout;
  const fld = man.fields.find((f) => f.name === 'n_e');
  if (!fld || !fl) return null;
  const fr = man.frames[fi];
  const [nr, nz] = fld.shape;
  const fd = fs.openSync(path.join(dir, fl.file), 'r');
  const buf = Buffer.alloc(fld.byteLength);
  fs.readSync(fd, buf, 0, fld.byteLength, fr.byteOffset + fld.offsetInFrame);
  fs.closeSync(fd);
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const codec = fr.fields.n_e;
  const prof = new Float64Array(nr);
  const w = fld.dtype === 'u8' || fld.dtype === 'i8' ? 1 : 2;
  const rd = fld.dtype === 'u8' ? (o) => dv.getUint8(o)
    : fld.dtype === 'i8' ? (o) => dv.getInt8(o)
      : fld.dtype === 'i16' ? (o) => dv.getInt16(o, true) : (o) => dv.getUint16(o, true);
  for (let i = 0; i < nr; i++) {
    let m = 0;
    for (let j = 0; j < nz; j++) {
      const q = rd((i * nz + j) * w);
      const v = decode(codec, q, fld.map, fld.decades);
      if (v > m) m = v;
    }
    prof[i] = m;
  }
  return prof;
}

function du(dir) {
  let s = 0;
  for (const f of fs.readdirSync(dir)) {
    const st = fs.statSync(path.join(dir, f));
    if (st.isFile()) s += st.size;
  }
  return s;
}

const rows = [];
for (const run of fs.readdirSync(dataDir).sort()) {
  const dir = path.join(dataDir, run);
  const mf = path.join(dir, 'manifest.json');
  if (!fs.existsSync(mf)) continue;
  const man = JSON.parse(fs.readFileSync(mf, 'utf8'));
  const st = man.stats || {};
  const V = man.validation || {};
  const ser = readSeries(dir, man);

  // пик тока и мощность по рядам
  let Ipeak = 0, kPeak = -1, tPeak = NaN, ENmax = 0, o3 = 0, W = 0, tPrev = null, uPrev = 0, iPrev = 0;
  if (ser) {
    for (let k = 0; k < ser.count; k++) {
      const t = ser.get(k, 't'), Ic = Math.abs(ser.get(k, 'Icond'));
      if (Ic > Ipeak) { Ipeak = Ic; kPeak = k; tPeak = t; }
      const en = ser.get(k, 'maxEN'); if (en > ENmax) ENmax = en;
      const o = ser.get(k, 'o3ppm'); if (o > o3) o3 = o;
      const U = ser.get(k, 'Uapp'), I = ser.get(k, 'Itot');
      if (tPrev !== null) W += 0.5 * (U * I + uPrev * iPrev) * (t - tPrev);
      tPrev = t; uPrev = U; iPrev = I;
    }
  }
  const tSpan = st.tLast - st.tFirst;
  const P = tSpan > 0 ? W / tSpan : NaN;   // средняя мощность по СЧИТАННОМУ интервалу

  // филамент: кадр, ближайший к моменту пика тока
  let ratio = NaN, rHalfUm = NaN, neAxis = NaN;
  if (man.frames && man.frames.length && Number.isFinite(tPeak)) {
    let best = 0, bd = Infinity;
    for (let i = 0; i < man.frames.length; i++) {
      const d = Math.abs(man.frames[i].t - tPeak);
      if (d < bd) { bd = d; best = i; }
    }
    const prof = frameNeProfile(dir, man, best);
    if (prof && prof.length) {
      neAxis = prof[0];
      ratio = prof[0] / Math.max(prof[prof.length - 1], 1e-300);
      const rc = man.grid.rCenters;
      for (let i = 1; i < prof.length; i++) {
        if (prof[i] < 0.5 * prof[0]) {
          const f = (0.5 * prof[0] - prof[i - 1]) / (prof[i] - prof[i - 1]);
          rHalfUm = (rc[i - 1] + f * (rc[i] - rc[i - 1])) * 1e3;   // rCenters в мм
          break;
        }
      }
    }
  }

  // Кадры прорежены под бюджет, поэтому ближайший к пику кадр может быть в стороне.
  // Если рядом лежит summary.json — там снимок n_e(r) взят ТОЧНО в шаге максимума
  // |I_cond| (Analyzer.peakSnap), он и есть авторитетный источник для критерия.
  const sumPath = path.join(dir, 'summary.json');
  if (fs.existsSync(sumPath)) {
    try {
      const sm = JSON.parse(fs.readFileSync(sumPath, 'utf8'));
      const f = sm.filamentAtPeak;
      if (f && f.profile_ne && f.profile_ne.length) {
        const p = f.profile_ne;
        ratio = p[0] / Math.max(p[p.length - 1], 1e-300);
        rHalfUm = f.rHalf_um;
        neAxis = f.neAxisMax;
      }
    } catch (e) { /* summary может отсутствовать/быть битым — тогда остаётся оценка по кадру */ }
  }

  rows.push({
    run, level: man.level, nr: man.grid.nr, nz: man.grid.nz,
    U0kV: man.params.U0kV, fkHz: man.params.freqKHz,
    photo: `${+man.params.photoIonization}${+man.params.photoEmission}${+man.params.photoDetachment}`,
    frames: st.frames, bytes: du(dir),
    t0: st.tFirst, t1: st.tLast, periods: tSpan * man.params.freqKHz * 1e3,
    steps: st.steps, Ipeak, tPeak, P, ENmax, o3,
    ratio, rHalfUm, neAxis,
    crashed: man.crashed ? man.crashed.t : null,
    pulses: (V.pulses || []).length,
  });
}

const F = (x, d = 3) => (Number.isFinite(x) ? (Math.abs(x) >= 1e5 || (x !== 0 && Math.abs(x) < 1e-3) ? x.toExponential(d) : x.toFixed(d)) : '—');
const MB = (b) => (b / 1e6).toFixed(1);

console.log('# Сводка прогонов в', path.resolve(dataDir));
console.log('');
const hdr = ['run', 'lvl', 'сетка', 'U0', 'f', 'photo', 'кадры', 'МБ', 't, мкс', 'пер.',
  'I_peak, А', 'P, Вт', 'E/N max', 'O3 ppm', 'n_e(0)/n_e(R)', 'r½, мкм', 'филамент'];
console.log(hdr.join(' | '));
console.log(hdr.map(() => '---').join(' | '));
for (const r of rows) {
  console.log([
    r.run, r.level, `${r.nr}x${r.nz}`, r.U0kV, r.fkHz, r.photo, r.frames, MB(r.bytes),
    `${F(r.t0 * 1e6, 2)}…${F(r.t1 * 1e6, 2)}`, F(r.periods, 3),
    F(r.Ipeak, 4), F(r.P, 4), F(r.ENmax, 0), F(r.o3, 2),
    Number.isFinite(r.ratio) ? F(r.ratio, 2) : '∞', F(r.rHalfUm, 1), r.ratio > 3 ? 'ДА' : 'нет',
  ].join(' | '));
}
console.log('');
for (const r of rows) if (r.crashed) console.log(`! ${r.run}: аварийная остановка на t = ${r.crashed.toExponential(4)} с`);
const anyFil = rows.some((r) => r.ratio > 3);
console.log('');
console.log(anyFil ? 'ФИЛАМЕНТ ЕСТЬ хотя бы в одном прогоне (n_e ось/край > 3).'
  : 'ФИЛАМЕНТА НЕТ НИ В ОДНОМ ПРОГОНЕ — см. ERRATA B6 (затравочное пятно / сетка).');
console.log(`Суммарный объём data/: ${MB(rows.reduce((a, r) => a + r.bytes, 0))} МБ`);
