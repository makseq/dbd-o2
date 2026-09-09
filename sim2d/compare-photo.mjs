#!/usr/bin/env node
// compare-photo.mjs — прямое сравнение двух прогонов (фото ВКЛ / ВЫКЛ) по series.bin.
// Печатает то, что нужно для раздела «Влияние фотопроцессов» в VALIDATION.md.
//
//   node sim2d/compare-photo.mjs data/run-default data/run-nophoto

import fs from 'node:fs';
import path from 'node:path';

function load(dir) {
  const man = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'), 'utf8'));
  const s = man.series;
  const buf = fs.readFileSync(path.join(dir, s.file));
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const idx = Object.fromEntries(s.names.map((n, i) => [n, i]));
  const rows = [];
  for (let k = 0; k < s.count; k++) {
    const o = s.headerBytes + k * s.recordStride;
    const g = (n) => {
      const i = idx[n];
      if (i === undefined) return NaN;
      return s.dtypes[i] === 'f64' ? dv.getFloat64(o + s.offsets[i], true) : dv.getFloat32(o + s.offsets[i], true);
    };
    rows.push({ t: g('t'), Uapp: g('Uapp'), Ugap: g('Ugap'), Icond: g('Icond'), Itot: g('Itot'),
                Q: g('Q'), maxEN: g('maxEN'), o3ppm: g('o3ppm'), sigmaMax: g('sigmaMax') });
  }
  const sum = JSON.parse(fs.readFileSync(path.join(dir, 'summary.json'), 'utf8'));
  return { dir, man, rows, sum };
}

/** Момент зажигания: первое превышение |Icond| порога thr. */
function ignition(rows, thr) {
  for (let k = 1; k < rows.length; k++) {
    if (Math.abs(rows[k].Icond) > thr) {
      const a = rows[k - 1], b = rows[k];
      const f = (thr - Math.abs(a.Icond)) / (Math.abs(b.Icond) - Math.abs(a.Icond));
      return { t: a.t + f * (b.t - a.t), Ugap: a.Ugap, Uapp: a.Uapp, k };
    }
  }
  return null;
}

function stats(r) {
  const rows = r.rows;
  let Ipk = 0, tPk = NaN, W = 0, ENmax = 0, sigMax = 0;
  for (let k = 0; k < rows.length; k++) {
    const a = rows[k];
    if (Math.abs(a.Icond) > Ipk) { Ipk = Math.abs(a.Icond); tPk = a.t; }
    if (a.maxEN > ENmax) ENmax = a.maxEN;
    if (a.sigmaMax > sigMax) sigMax = a.sigmaMax;
    if (k > 0) W += 0.5 * (a.Uapp * a.Itot + rows[k - 1].Uapp * rows[k - 1].Itot) * (a.t - rows[k - 1].t);
  }
  const fil = r.sum.filamentAtPeak || {};
  return {
    tEnd: rows[rows.length - 1].t, Ipk, tPk, W, ENmax, sigMax,
    o3: rows[rows.length - 1].o3ppm,
    ign1: ignition(rows, 1e-3), ign10: ignition(rows, 1e-2),
    rHalf: fil.rHalf_um, rRms: fil.rRms_um, neAxis: fil.neAxisMax,
    UgapAtPeak: fil.Ugap, crashed: r.sum.run.crashed ? r.sum.run.crashed.t : null,
  };
}

const [dA, dB] = process.argv.slice(2, 4);
const A = load(dA), B = load(dB);
const sa = stats(A), sb = stats(B);
const rel = (x, y) => (Number.isFinite(x) && Number.isFinite(y) && y !== 0 ? ((x - y) / y) * 100 : NaN);
const F = (x, d = 4) => (Number.isFinite(x) ? (Math.abs(x) >= 1e5 || (x !== 0 && Math.abs(x) < 1e-3) ? x.toExponential(d) : x.toFixed(d)) : '—');

const rowsOut = [
  ['момент зажигания (|I_cond| > 1 мА), с', sa.ign1 && sa.ign1.t, sb.ign1 && sb.ign1.t],
  ['U_gap в этот момент, В', sa.ign1 && sa.ign1.Ugap, sb.ign1 && sb.ign1.Ugap],
  ['момент |I_cond| > 10 мА, с', sa.ign10 && sa.ign10.t, sb.ign10 && sb.ign10.t],
  ['I_peak (|I_cond|), А', sa.Ipk, sb.Ipk],
  ['t(I_peak), с', sa.tPk, sb.tPk],
  ['U_gap в пике тока, В', sa.UgapAtPeak, sb.UgapAtPeak],
  ['радиус филамента r½, мкм', sa.rHalf, sb.rHalf],
  ['sqrt(<r²>) по весу n_e, мкм', sa.rRms, sb.rRms],
  ['n_e на оси в пике, м⁻³', sa.neAxis, sb.neAxis],
  ['max E/N за прогон, Тд', sa.ENmax, sb.ENmax],
  ['max |σ|, Кл/м²', sa.sigMax, sb.sigMax],
  ['энергия ∫U·I dt до конца, Дж', sa.W, sb.W],
  ['O₃ в конце, ppm', sa.o3, sb.o3],
  ['конец прогона t, с', sa.tEnd, sb.tEnd],
];

console.log(`# ${path.basename(dA)} (A) против ${path.basename(dB)} (B)`);
console.log(`A: photo=${+A.man.params.photoIonization}${+A.man.params.photoEmission}${+A.man.params.photoDetachment}`
  + `  B: photo=${+B.man.params.photoIonization}${+B.man.params.photoEmission}${+B.man.params.photoDetachment}`
  + `  сетка ${A.man.grid.nr}x${A.man.grid.nz}, U0=${A.man.params.U0kV} кВ, f=${A.man.params.freqKHz} кГц`);
console.log('');
console.log('| Показатель | A (фото ВКЛ) | B (фото ВЫКЛ) | A/B − 1 |');
console.log('|---|---|---|---|');
for (const [name, a, b] of rowsOut) {
  console.log(`| ${name} | ${F(a)} | ${F(b)} | ${Number.isFinite(rel(a, b)) ? rel(a, b).toFixed(1) + ' %' : '—'} |`);
}
console.log('');
console.log(`авария A: ${sa.crashed ? sa.crashed.toExponential(4) + ' с' : 'нет'}; авария B: ${sb.crashed ? sb.crashed.toExponential(4) + ' с' : 'нет'}`);
