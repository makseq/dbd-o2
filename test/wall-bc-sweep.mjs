// wall-bc-sweep.mjs — развёртки по пристеночному ГУ (ГУ Хагелаара, src/solver.js
// _wallCoeffs) и по амплитуде. Печатает по каждому варианту: досчитал/сорвался,
// момент срыва, |sigma|max, max E/N, max n_e, ёмкости из Лиссажу, озон.
//
// Запуск:
//     node test/wall-bc-sweep.mjs [периодов] [режим] [набор]
//   набор: refl (по reflE), u0 (по амплитуде), bc (hagelaar vs legacy), all
//
// Потолки физичности — те же, что в test/run-period.mjs (E/N < 5000 Td, n_e < N).

import { DBDSolver } from '../src/solver.js';

const EN_CEIL = 5000;
const NE_CEIL = 2.4463e25;

const NPER = Number(process.argv[2] || 2);
const MODE = process.argv[3] || 'demo';
const SET = process.argv[4] || 'all';

function run(tag, params) {
  const t0 = Date.now();
  const solver = new DBDSolver(Object.assign({
    gapMM: 1.0, dielMM1: 0.5, dielMM2: 0.5, epsR: 9, areaCM2: 1.0,
    freqKHz: 10, mode: MODE,
  }, params));
  const period = solver.period;
  const tEnd = NPER * period;
  let diverged = null, sigMax = 0, enMax = 0, neMax = 0, ipk = 0, steps = 0;
  while (solver.t < tEnd) {
    solver.step();
    steps++;
    const s = solver.state;
    if (Math.abs(s.sigmaL) > sigMax) sigMax = Math.abs(s.sigmaL);
    if (Math.abs(s.sigmaR) > sigMax) sigMax = Math.abs(s.sigmaR);
    if (s.maxEN > enMax) enMax = s.maxEN;
    let ne = 0;
    for (let i = 0; i < s.n.e.length; i++) if (s.n.e[i] > ne) ne = s.n.e[i];
    if (ne > neMax) neMax = ne;
    if (Math.abs(s.current) > ipk) ipk = Math.abs(s.current);
    if (!(s.maxEN < EN_CEIL)) { diverged = `E/N=${s.maxEN.toFixed(0)}Td`; break; }
    if (!(ne < NE_CEIL)) { diverged = `n_e=${ne.toExponential(2)}`; break; }
  }
  const ps = solver.periodStats, s = solver.state;
  const wall = (Date.now() - t0) / 1000;
  return {
    tag, diverged, t: solver.t, steps, sigMax, enMax, neMax, ipk,
    sigL: s.sigmaL, sigR: s.sigmaR,
    Ccell: ps.Ccell, Cdiel: ps.CdielIdentity, CdielSlope: ps.Cdiel,
    o3: gasAvg(s, s.n.O3), o3ppm: s.o3ppm, power: ps.powerW, wall,
  };
}

function gasAvg(s, arr) {
  let sum = 0, w = 0;
  for (let i = 0; i < s.x.length; i++) {
    if (!s.gasMask[i]) continue;
    const dx = s.xFaces[i + 1] - s.xFaces[i];
    sum += arr[i] * dx; w += dx;
  }
  return w > 0 ? sum / w : 0;
}

function table(title, rows) {
  console.log('');
  console.log(`=== ${title} ${'='.repeat(Math.max(0, 62 - title.length))}`);
  console.log('вариант          | исход           | t, мкс | |σ|max Кл/м² | E/N max | n_e max | C_cell пФ | C_diel пФ | O3, м^-3 | P, Вт');
  for (const r of rows) {
    console.log(
      `${r.tag.padEnd(16)} | ${(r.diverged ? 'СРЫВ ' + r.diverged : 'досчитал').padEnd(15)} |` +
      ` ${(r.t * 1e6).toFixed(3).padStart(6)} | ${r.sigMax.toExponential(3).padStart(12)} |` +
      ` ${r.enMax.toFixed(0).padStart(7)} | ${r.neMax.toExponential(1).padStart(7)} |` +
      ` ${(r.Ccell * 1e12).toFixed(3).padStart(9)} | ${(r.Cdiel * 1e12).toFixed(2).padStart(9)} |` +
      ` ${r.o3.toExponential(2).padStart(8)} | ${r.power.toFixed(2).padStart(5)}`);
  }
}

const U0 = Number(process.env.U0KV || 10);
console.log(`развёртка: ${NPER} период(ов), режим ${MODE}`);

if (SET === 'bc' || SET === 'all') {
  const rows = [];
  for (const bc of ['legacy', 'hagelaar']) {
    rows.push(run(bc, { U0kV: U0, wallBC: bc }));
    console.log(`  ...${bc} готов (${rows[rows.length - 1].wall.toFixed(0)} с)`);
  }
  table(`ГУ: legacy vs Хагелаар, U0 = ${U0} кВ`, rows);
}

if (SET === 'refl' || SET === 'all') {
  const rows = [];
  for (const r of [0, 0.1, 0.2]) {
    rows.push(run(`reflE=${r}`, { U0kV: U0, reflE: r }));
    console.log(`  ...reflE=${r} готов (${rows[rows.length - 1].wall.toFixed(0)} с)`);
  }
  table(`развёртка по reflE, U0 = ${U0} кВ`, rows);
}

if (SET === 'u0' || SET === 'all') {
  const rows = [];
  for (const u of [6, 8, 10, 12]) {
    rows.push(run(`U0=${u} кВ`, { U0kV: u }));
    console.log(`  ...U0=${u} готов (${rows[rows.length - 1].wall.toFixed(0)} с)`);
  }
  table('развёртка по U0 (ГУ Хагелаара, reflE=0)', rows);
}
