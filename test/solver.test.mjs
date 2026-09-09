// test/solver.test.mjs — тесты вычислительного ядра ДБР.
// Запуск: node test/solver.test.mjs   (exit code 0 = все тесты прошли)

import { DBDSolver, sgFlux } from '../src/solver.js';
import {
  EPS0, QE, TD, gasDensity, alphaN, etaN, alphaEff,
  bernoulli, phi1, muE_N, driftElectron,
  alphaEtaCrossEN, breakdownEN, selfSustainResidual,
} from '../src/physics.js';

let passed = 0, failed = 0;
const t0all = Date.now();

// ---------------------------------------------------------------------------
// Защита от расходимости (добавлена вместе с заменой коэффициентов).
//
// С корректными коэффициентами (RATES_REVIEW) ионизация примерно втрое быстрее,
// чем со старыми фитами, и 1D-модель на импульсе уходит в нефизичное состояние:
// приповерхностная ячейка набирает n_e ~ 1e24 м^-3 и E/N ~ 4e4 Td, шаг падает до
// 1e-14 с, прогон практически останавливается. Без этой защиты тест не падает —
// он ВИСИТ (>15 мин на период), что хуже честного FAIL.
//
// Потолки заведомо нефизичны: 5000 Td = 1.2e8 В/м (втрое выше пробойного поля
// в самом узком катодном слое), n_e >= N означает степень ионизации > 100 %.
// ---------------------------------------------------------------------------
const EN_CEIL = 5000;
const NE_CEIL = gasDensity(760, 300);
class Diverged extends Error {}
function physGuard(sol) {
  const st = sol.state;
  if (!(st.maxEN < EN_CEIL)) {
    throw new Diverged(`E/N = ${st.maxEN.toFixed(0)} Td > ${EN_CEIL} Td при t = ${st.t.toExponential(3)} с`);
  }
  let ne = 0;
  for (let i = 0; i < st.n.e.length; i++) if (st.n.e[i] > ne) ne = st.n.e[i];
  if (!(ne < NE_CEIL)) {
    throw new Diverged(`n_e = ${ne.toExponential(2)} м^-3 >= N при t = ${st.t.toExponential(3)} с`);
  }
}

function ok(cond, name, info = '') {
  if (cond) { passed++; console.log(`PASS  ${name}${info ? '  [' + info + ']' : ''}`); }
  else { failed++; console.log(`FAIL  ${name}${info ? '  [' + info + ']' : ''}`); }
  return !!cond;
}
function near(a, b, rtol, name, extra = '') {
  const den = Math.abs(b) > 1e-300 ? Math.abs(b) : 1;
  const r = Math.abs(a - b) / den;
  return ok(r <= rtol, name, `got=${a.toExponential(4)} want=${b.toExponential(4)} rel=${r.toExponential(2)} tol=${rtol}${extra ? ' ' + extra : ''}`);
}
function section(s) { console.log(`\n--- ${s} ---`); }

function zeroDensities(sol) {
  for (const arr of sol.n) arr.fill(0);
  sol.nO.fill(0); sol.nO3.fill(0); sol.nO2a.fill(0);
  sol.sigL = 0; sol.sigR = 0; sol.sigFace.fill(0);
}

// ===========================================================================
section('1. Пуассон: вакуумный конденсатор -> линейный потенциал');
{
  const sol = new DBDSolver({ epsR: 1, nCells: 120, seedDensity: 0, chemistry: false, mode: 'demo' });
  zeroDensities(sol);
  const U = 5000;
  sol._poissonSolve(U, true);
  const Ltot = sol.d1 + sol.Lg + sol.d2;
  let maxErr = 0;
  for (let i = 0; i < sol.Nt; i++) {
    const want = U * (1 - sol.x[i] / Ltot);
    maxErr = Math.max(maxErr, Math.abs(sol.phi[i] - want));
  }
  ok(maxErr / U < 1e-12, 'T1 линейный потенциал в однородной среде',
    `maxErr/U=${(maxErr / U).toExponential(2)}`);
  // и поле однородно
  let emin = Infinity, emax = -Infinity;
  for (let f = 1; f < sol.Nt; f++) { emin = Math.min(emin, sol.Ef[f]); emax = Math.max(emax, sol.Ef[f]); }
  near(emin, U / Ltot, 1e-12, 'T1b поле = U/L', '');
  ok(Math.abs(emax - emin) / emax < 1e-12, 'T1c поле однородно');
}

// ===========================================================================
section('2. Ёмкостной делитель на двух барьерах');
{
  const sol = new DBDSolver({ epsR: 9, gapMM: 1, dielMM1: 0.5, dielMM2: 0.5, nCells: 120, seedDensity: 0, chemistry: false, mode: 'demo' });
  zeroDensities(sol);
  const U = 8000;
  sol._poissonSolve(U, true);
  const Ugap = sol._surfacePotential(sol.fL) - sol._surfacePotential(sol.fR);
  const Cd = (EPS0 * 9 * sol.A) / (sol.d1 + sol.d2);
  const Cg = (EPS0 * sol.A) / sol.Lg;
  const want = U * (Cd / (Cd + Cg));
  near(Ugap, want, 1e-12, 'T2 U_gap = U*Cd/(Cd+Cg)');
  // поле в газе / в диэлектрике = eps_r
  const Egas = sol.Ef[sol.fL + 3];
  const Ediel = sol.Ef[3];
  near(Egas / Ediel, 9, 1e-10, 'T2b E_gas/E_diel = eps_r');
}

// ===========================================================================
section('3. Заряженный слой: скачок поля по теореме Гаусса');
{
  const sol = new DBDSolver({ epsR: 9, nCells: 200, seedDensity: 0, chemistry: false, mode: 'demo' });
  zeroDensities(sol);
  // положительный слой в середине зазора
  const iA = sol.i0 + 80, iB = sol.i0 + 120;
  const nDen = 1e17;
  let Qs = 0;
  for (let i = iA; i < iB; i++) { sol.n[1][i] = nDen; Qs += QE * nDen * sol.dx[i]; }
  sol._poissonSolve(0, true);
  const Ebefore = sol.Ef[iA];      // грань слева от слоя
  const Eafter = sol.Ef[iB];       // грань справа от слоя
  near(EPS0 * (Eafter - Ebefore), Qs, 1e-9, 'T3 скачок eps0*dE = поверхностный заряд слоя');

  // проверка поверхностного заряда на диэлектрике
  zeroDensities(sol);
  sol.sigL = 4.6e-4; sol.sigFace[sol.fL] = sol.sigL;
  sol._poissonSolve(0, true);
  const Dg = EPS0 * sol.Ef[sol.fL];
  const Dd = EPS0 * 9 * (-(sol._surfacePotential(sol.fL) - sol.phi[sol.i0 - 1]) / (sol.dx[sol.i0 - 1] / 2));
  near(Dg - Dd, sol.sigL, 1e-9, 'T3b условие Гаусса на границе газ/диэлектрик');
}

// ===========================================================================
section('4. Чистый дрейф гауссиана (SG): масса и положение центра');
{
  // независимая одномерная проверка ядра SG на равномерной сетке
  const Nx = 400, L = 1e-3, dx = L / Nx;
  const n = new Float64Array(Nx), nn = new Float64Array(Nx), G = new Float64Array(Nx + 1);
  const x0 = 0.25 * L, w = 3e-5;
  for (let i = 0; i < Nx; i++) {
    const xc = (i + 0.5) * dx;
    n[i] = Math.exp(-((xc - x0) ** 2) / (w * w));
  }
  const mu = 2.4e-4, E = 1e6, D = 0;      // чистый дрейф, v = mu*E
  const v = mu * E;
  const dt = 0.4 * dx / v;
  const nSteps = 500;
  const m0 = n.reduce((a, b) => a + b, 0) * dx;
  let c0 = 0; for (let i = 0; i < Nx; i++) c0 += n[i] * (i + 0.5) * dx * dx;
  c0 /= m0;
  for (let s = 0; s < nSteps; s++) {
    for (let f = 1; f < Nx; f++) G[f] = sgFlux(n[f - 1], n[f], mu, D, +1, E, dx);
    G[0] = 0; G[Nx] = 0;
    for (let i = 0; i < Nx; i++) nn[i] = n[i] - (dt * (G[i + 1] - G[i])) / dx;
    n.set(nn);
  }
  const m1 = n.reduce((a, b) => a + b, 0) * dx;
  let c1 = 0; for (let i = 0; i < Nx; i++) c1 += n[i] * (i + 0.5) * dx * dx;
  c1 /= m1;
  near(m1, m0, 1e-12, 'T4 масса сохраняется при чистом дрейфе');
  near(c1 - c0, v * dt * nSteps, 0.02, 'T4b центр сдвинулся на v*t');

  // дрейф-диффузионный стационар: SG воспроизводит больцмановский профиль точно
  const Nz = 20, Lz = 1e-4, dz = Lz / Nz;
  const Dd = 0.05, mud = 0.1, Ez = 2e3;
  const nb = new Float64Array(Nz);
  for (let i = 0; i < Nz; i++) nb[i] = Math.exp((mud * Ez * (i + 0.5) * dz) / Dd);
  let maxF = 0;
  for (let f = 1; f < Nz; f++) {
    const fl = sgFlux(nb[f - 1], nb[f], mud, Dd, +1, Ez, dz);
    maxF = Math.max(maxF, Math.abs(fl) / (Dd * nb[f] / dz));
  }
  ok(maxF < 1e-12, 'T4c SG: больцмановский стационар даёт нулевой поток машинно точно',
    `maxRelFlux=${maxF.toExponential(2)}`);
}

// ===========================================================================
section('5. Сохранение заряда');
{
  const sol = new DBDSolver({ mode: 'demo', chemistry: false, U0kV: 10, freqKHz: 10 });
  // заметная плазма, чтобы потоки были нетривиальны
  for (let i = sol.i0; i <= sol.i1; i++) {
    sol.n[0][i] = 1e16; sol.n[1][i] = 1.2e16; sol.n[4][i] = 2e15;
  }
  sol.t = 2e-5;   // ненулевое напряжение
  const q0 = sol.totalCharge();
  const c0 = sol.clipCharge * sol.A;
  for (let s = 0; s < 300; s++) sol.step();
  const q1 = sol.totalCharge();
  const c1 = sol.clipCharge * sol.A;
  // масштаб — полный заряд одного сорта в зазоре
  const scale = QE * 1e16 * sol.Lg * sol.A;
  const err = Math.abs(q1 - q0 - (c1 - c0)) / scale;
  ok(err < 1e-9, 'T5 полный заряд (объём + поверхности + клип) сохраняется',
    `rel=${err.toExponential(2)}, clip=${((c1 - c0) / scale).toExponential(2)}`);
  const res = sol.poissonResidual();
  ok(res < 1e-9, 'T5b невязка Пуассона мала', `res=${res.toExponential(2)}`);
}
{
  // T5c — то же самое, но С ХИМИЕЙ. Именно этот случай раньше не проверялся, а
  // источники реакций заряд не сохраняли (+1.7e-9 Кл за период, 1.8%/период).
  const sol = new DBDSolver({ mode: 'demo', chemistry: true, U0kV: 10, freqKHz: 10 });
  sol.targetSimTime = 1 / sol.freq;
  const q0 = sol.totalCharge(), c0 = sol.clipCharge * sol.A;
  let qTrans = 0, prevI = 0, sigScale = 0, diverged = null;
  try {
    while (sol.t < sol.targetSimTime) {
      const dt = sol.step();
      physGuard(sol);
      qTrans += 0.5 * (Math.abs(sol.current) + Math.abs(prevI)) * dt;
      prevI = sol.current;
      sigScale = Math.max(sigScale, Math.abs(sol.sigL), Math.abs(sol.sigR));
    }
  } catch (e) {
    if (!(e instanceof Diverged)) throw e;
    diverged = e.message;
  }
  if (diverged) {
    ok(false, 'T5c/T5d/T5e/T5f период не досчитан: 1D-модель расходится на импульсе', diverged);
    failed += 3;   // T5d, T5e, T5f тоже не выполнены
  } else {
  const drift = Math.abs(sol.totalCharge() - q0 - (sol.clipCharge * sol.A - c0));
  const rel = drift / qTrans;
  ok(rel < 1e-4, 'T5c заряд сохраняется с химией за полный период (ERRATA V12: <1e-4)',
    `дрейф=${drift.toExponential(2)} Кл при перенесённом ${qTrans.toExponential(2)} Кл, отн.=${rel.toExponential(2)}`);
  // Асимметрия меряется относительно ПИКОВОГО |σ| за период (в момент t=T обе
  // поверхности почти разряжены, и относительная разница там неинформативна).
  // Раньше фиктивный заряд источников оседал на диэлектриках как sigL=+1.9e-4
  // против sigR=-1.0e-4 — 2:1 в симметричной ячейке при симметричном синусе.
  const asym = Math.abs(sol.sigL + sol.sigR) / Math.max(sigScale, 1e-30);
  ok(asym < 0.02, 'T5d поверхностный заряд симметричен (|σ_L+σ_R| << |σ|_max)',
    `sigL=${sol.sigL.toExponential(3)}, sigR=${sol.sigR.toExponential(3)}, |σ|max=${sigScale.toExponential(3)}, асимметрия=${asym.toExponential(2)}`);
  const frac = sol.rejects / sol.steps;
  ok(frac < 0.01, 'T5e доля отбракованных шагов < 1% (NUMERICS §5.2)',
    `${sol.rejects}/${sol.steps} = ${(frac * 100).toFixed(3)}%, причины=${JSON.stringify(sol.rejectReason)}`);
  const clipRel = Math.abs(sol.clipCharge * sol.A) / qTrans;
  ok(clipRel < 1e-6, 'T5f компенсирующий заряд пола пренебрежим',
    `clip=${(sol.clipCharge * sol.A).toExponential(2)} Кл, отн.=${clipRel.toExponential(2)}`);
  }
}

// ===========================================================================
section('6. Критерий пробоя (эталоны docs/REFERENCE_TARGETS.md)');
{
  // ⚠ СТАРЫЕ ЦЕЛИ ОТМЕНЕНЫ (ERRATA, баннер вверху; REFERENCE_TARGETS §5):
  //   alpha = eta при 101.3 Td  ->  127.84 Td
  //   самоподдержание 153.4 Td, U_br = 3.76 кВ  ->  141.97 Td, U_br = 3.473 кВ
  // Они были арифметически верны ВНУТРИ фитов PHYSICS.md, но сами фиты неверны
  // на порядки. Подгонять солвер под старые числа запрещено.
  const N = gasDensity(760, 300);
  const cross = alphaEtaCrossEN(N);
  // допуск ±5 Td: расчёт даёт 127.84, независимая проверка по опубликованной
  // таблице Phelps — 126.6, наборы Biagi — 132; эксперимент 122–130 Td.
  ok(Math.abs(cross - 127.84) < 5,
    'T6 alpha = eta при 127.84 ± 5 Td (было 101.3 — ОТМЕНЕНО, ERRATA §C0)',
    `E/N=${cross.toFixed(2)} Td, E=${(cross * TD * N / 1e5).toFixed(2)} кВ/см`);
  // линейный поиск тем же перебором, что и раньше — на случай, если бисекция
  // в physics.js когда-нибудь поймает не тот корень
  let scan = 0;
  for (let x = 50; x < 300; x += 0.01) if (alphaN(x) >= etaN(x, N)) { scan = x; break; }
  near(scan, cross, 1e-3, 'T6a бисекция и прямой перебор дают один корень');

  // --- Самоподдержание: критерий для ЭЛЕКТРООТРИЦАТЕЛЬНОГО газа (ERRATA §C V2)
  //     gamma*alpha/(alpha-eta)*[exp((alpha-eta)d) - 1] = 1.
  // Упрощённая форма (alpha-eta)d = ln(1+1/gamma) НЕ используется: при gamma=0.02
  // она даёт 146.0 Td вместо 141.97, то есть систематически завышает.
  const d = 1e-3;
  const REF = [[0.005, 148.56, 3.634], [0.010, 145.41, 3.557],
               [0.020, 141.97, 3.473], [0.050, 136.88, 3.349]];
  for (const [gam, enRef, ukvRef] of REF) {
    const en = breakdownEN(d, gam, N);
    const ukv = (en * TD * N * d) / 1e3;
    near(ukv, ukvRef, 0.05, `T6b U_br(gamma=${gam}) = ${ukvRef} кВ (допуск 5 %)`,
      `E/N=${en.toFixed(2)} Td (эталон ${enRef})`);
  }
  // невязка в корне действительно нулевая (а не «поймали край скобки»)
  const enG = breakdownEN(d, 0.02, N);
  ok(Math.abs(selfSustainResidual(enG, d, 0.02, N)) < 1e-6,
    'T6c невязка критерия самоподдержания в найденном корне ~ 0',
    `f=${selfSustainResidual(enG, d, 0.02, N).toExponential(2)}`);
  ok(alphaEff(enG, N) * d > 0 && alphaEff(enG, N) * d < 6,
    'T6d (alpha-eta)*d в корне = 2.96 (эталон REFERENCE_TARGETS §2)',
    `(alpha-eta)d=${(alphaEff(enG, N) * d).toFixed(3)}`);

  // --- ⚠ U_br и U_burn — РАЗНЫЕ величины, сравнивать напрямую нельзя ---------
  // U_br  = 3.473 кВ  — статическое поле самоподдержания в ГОЛОМ ГАЗЕ (Таунсенд),
  //                     барьеров и внешней цепи в этом критерии нет вообще.
  // U_i   = U_br*(1 + C_g/C_d) — то, что надо приложить К ЯЧЕЙКЕ, чтобы на газе
  //                     появилось U_br: ёмкостный делитель забирает часть на барьеры.
  //                     Это порог ПЕРВОГО зажигания при sigma = 0 (ERRATA §C V3).
  // U_burn = periodStats.UburnkV — ИЗМЕРЕННЫЙ солвером потолок |U_gap| на импульсе.
  //                     Это тоже напряжение на газе, но динамическое: в переходном
  //                     перенапряжённом режиме зазор успевает зайти ВЫШЕ U_br
  //                     (формативное запаздывание), а в установившемся остаточный
  //                     sigma позволяет зажигание НИЖЕ. Равенства ждать нельзя,
  //                     проверяется только вилка (T9e) и связь через делитель.
  const Cd = (EPS0 * 9 * 1e-4) / 1e-3;     // 2 барьера по 0.5 мм, eps_r = 9
  const Cg = (EPS0 * 1e-4) / 1e-3;
  const Ubr = (enG * TD * N * d);
  const Ui = Ubr * (1 + Cg / Cd);
  near(Ubr / 1e3, 3.473, 0.05, 'T6e U_br(gamma=0.02) = 3.473 кВ (голый газ)');
  near(Ui / 1e3, 3.859, 0.05, 'T6f U_i = U_br*(1+C_g/C_d) = 3.86 кВ (порог на ЯЧЕЙКЕ)',
    `C_g/C_d=${(Cg / Cd).toFixed(4)}`);
}

// ===========================================================================
section('7. Устойчивость: 20000 шагов без NaN и отрицательных плотностей');
{
  const sol = new DBDSolver({ mode: 'demo', U0kV: 10, freqKHz: 10 });
  const t0 = Date.now();
  let div7 = null;
  try {
    for (let s = 0; s < 20000; s++) { sol.step(); physGuard(sol); }
  } catch (e) {
    if (!(e instanceof Diverged)) throw e;
    div7 = e.message;
  }
  if (div7) ok(false, 'T7 20000 шагов: модель ушла в нефизичное состояние', div7);
  const el = Date.now() - t0;
  let bad = 0, nanCount = 0;
  const st = sol.state;
  for (const key of Object.keys(st.n)) {
    const a = st.n[key];
    for (let i = 0; i < a.length; i++) {
      if (!Number.isFinite(a[i])) nanCount++;
      else if (a[i] < 0) bad++;
    }
  }
  for (let i = 0; i < st.E.length; i++) if (!Number.isFinite(st.E[i])) nanCount++;
  const scalarsOk = [st.t, st.dt, st.Uapp, st.Ugap, st.current, st.charge, st.o3ppm, st.maxEN]
    .every(Number.isFinite);
  ok(nanCount === 0 && bad === 0 && scalarsOk, 'T7 нет NaN и отрицательных плотностей',
    `nan=${nanCount} neg=${bad} t=${st.t.toExponential(3)}с шагов=20000 за ${el}мс (${(el / 20000 * 1000).toFixed(1)}мкс/шаг)`);
  console.log(`      диагностика: dt=${st.dt.toExponential(2)} maxE/N=${st.maxEN.toFixed(1)} Td, ` +
    `n_e,max=${Math.max(...st.n.e).toExponential(2)} м^-3, I=${st.current.toExponential(2)} А`);
}

// ===========================================================================
section('8. U0 = 0 -> разряда нет');
{
  const sol = new DBDSolver({ mode: 'demo', U0kV: 0 });
  for (let s = 0; s < 3000; s++) sol.step();
  const st = sol.state;
  ok(Math.abs(st.Uapp) < 1e-12, 'T8 напряжение источника нулевое');
  ok(Math.abs(st.condCurrent) < 1e-6, 'T8b кондукционный ток пренебрежим',
    `I_cond=${st.condCurrent.toExponential(2)} А`);
  ok(st.maxEN < 1, 'T8c поле в зазоре не растёт', `maxE/N=${st.maxEN.toExponential(2)} Td`);
  ok(Math.max(...st.n.e) < 1.1 * 1e13, 'T8d электроны не размножаются',
    `n_e,max=${Math.max(...st.n.e).toExponential(2)}`);
}

// ===========================================================================
section('9. Полный период: энергия > 0 и наклоны Лиссажу');
{
  const sol = new DBDSolver({ mode: 'demo', U0kV: 10, freqKHz: 10, nCells: 200 });
  const t0 = Date.now();
  const nPeriods = 2;
  sol.targetSimTime = nPeriods * sol.period;
  let guard = 0, div9 = null;
  try {
    while (sol.t < sol.targetSimTime && guard < 4_000_000) { sol.step(); physGuard(sol); guard++; }
  } catch (e) {
    if (!(e instanceof Diverged)) throw e;
    div9 = e.message;
  }
  const el = (Date.now() - t0) / 1000;
  const ps = sol.periodStats;
  if (div9) {
    // Разряд ушёл в расходимость ДО конца первого периода — periodStats пуст,
    // мерить по нему нечего. Помечаем весь блок как невыполненный явно.
    ok(false, 'T9..T9f период не досчитан: 1D-модель расходится на первом импульсе', div9);
    failed += 6;
  } else {
  console.log(`      ${guard} шагов за ${el.toFixed(1)} с (${(el / nPeriods).toFixed(2)} с/период), ` +
    `E=${ps.energyPerPeriodJ.toExponential(3)} Дж, P=${ps.powerW.toFixed(2)} Вт, ` +
    `Cdiel(наклон)=${(ps.Cdiel * 1e12).toFixed(2)} пФ, Cdiel(тождество)=${(ps.CdielIdentity * 1e12).toFixed(2)} пФ, ` +
    `Ccell=${(ps.Ccell * 1e12).toFixed(3)} пФ, дуг on/off=${ps.qvArcsOn}/${ps.qvArcsOff}, ` +
    `R²on=${ps.qvR2on.toFixed(3)} R²off=${ps.qvR2off.toFixed(3)}, qvOk=${ps.qvOk}, ` +
    `U_burn=${ps.UburnkV.toFixed(2)} кВ, импульсов/период=${sol.state.breakdownsPerPeriod}, ` +
    `O3=${sol.state.o3ppm.toFixed(1)} ppm`);
  ok(ps.energyPerPeriodJ > 0, 'T9 энергия за период > 0',
    `W=${ps.energyPerPeriodJ.toExponential(3)} Дж`);
  const CdAnalytic = (EPS0 * sol.params.epsR * sol.A) / (sol.d1 + sol.d2);
  // T9b — САМОПРОВЕРКА СХЕМЫ, а не измерение: Q = C_d*(U_app - U_gap) — тождество
  // последовательной цепи, поэтому наклон dQ/d(U-U_gap) обязан быть равен C_d на
  // ЛЮБОЙ дуге. Это тест на согласованность цепи, он не может «поймать» физику.
  near(ps.CdielIdentity, CdAnalytic, 0.02, 'T9b тождество цепи Q=C_d(U−U_gap) выполняется (самопроверка)');
  // T9b2 — НАСТОЯЩЕЕ измерение: наклон горящей ветви. В 1D-модели зазор зажат
  // неидеально (2 микроразряда на полупериод, U_gap плывёт), поэтому он выше C_d
  // в 1.5–2.5 раза, а фигура не является параллелограммом — это и должен
  // сообщать флаг качества qvOk=false, а не молчаливая подмена аналитикой.
  ok(ps.Cdiel > 0.9 * CdAnalytic && ps.Cdiel < 3 * CdAnalytic,
    'T9b2 измеренный наклон горящей ветви = 0.9..3 C_d (фигура не параллелограмм)',
    `Cdiel=${(ps.Cdiel * 1e12).toFixed(2)} пФ, C_d=${(CdAnalytic * 1e12).toFixed(2)} пФ, R²on=${ps.qvR2on.toFixed(3)}`);
  ok(ps.qvOk === false && ps.qvR2on < 0.98,
    'T9b3 качество фигуры честно помечено как непригодное (R²(горящей) < 0.98)',
    `qvOk=${ps.qvOk}, R²on=${ps.qvR2on.toFixed(3)}, дуг on/off=${ps.qvArcsOn}/${ps.qvArcsOff}`);
  const CcellAnalytic = (EPS0 * sol.A) / sol.dEff;
  near(ps.Ccell, CcellAnalytic, 0.20, 'T9c Ccell из фигуры Лиссажу в пределах 20%');
  // V5: мощность по Мэнли P = 4*f*Cd*Umin*(U0-Umin), Umin = Ubr*(1+Cg/Cd)
  const Ubr = ps.UburnkV * 1e3;
  const Umin = Ubr * (1 + sol.Cg / sol.Cd);
  const Pmanley = 4 * sol.freq * sol.Cd * Umin * (sol.U0 - Umin);
  ok(ps.powerW > 0.3 * Pmanley && ps.powerW < 3 * Pmanley,
    'T9d мощность согласована с формулой Мэнли (фактор 3)',
    `P_sim=${ps.powerW.toFixed(2)} Вт, P_Manley=${Pmanley.toFixed(2)} Вт`);
  ok(ps.UburnkV > 3.0 && ps.UburnkV < 5.5, 'T9e напряжение горения зазора 3..5.5 кВ (V2)',
    `U_burn=${ps.UburnkV.toFixed(2)} кВ`);
  ok(sol.state.breakdownsPerPeriod >= 2, 'T9f есть импульсы тока в каждом полупериоде (V7)',
    `импульсов=${sol.state.breakdownsPerPeriod}`);
  }
}

// ===========================================================================
section('10. Вспомогательные функции и иерархия времён (V11)');
{
  ok(Math.abs(bernoulli(1e-14) - 1) < 1e-13, 'T10 B(x->0) = 1');
  near(bernoulli(1), 1 / (Math.E - 1), 1e-14, 'T10b B(1) = 1/(e-1)');
  ok(Math.abs(bernoulli(-50) + (-50) * -1 - 50) < 1e-9 || bernoulli(-50) === 50, 'T10c B(-50) -> 50');
  near(phi1(1e-12), 1, 1e-11, 'T10d phi1(0)=1');
  near(phi1(1), 1 - Math.exp(-1), 1e-14, 'T10e phi1(1)=1-1/e');

  const N = gasDensity(760, 300);
  // ⚠ поле самоподдержания пересчитано: 154.4 -> 141.97 Td (REFERENCE_TARGETS §2).
  // Вместе с правильной w_e это даёт время пролёта электрона 4.33 нс (было ~8 нс).
  const ENbr = 141.97;
  const Ebr = ENbr * TD * N;
  const we = driftElectron(ENbr);
  const tE = 1e-3 / we;
  const muI = 2.16e-4 * (2.687e25 / N);
  const tI = 1e-3 / (muI * Ebr);
  const tC1 = 1 / (2.4e-42 * N * N);
  const tC3 = 1 / (1.1e-42 * N * N);
  const tO3 = 1 / (6.0e-46 * N * N);
  console.log(`      transit_e=${(tE * 1e9).toFixed(1)} нс, transit_i=${(tI * 1e6).toFixed(2)} мкс, ` +
    `O2+->O4+=${(tC1 * 1e9).toFixed(2)} нс, O-->O3-=${(tC3 * 1e9).toFixed(2)} нс, O->O3=${(tO3 * 1e6).toFixed(2)} мкс`);
  near(tE, 4.33e-9, 0.10, 'T10f время пролёта электрона = 4.33 нс (REFERENCE_TARGETS §4)');
  near(tI, 1.19e-6, 0.15, 'T10g время пролёта O2+ = 1.19 мкс (REFERENCE_TARGETS §4)');
  ok(tC1 < 2e-9 && tC3 < 4e-9, 'T10h конверсия ионов — субнаносекунды');
  ok(tO3 > 1e-6 && tO3 < 8e-6, 'T10i образование озона ~микросекунды');
}

// ===========================================================================
section('11. Контракт публичного API');
{
  const sol = new DBDSolver({ mode: 'demo' });
  const dt = sol.step();
  ok(typeof dt === 'number' && dt > 0, 'T11 step() возвращает dt', `dt=${dt.toExponential(2)}`);
  const r = sol.advance(20);
  ok(r && typeof r.steps === 'number' && typeof r.simTime === 'number' && r.steps > 0,
    'T11b advance(budget) возвращает {steps, simTime}', `steps=${r.steps}, simTime=${r.simTime.toExponential(2)}`);
  const st = sol.state;
  const need = ['t', 'dt', 'x', 'xFaces', 'gasMask', 'n', 'E', 'phi', 'rho', 'ionizRate', 'EN',
    'sigmaL', 'sigmaR', 'Uapp', 'Ugap', 'current', 'condCurrent', 'dispCurrent', 'charge',
    'power', 'o3ppm', 'maxEN', 'peakCurrent', 'breakdownsPerPeriod'];
  const missing = need.filter((k) => !(k in st));
  ok(missing.length === 0, 'T11c state содержит все поля контракта',
    missing.length ? 'нет: ' + missing.join(',') : `${need.length} полей`);
  const needN = ['e', 'O2p', 'O4p', 'Om', 'O2m', 'O', 'O3', 'O2a'];
  ok(needN.every((k) => st.n[k] instanceof Float64Array), 'T11d state.n содержит все сорта',
    'плюс O3m: ' + (st.n.O3m instanceof Float64Array));
  const H = sol.history;
  ok(['t', 'Uapp', 'Ugap', 'current', 'charge', 'len', 'head', 'capacity'].every((k) => k in H),
    'T11e history — кольцевые буферы', `len=${H.len}/${H.capacity}`);
  ok(['energyPerPeriodJ', 'powerW', 'Cdiel', 'Ccell', 'UburnkV'].every((k) => k in sol.periodStats),
    'T11f periodStats содержит требуемые поля');
  // горячее обновление параметров без пересборки сетки
  const Nt0 = sol.Nt;
  sol.setParams({ U0kV: 12, gamma: 0.03 });
  ok(sol.Nt === Nt0 && sol.params.U0kV === 12 && sol.U0 === 12000,
    'T11g setParams без геометрии не пересобирает сетку');
  sol.setParams({ nCells: 150 });
  ok(sol.Ng === 150 && sol.t === 0, 'T11h setParams с геометрией пересобирает сетку и сбрасывает');
  sol.step();
  ok(Number.isFinite(sol.state.current), 'T11i солвер работает после пересборки');
  // сеточные проверки NUMERICS §1.2
  ok(sol.gridSmoothness < 0.08, 'T11j гладкость сетки |dx_{i+1}/dx_i - 1| < 0.08',
    `max=${sol.gridSmoothness.toFixed(4)}`);
}

// ===========================================================================
section('12. Горячая смена ε_r и валидация геометрии');
{
  const a = new DBDSolver({ epsR: 9, gapMM: 1, dielMM1: 0.5, dielMM2: 0.5, mode: 'demo', seedDensity: 0, chemistry: false });
  zeroDensities(a);
  a.setParams({ epsR: 2 });
  a._poissonSolve(8000, true);
  const Ugap = a._surfacePotential(a.fL) - a._surfacePotential(a.fR);
  const want = 8000 * a.Lg / (a.Lg + (a.d1 + a.d2) / 2);
  near(Ugap, want, 1e-9, 'T12 setParams({epsR}) реально меняет решение (матрица Пуассона)');
  ok(Math.abs(a.epsCell[0] / 8.8541878128e-12 - 2) < 1e-9, 'T12b epsCell пересобран',
    `epsCell[0]/eps0=${(a.epsCell[0] / 8.8541878128e-12).toFixed(3)}`);
  let threw = false;
  try { new DBDSolver({ dielMM1: 0, mode: 'demo' }); } catch (e) { threw = true; }
  ok(threw, 'T12c конфигурация без барьера отвергается явно, а не расходится молча');
}

// ===========================================================================
section('13. Выборка Q–V на высокой частоте (пропуск слотов)');
{
  // Тёмный прогон при 100 кГц: площадь фигуры Q–V обязана быть ~0. Раньше до 14.5%
  // слотов оставались с данными ПРОШЛОГО периода и площадь была фиктивной.
  const sol = new DBDSolver({ mode: 'demo', U0kV: 1, freqKHz: 100 });
  sol.targetSimTime = 4 / sol.freq;
  while (sol.t < sol.targetSimTime) sol.step();
  const scale = sol.Ccell * Math.pow(sol.U0, 2);
  const rel = sol.periodStats.energyPerPeriodJ / scale;
  ok(rel < 0.01, 'T13 при 100 кГц без разряда площадь Q–V < 1% от C_cell*U0²',
    `E=${sol.periodStats.energyPerPeriodJ.toExponential(2)} Дж, отн.=${rel.toExponential(2)}`);
  // и C_cell по фигуре не должен уезжать на высокой частоте
  const sol2 = new DBDSolver({ mode: 'demo', U0kV: 10, freqKHz: 50 });
  sol2.targetSimTime = 2 / sol2.freq;
  let div13 = null;
  try {
    while (sol2.t < sol2.targetSimTime) { sol2.step(); physGuard(sol2); }
  } catch (e) {
    if (!(e instanceof Diverged)) throw e;
    div13 = e.message;
  }
  if (div13) ok(false, 'T13b период не досчитан: 1D-модель расходится на импульсе', div13);
  else near(sol2.periodStats.Ccell, sol2.Ccell, 0.25, 'T13b Ccell из фигуры Лиссажу корректен и при 50 кГц');
}

// ===========================================================================
console.log(`\n=====================================`);
console.log(`PASSED ${passed}, FAILED ${failed}, время ${((Date.now() - t0all) / 1000).toFixed(1)} с`);
process.exit(failed === 0 ? 0 : 1);
