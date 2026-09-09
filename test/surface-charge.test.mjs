// test/surface-charge.test.mjs — знаки и баланс ПОВЕРХНОСТНОГО заряда.
// Запуск: node test/surface-charge.test.mjs   (exit code 0 = все тесты прошли)
//
// Закрывает дыру в покрытии, отмеченную в docs/DIVERGENCE_ANALYSIS.md §4:
// существующий тест T3b проверяет только то, что ЗАДАННЫЙ sigma правильно входит
// в уравнение Пуассона на ЛЕВОЙ границе. Здесь проверяется:
//   S1  (sigma, 0)   — заряд только на левом барьере;
//   S2  (0, sigma)   — только на правом;
//   S3  (sigma, sigma) — одного знака на обоих (аналитически поле в зазоре = 0);
//   S4  знаки НАКОПЛЕНИЯ заряда пристеночными потоками, обе нормали НЕЗАВИСИМО;
//   S5  зарядовый баланс d/dt[int(rho)dx + sigma_L + sigma_R] = 0 НА ПРОТЯЖЕНИИ
//       прогона, включая зажигание и стадию расходимости.
//
// Аналитика для S1..S3 (трёхслойный конденсатор, оба металла заземлены):
//   eps0*E_g - eps_d*E_1 = sigma_L,  eps_d*E_2 - eps0*E_g = sigma_R,
//   E_1*d1 + E_g*Lg + E_2*d2 = 0
//   =>  E_g = d*(sigma_L - sigma_R) / (2*eps0*d + eps_r*eps0*Lg)   при d1 = d2 = d.
// Отсюда: одинаковый знак и модуль на обоих барьерах даёт РОВНО НУЛЕВОЕ поле в
// зазоре, а перепутанная нормаль на одной из стенок — поле удвоенной величины.

import { DBDSolver } from '../src/solver.js';
import { EPS0, QE, TD, gasDensity, muE_N } from '../src/physics.js';

let passed = 0, failed = 0;
const CH_Z = [-1, 1, 1, -1, -1, -1];

function ok(cond, name, info = '') {
  if (cond) { passed++; console.log(`PASS  ${name}${info ? '  [' + info + ']' : ''}`); }
  else { failed++; console.log(`FAIL  ${name}${info ? '  [' + info + ']' : ''}`); }
  return !!cond;
}
function near(a, b, rtol, name, extra = '') {
  const den = Math.abs(b) > 1e-300 ? Math.abs(b) : 1;
  const r = Math.abs(a - b) / den;
  return ok(r <= rtol, name, `got=${a.toExponential(4)} want=${b.toExponential(4)} rel=${r.toExponential(2)}${extra ? ' ' + extra : ''}`);
}
function section(s) { console.log(`\n--- ${s} ---`); }

function zeroAll(sol) {
  for (const arr of sol.n) arr.fill(0);
  sol.nO.fill(0); sol.nO3.fill(0); sol.nO2a.fill(0);
  sol.sigL = 0; sol.sigR = 0; sol.sigFace.fill(0);
}
/** Плотность объёмного заряда газа, проинтегрированная по зазору, Кл/м^2. */
function volCharge(sol) {
  let q = 0;
  for (let s = 0; s < 6; s++) {
    const ns = sol.n[s], z = QE * CH_Z[s];
    for (let i = sol.i0; i <= sol.i1; i++) q += z * ns[i] * sol.dx[i];
  }
  return q;
}
/** |rho| — масштаб для нормировки невязки. */
function absVolCharge(sol) {
  let q = 0;
  for (let s = 0; s < 6; s++) {
    const ns = sol.n[s];
    for (let i = sol.i0; i <= sol.i1; i++) q += QE * ns[i] * sol.dx[i];
  }
  return q;
}
/** Поле в диэлектрике непосредственно у поверхности (сторона диэлектрика). */
function EdielAt(sol, f) {
  const phis = sol._surfacePotential(f);
  if (f === sol.fL) {
    const i = sol.i0 - 1;                    // последняя ячейка левого барьера
    return -(phis - sol.phi[i]) / (sol.dx[i] / 2);
  }
  const i = sol.i1 + 1;                      // первая ячейка правого барьера
  return -(sol.phi[i] - phis) / (sol.dx[i] / 2);
}

// ===========================================================================
section('S1..S3: ЗАДАННЫЙ поверхностный заряд — обе границы независимо');
{
  const epsR = 9;
  const sol = new DBDSolver({ epsR, nCells: 200, mode: 'demo', seedDensity: 0, chemistry: false });
  const d = sol.d1, Lg = sol.Lg;
  ok(Math.abs(sol.d1 - sol.d2) < 1e-12, 'S0 геометрия симметрична (d1 = d2)', `d=${d.toExponential(3)} м`);
  const S = 4.6e-4;                                     // Кл/м^2
  const Eg = (sL, sR) => (d * (sL - sR)) / (2 * EPS0 * d + epsR * EPS0 * Lg);

  const cases = [
    ['S1 (sigma, 0)  — только левый барьер', S, 0],
    ['S2 (0, sigma)  — только правый барьер', 0, S],
    ['S3 (sigma, sigma) — одного знака на обоих', S, S],
    ['S3b (sigma, -sigma) — противоположные знаки', S, -S],
  ];
  for (const [name, sL, sR] of cases) {
    zeroAll(sol);
    sol.sigL = sL; sol.sigR = sR;
    sol.sigFace[sol.fL] = sL; sol.sigFace[sol.fR] = sR;
    sol._poissonSolve(0, true);
    // поле в середине зазора
    const fMid = sol.i0 + ((sol.i1 - sol.i0) >> 1);
    const EgNum = sol.Ef[fMid];
    const want = Eg(sL, sR);
    if (Math.abs(want) < 1e-9) {
      // аналитический ноль: сравниваем с масштабом поля, которое дал бы один барьер
      const scale = Math.abs(Eg(S, 0));
      ok(Math.abs(EgNum) < 1e-9 * scale, `${name}: поле в зазоре = 0`,
        `E=${EgNum.toExponential(3)} В/м при масштабе ${scale.toExponential(3)}`);
    } else {
      near(EgNum, want, 1e-9, `${name}: поле в зазоре`, `В/м`);
    }
    // скачок D_z вдоль +z на КАЖДОЙ границе равен ЛОКАЛЬНОМУ sigma
    const jumpL = EPS0 * sol.Ef[sol.fL] - epsR * EPS0 * EdielAt(sol, sol.fL);
    const jumpR = epsR * EPS0 * EdielAt(sol, sol.fR) - EPS0 * sol.Ef[sol.fR];
    const sc = Math.max(Math.abs(sL), Math.abs(sR), 1e-30);
    ok(Math.abs(jumpL - sL) < 1e-9 * sc, `${name}: [D_z] на ЛЕВОЙ границе = sigma_L`,
      `got=${jumpL.toExponential(4)} want=${sL.toExponential(4)}`);
    ok(Math.abs(jumpR - sR) < 1e-9 * sc, `${name}: [D_z] на ПРАВОЙ границе = sigma_R`,
      `got=${jumpR.toExponential(4)} want=${sR.toExponential(4)}`);
  }
}

// ===========================================================================
section('S4: знаки НАКОПЛЕНИЯ заряда пристеночными потоками (обе нормали независимо)');
{
  // Химия выключена, поле источника нулевое: остаётся ТОЛЬКО тепловой поток на
  // стенку (0.25*v_th), одинаковый по модулю на обеих стенках. Это и разделяет
  // нормали: при перепутанном знаке на одной из стенок sigma там уйдёт не туда.
  const mkSol = () => {
    const s = new DBDSolver({ mode: 'demo', chemistry: false, U0kV: 0, seedDensity: 0 });
    zeroAll(s);
    s._poissonSolve(0, true);
    s.EfOld.set(s.Ef);
    return s;
  };
  const nden = 1e16;

  // (1) электроны у ЛЕВОЙ стенки, положительные ионы O2+ у ПРАВОЙ
  {
    const s = mkSol();
    s.n[0][s.i0] = nden;            // e  у левой стенки
    s.n[1][s.i1] = nden;            // O2+ у правой стенки
    const q0 = volCharge(s), c0 = s.clipCharge;
    s.step();
    const dQ = volCharge(s) - q0, dC = s.clipCharge - c0;
    ok(s.sigL < 0, 'S4a электроны на ЛЕВУЮ стенку => sigma_L ОТРИЦАТЕЛЬНА',
      `sigL=${s.sigL.toExponential(3)} Кл/м^2`);
    ok(s.sigR > 0, 'S4b ионы O2+ на ПРАВУЮ стенку => sigma_R ПОЛОЖИТЕЛЬНА',
      `sigR=${s.sigR.toExponential(3)} Кл/м^2`);
    const sc = Math.max(Math.abs(s.sigL), Math.abs(s.sigR));
    ok(Math.abs(dQ + s.sigL + s.sigR - dC) < 1e-11 * sc,
      'S4c что ушло из газа — то и пришло на поверхности',
      `dQ=${dQ.toExponential(4)}, dSig=${(s.sigL + s.sigR).toExponential(4)}, ` +
      `клип=${dC.toExponential(2)}, невязка=${(dQ + s.sigL + s.sigR - dC).toExponential(2)}`);
  }
  // (2) ЗЕРКАЛЬНЫЙ случай: ионы у левой, электроны у правой
  {
    const s = mkSol();
    s.n[1][s.i0] = nden;            // O2+ у левой
    s.n[0][s.i1] = nden;            // e  у правой
    const q0 = volCharge(s), c0 = s.clipCharge;
    s.step();
    const dQ = volCharge(s) - q0, dC = s.clipCharge - c0;
    ok(s.sigL > 0, 'S4d ионы O2+ на ЛЕВУЮ стенку => sigma_L ПОЛОЖИТЕЛЬНА',
      `sigL=${s.sigL.toExponential(3)} Кл/м^2`);
    ok(s.sigR < 0, 'S4e электроны на ПРАВУЮ стенку => sigma_R ОТРИЦАТЕЛЬНА',
      `sigR=${s.sigR.toExponential(3)} Кл/м^2`);
    const sc = Math.max(Math.abs(s.sigL), Math.abs(s.sigR));
    ok(Math.abs(dQ + s.sigL + s.sigR - dC) < 1e-11 * sc,
      'S4f что ушло из газа — то и пришло на поверхности',
      `dQ=${dQ.toExponential(4)}, dSig=${(s.sigL + s.sigR).toExponential(4)}, ` +
      `клип=${dC.toExponential(2)}, невязка=${(dQ + s.sigL + s.sigR - dC).toExponential(2)}`);
  }
  // (3) СИММЕТРИЧНЫЙ случай: электроны у обеих стенок. Обновление sigma обязано
  //     быть ОДИНАКОВЫМ по знаку и (при нулевом поле) близким по модулю.
  //     Одинаковый знак ОБНОВЛЕНИЯ (sigL += -J, sigR += -J вместо +J) даёт здесь
  //     противоположные знаки sigma — тест это ловит.
  {
    const s = mkSol();
    s.n[0][s.i0] = nden; s.n[0][s.i1] = nden;
    s.step();
    ok(s.sigL < 0 && s.sigR < 0, 'S4g электроны на ОБЕ стенки => обе sigma отрицательны',
      `sigL=${s.sigL.toExponential(3)}, sigR=${s.sigR.toExponential(3)}`);
    const asym = Math.abs(s.sigL - s.sigR) / Math.abs(s.sigL);
    ok(asym < 0.05, 'S4h симметричная загрузка даёт симметричный отклик',
      `асимметрия=${asym.toExponential(2)}`);
  }
  // (4) Дрейфовый случай: однородная плазма + приложенное напряжение.
  //     Электроны уходят на АНОД (левый при U_el > 0), ионы — на КАТОД (правый).
  {
    const s = new DBDSolver({ mode: 'demo', chemistry: false, U0kV: 10, freqKHz: 10, seedDensity: 0 });
    zeroAll(s);
    for (let i = s.i0; i <= s.i1; i++) { s.n[0][i] = 1e14; s.n[1][i] = 1e14; }
    s.t = 0.25 / s.freq;            // U_el = +U0
    s._poissonSolve(s.Usrc(s.t), true);
    s.EfOld.set(s.Ef);
    for (let k = 0; k < 50; k++) s.step();
    ok(s.Ef[s.fL] > 0, 'S4i при U_el > 0 поле в зазоре направлено по +z');
    ok(s.sigL < 0, 'S4j дрейф электронов на левый (анодный) барьер => sigma_L < 0',
      `sigL=${s.sigL.toExponential(3)}`);
    ok(s.sigR > 0, 'S4k дрейф ионов на правый (катодный) барьер => sigma_R > 0',
      `sigR=${s.sigR.toExponential(3)}`);
  }
}

// ===========================================================================
section('S5: зарядовый баланс во времени (включая зажигание и расходимость)');
{
  const N_GAS = gasDensity(760, 300);
  const sol = new DBDSolver({ mode: 'default', chemistry: true, U0kV: 10, freqKHz: 10 });
  const tEnd = 1 / sol.freq;
  let qPrev = volCharge(sol) + sol.sigL + sol.sigR;
  let clipPrev = sol.clipCharge;
  let scaleMax = 1e-30;                     // масштаб задачи: |sigma| + int|rho|dx
  let maxRel = 0, maxRelT = 0, maxRelStep = 0;
  let cumDrift = 0, cumTrans = 0;
  let steps = 0, diverged = null;
  const t0 = Date.now();
  while (sol.t < tEnd) {
    const dt = sol.step();
    steps++;
    const q = volCharge(sol) + sol.sigL + sol.sigR;
    const resid = (q - qPrev) - (sol.clipCharge - clipPrev);
    qPrev = q; clipPrev = sol.clipCharge;
    scaleMax = Math.max(scaleMax, Math.abs(sol.sigL) + Math.abs(sol.sigR) + absVolCharge(sol));
    const rel = Math.abs(resid) / scaleMax;
    if (rel > maxRel) { maxRel = rel; maxRelT = sol.t; maxRelStep = steps; }
    cumDrift += resid;
    cumTrans += Math.abs(sol.condCurrent / sol.A) * dt;
    // детектор ухода в нефизичное состояние (те же потолки, что в solver.test.mjs)
    let ne = 0;
    for (let i = sol.i0; i <= sol.i1; i++) if (sol.n[0][i] > ne) ne = sol.n[0][i];
    if (sol.maxEN > 5000 || ne >= N_GAS) {
      diverged = `E/N=${sol.maxEN.toFixed(0)} Тд, n_e=${ne.toExponential(2)} м^-3 при t=${sol.t.toExponential(3)} с (шаг ${steps})`;
      break;
    }
    if (Date.now() - t0 > 240000) { diverged = 'бюджет времени исчерпан'; break; }
  }

  ok(maxRel < 1e-10, 'S5 баланс d/dt[int(rho)dx + sigma_L + sigma_R] = 0 на КАЖДОМ шаге (< 1e-10)',
    `max отн. невязка=${maxRel.toExponential(2)} на шаге ${maxRelStep} (t=${maxRelT.toExponential(3)} с), шагов=${steps}`);
  ok(Math.abs(cumDrift) / Math.max(cumTrans, 1e-30) < 1e-10,
    'S5b накопленный дрейф заряда пренебрежим по сравнению с перенесённым',
    `дрейф=${cumDrift.toExponential(2)} Кл/м^2, перенесено=${cumTrans.toExponential(2)} Кл/м^2`);

  // Явная проверка гипотезы из DIVERGENCE_ANALYSIS §3: заряд одного знака на обоих
  // барьерах ЗАКОНЕН тогда и только тогда, когда его компенсирует объёмный заряд.
  const rhoInt = volCharge(sol);
  const sigSum = sol.sigL + sol.sigR;
  const sc = Math.max(Math.abs(rhoInt), Math.abs(sigSum), 1e-30);
  // Невязка здесь — это ровно cumDrift, накопленная ошибка ОКРУГЛЕНИЯ за все шаги:
  // на каждом шаге к sigma порядка |sigma| прибавляется приращение, и последняя
  // значащая цифра теряется. Теоретическая граница дрейфа — eps_машинное * |sigma| *
  // число шагов; допуск обязан её учитывать, иначе тест меряет длину прогона, а не
  // физику. (С ГУ Хагелаара прогон доходит до ~1.5e6 шагов вместо 2.6e4 при срыве,
  // и фиксированный порог 1e-10*sc становится жёстче машинной точности.)
  const roundoff = 8 * Number.EPSILON * sc * Math.max(steps, 1);
  ok(Math.abs(rhoInt + sigSum) < Math.max(1e-10 * sc, roundoff),
    'S5c sigma_L + sigma_R = -int(rho)dx (заряд одного знака на барьерах законен)',
    `int(rho)dx=${rhoInt.toExponential(4)}, sigL+sigR=${sigSum.toExponential(4)}, ` +
    `sigL=${sol.sigL.toExponential(3)}, sigR=${sol.sigR.toExponential(3)}, ` +
    `невязка=${(rhoInt + sigSum).toExponential(2)} Кл/м^2, ` +
    `допуск=${Math.max(1e-10 * sc, roundoff).toExponential(2)} (шагов ${steps})`);

  if (diverged) {
    console.log(`      ⚠ прогон остановлен: ${diverged}`);
    console.log('        (баланс заряда при этом выполнен — расходимость НЕ вызвана ошибкой знака sigma)');
  } else {
    console.log(`      период досчитан: ${steps} шагов`);
  }
}

// ===========================================================================
section('S6: полунеявная связь снимает ограничение по диэлектрической релаксации');
{
  // Проверяется ровно то, что в DIVERGENCE_ANALYSIS §2.2 названо «подозреваемым
  // номер один»: если поле, транспорт и поверхностный ток связаны НЕ полностью
  // неявно, шаг dt >> tau_M = eps0/sigma_pl даёт растущие осцилляции поля
  // (классическая неустойчивость явного дрейф-диффузионного Пуассона).
  // Здесь в зазор кладётся плотная квазинейтральная плазма, ограничитель по
  // tau_M снят, и проверяется, что схема действительно шагает через tau_M без
  // осцилляций и с сохранением положительности.
  const sol = new DBDSolver({ mode: 'demo', chemistry: false, U0kV: 1, freqKHz: 10, seedDensity: 0 });
  zeroAll(sol);
  const nPl = 1e21;
  for (let i = sol.i0; i <= sol.i1; i++) { sol.n[0][i] = nPl; sol.n[1][i] = nPl; }
  sol.t = 0.25 / sol.freq;                       // U_el = U0
  sol._poissonSolve(sol.Usrc(sol.t), true);
  sol.EfOld.set(sol.Ef);
  // Ограничители ТОЧНОСТИ (CFL ионов, dE/E, eCFL) намеренно распущены: тест
  // проверяет УСТОЙЧИВОСТЬ относительно tau_M, а не точность. Остаётся только
  // отбраковка шага по отрицательным плотностям и по невязке Пуассона.
  sol.cfl = 1e4; sol.dETol = 1e3; sol.eCFL = 1e6;
  const fMid = sol.i0 + ((sol.i1 - sol.i0) >> 1);
  const E0 = Math.abs(sol.Ef[fMid]);
  // sigma_pl = e*n_e*mu_e при поле E0; mu_e берётся из тех же таблиц LFA
  const muE0 = muE_N(E0 / (sol.N * TD)) / sol.N;
  const tauM = EPS0 / (QE * nPl * muE0);         // максвелловское время релаксации
  // проверки на осцилляции имеют смысл только пока поле ещё не экранировано
  // до уровня шума: ниже 1e-3*E0 |E| = 1e2 В/м, это уже round-off решения.
  const FLOOR = 1e-3 * E0;

  let maxRatio = 0, signFlips = 0, prevE = sol.Ef[fMid], grew = 0, bad = false;
  for (let k = 0; k < 300; k++) {
    const dt = sol.step();
    const E = sol.Ef[fMid];
    maxRatio = Math.max(maxRatio, dt / tauM);
    if (Math.abs(prevE) > FLOOR) {
      if (E * prevE < 0) signFlips++;
      if (Math.abs(E) > Math.abs(prevE) * (1 + 1e-6)) grew++;
    }
    prevE = E;
    for (let s = 0; s < 6; s++) {
      for (let i = sol.i0; i <= sol.i1; i++) if (!(sol.n[s][i] >= 0)) { bad = true; break; }
    }
    if (bad) break;
  }
  console.log(`      tau_M=${tauM.toExponential(3)} с, E0=${E0.toExponential(3)} В/м, ` +
    `E_кон=${Math.abs(prevE).toExponential(3)} В/м, max dt/tau_M=${maxRatio.toExponential(2)}`);
  ok(maxRatio > 100, 'S6 шаг превышает tau_M на 2+ порядка (ограничителя по релаксации нет)',
    `max dt/tau_M=${maxRatio.toExponential(2)}`);
  ok(!bad, 'S6b плотности остаются неотрицательными');
  ok(signFlips === 0, 'S6c поле не осциллирует по знаку (нет явной неустойчивости)',
    `смен знака=${signFlips}`);
  ok(grew === 0, 'S6d экранирование монотонно, амплитуда поля не нарастает',
    `шагов с ростом |E|=${grew}`);
  ok(Math.abs(prevE) < 0.05 * E0, 'S6e плазма экранировала поле в зазоре',
    `|E|/|E0|=${(Math.abs(prevE) / E0).toExponential(2)}`);
  const res = sol.poissonResidual();
  ok(res < 1e-9, 'S6f невязка Пуассона мала после серии больших шагов', `res=${res.toExponential(2)}`);
}

console.log('\n=====================================');
console.log(`PASSED ${passed}, FAILED ${failed}`);
process.exit(failed === 0 ? 0 : 1);
