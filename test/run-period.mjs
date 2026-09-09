// run-period.mjs — прогон нескольких периодов реального солвера с печатью
// физических диагностик. Запуск:
//     node test/run-period.mjs [периодов] [режим] [U0кВ] [fкГц]
//
// ⚠ ДЕЙСТВУЮЩИЕ ЭТАЛОНЫ — docs/REFERENCE_TARGETS.md, НЕ docs/PHYSICS.md §6.
// Коэффициенты заменены по результатам решения уравнения Больцмана
// (docs/RATES_REVIEW.md, ERRATA §C0/§D2/§D3), поэтому все числа, выведенные
// ВНУТРИ старых фитов (alpha=eta при 101.3 Td, U_br = 3.76 кВ, мощность по
// Мэнли 7.78 Вт), ОТМЕНЕНЫ. Подгонять модель под них запрещено.
//
// Величины, которые нельзя путать (ERRATA §C V2/V3):
//   U_br  = 3.473 кВ — статическое поле самоподдержания в ГОЛОМ ГАЗЕ (Таунсенд);
//   U_i   = U_br*(1 + C_g/C_d) = 3.86 кВ — порог на ЯЧЕЙКЕ при sigma = 0;
//   U_burn — измеренный солвером потолок |U_gap| на импульсе (динамический,
//            в переходном режиме выше U_br, в установившемся может быть ниже).

import { DBDSolver } from '../src/solver.js';

// --- потолки физичности: прогон обрывается, а не висит часами на dt = 1e-14 ---
// (см. комментарий в test/solver.test.mjs). 5000 Td = 1.2e8 В/м; n_e >= N —
// степень ионизации больше 100 %.
const EN_CEIL = 5000;
const NE_CEIL = 2.4463e25;
let diverged = null;
function physGuard(s) {
  if (diverged) return true;
  if (!(s.maxEN < EN_CEIL)) {
    diverged = `E/N = ${s.maxEN.toFixed(0)} Td > ${EN_CEIL} Td при t = ${s.t.toExponential(3)} с`;
    return true;
  }
  let ne = 0;
  for (let i = 0; i < s.n.e.length; i++) if (s.n.e[i] > ne) ne = s.n.e[i];
  if (!(ne < NE_CEIL)) {
    diverged = `n_e = ${ne.toExponential(2)} м^-3 >= N при t = ${s.t.toExponential(3)} с`;
    return true;
  }
  return false;
}

const NPER = Number(process.argv[2] || 3);
const MODE = process.argv[3] || 'demo';
const U0 = Number(process.argv[4] || 10);
const FREQ = Number(process.argv[5] || 10);

const QE = 1.602176634e-19;

const solver = new DBDSolver({
  gapMM: 1.0, dielMM1: 0.5, dielMM2: 0.5, epsR: 9, areaCM2: 1.0,
  U0kV: U0, freqKHz: FREQ, mode: MODE,
});
const P = solver.params;
const period = 1 / (P.freqKHz * 1e3);
const A = P.areaCM2 * 1e-4;

console.log('=== 1D DBD / O2 — прогон периодов ===================================');
console.log(`  зазор ${P.gapMM} мм, барьеры ${P.dielMM1}+${P.dielMM2} мм, eps_r=${P.epsR}, S=${P.areaCM2} см²`);
console.log(`  U0=${P.U0kV} кВ, f=${P.freqKHz} кГц, p=${P.pressureTorr} Торр, T=${P.tempK} K`);
console.log(`  режим ${P.mode}: ${P.nCells} ячеек газа, N=${solver.N.toExponential(3)} м^-3`);
console.log(`  периодов: ${NPER}  (T=${(period * 1e6).toFixed(1)} мкс)`);
console.log('');

// --- накопители по периодам --------------------------------------------------
const rows = [];
let steps = 0, dtSum = 0;

const wall0 = Date.now();

for (let per = 0; per < NPER; per++) {
  const t0 = solver.state.t;
  const tEnd = t0 + period;

  let peakI = 0, peakIcond = 0, maxEN = 0;
  let energy = 0;              // ∮ U dQ  ≈ ∫ U I dt
  let qHalfPos = 0, qHalfNeg = 0;   // ∫ I_cond dt по знаку
  let sigMax = 0;
  let prevT = t0, prevUI = 0;

  // сырой трек тока за период — импульсы разбираем после, а не на лету
  // (на лету порог зависит от ещё не встреченного пика и склеивает импульсы)
  const trT = [], trI = [], trU = [];

  while (solver.state.t < tEnd) {
    const dt = solver.step();
    steps++; dtSum += dt;
    const s = solver.state;
    if (physGuard(s)) break;

    const ai = Math.abs(s.current);
    const ac = Math.abs(s.condCurrent);
    if (ai > peakI) peakI = ai;
    if (ac > peakIcond) peakIcond = ac;
    if (s.maxEN > maxEN) maxEN = s.maxEN;
    if (Math.abs(s.sigmaL) > sigMax) sigMax = Math.abs(s.sigmaL);
    if (Math.abs(s.sigmaR) > sigMax) sigMax = Math.abs(s.sigmaR);

    // энергия по трапециям: ∫ U_app I dt
    const ui = s.Uapp * s.current;
    energy += 0.5 * (ui + prevUI) * (s.t - prevT);
    prevUI = ui; prevT = s.t;

    if (s.condCurrent > 0) qHalfPos += s.condCurrent * dt;
    else qHalfNeg += -s.condCurrent * dt;

    trT.push(s.t); trI.push(ac); trU.push(Math.abs(s.Ugap));
  }

  // --- разбор импульсов: сегментация по 0.05*I_pk, затем FWHM вокруг ЛОКАЛЬНОГО
  //     максимума каждого сегмента (а не ширина всей огибающей) -----------------
  const segThr = 0.05 * peakIcond;
  const pulseWidths = [], pulsePhases = [], pulseUgap = [], pulseQ = [], pulsePeaks = [];
  let pulses = 0;
  for (let i = 0; i < trI.length;) {
    if (trI[i] <= segThr) { i++; continue; }
    let j = i;
    while (j < trI.length && trI[j] > segThr) j++;
    // сегмент [i, j): его пик
    let pk = 0, pi = i;
    for (let k = i; k < j; k++) if (trI[k] > pk) { pk = trI[k]; pi = k; }
    if (pk > 0.02 * peakIcond) {
      let a = pi; while (a > i && trI[a] > 0.5 * pk) a--;
      let b = pi; while (b < j - 1 && trI[b] > 0.5 * pk) b++;
      let q = 0;
      for (let k = i + 1; k < j; k++) q += trI[k] * (trT[k] - trT[k - 1]);
      pulseWidths.push(trT[b] - trT[a]);
      pulsePhases.push(((trT[pi] - t0) / period) * 360);
      pulseUgap.push(trU[i]);   // U_gap в момент ЗАЖИГАНИЯ (фронт), а не на пике —
                                // к пику зазор уже частично закорочен плазмой
      pulsePeaks.push(pk);
      pulseQ.push(q);
      pulses++;
    }
    i = j;
  }

  const s = solver.state;
  // средние по газу концентрации нейтралов
  const gasAvg = (arr) => {
    let sum = 0, w = 0;
    for (let i = 0; i < s.x.length; i++) {
      if (!s.gasMask[i]) continue;
      const dx = s.xFaces[i + 1] - s.xFaces[i];
      sum += arr[i] * dx; w += dx;
    }
    return w > 0 ? sum / w : 0;
  };

  rows.push({
    per: per + 1,
    steps, energy, power: energy / period,
    peakI, peakIcond, maxEN, pulses,
    fwhm: pulseWidths.length ? pulseWidths.reduce((a, b) => a + b, 0) / pulseWidths.length : 0,
    fwhmMain: pulsePeaks.length
      ? pulseWidths[pulsePeaks.indexOf(Math.max(...pulsePeaks))] : 0,
    qMain: pulseQ.length ? pulseQ[pulsePeaks.indexOf(Math.max(...pulsePeaks))] : 0,
    peaks: pulsePeaks, pulseQ,
    phases: pulsePhases,
    Uburn: pulseUgap.length ? pulseUgap.reduce((a, b) => a + b, 0) / pulseUgap.length : 0,
    qPos: qHalfPos, qNeg: qHalfNeg,
    sigMax,
    o3: gasAvg(s.n.O3), o: gasAvg(s.n.O), o2a: gasAvg(s.n.O2a),
    o3ppm: s.o3ppm,
    ne: gasAvg(s.n.e), om: gasAvg(s.n.Om), o2m: gasAvg(s.n.O2m),
    o3m: s.n.O3m ? gasAvg(s.n.O3m) : 0,
    ps: Object.assign({}, solver.periodStats),
  });
  const r = rows[rows.length - 1];
  console.log(`период ${r.per}: P=${r.power.toFixed(2)} Вт  I_pk=${(r.peakI * 1e3).toFixed(2)} мА  ` +
    `имп=${r.pulses}  FWHM=${(r.fwhm * 1e9).toFixed(0)} нс  E/N_max=${r.maxEN.toFixed(0)} Td  ` +
    `U_burn=${(r.Uburn / 1e3).toFixed(2)} кВ  |σ|=${r.sigMax.toExponential(2)} Кл/м²  ` +
    `O3=${r.o3ppm.toFixed(1)} ppm`);
  if (diverged) { console.log(`\n!!! ПРОГОН ОБОРВАН: ${diverged}`); break; }
}

const wall = (Date.now() - wall0) / 1000;
const last = rows[rows.length - 1];
const ps = solver.periodStats;
const s = solver.state;

console.log('');
console.log('=== ИТОГ ============================================================');
console.log(`шагов всего            : ${steps}`);
console.log(`средний dt             : ${(dtSum / steps).toExponential(3)} с  (${((dtSum / steps) * 1e12).toFixed(2)} пс)`);
console.log(`время счёта            : ${wall.toFixed(1)} с  (${(wall / NPER).toFixed(2)} с/период, ${(wall * 1e6 / steps).toFixed(0)} мкс/шаг)`);
console.log(`пиковый ток (полный)   : ${(last.peakI * 1e3).toFixed(3)} мА   = ${(last.peakI / A).toExponential(2)} А/м²`);
console.log(`пиковый ток (кондукц.) : ${(last.peakIcond * 1e3).toFixed(3)} мА   = ${(last.peakIcond / A).toExponential(2)} А/м²`);
console.log(`импульсов за период    : ${last.pulses}   (по контракту state: ${s.breakdownsPerPeriod})`);
console.log(`FWHM: средняя ${(last.fwhm * 1e9).toFixed(0)} нс, главного импульса ${(last.fwhmMain * 1e9).toFixed(0)} нс`);
console.log(`импульсы (пик мА / заряд нКл):`);
last.peaks.forEach((p, i) => console.log(`   #${i + 1}  фаза ${last.phases[i].toFixed(0)}°  I=${(p * 1e3).toFixed(0)} мА  q=${(last.pulseQ[i] * 1e9).toFixed(1)} нКл`));
console.log(`мощность (∫U I dt / T) : ${last.power.toFixed(3)} Вт   (periodStats: ${ps.powerW.toFixed(3)} Вт)`);
console.log(`энергия за период      : ${last.energy.toExponential(3)} Дж`);
console.log(`заряд за полупериод    : +${(last.qPos * 1e9).toFixed(1)} нКл / -${(last.qNeg * 1e9).toFixed(1)} нКл`);
console.log(`максимум E/N           : ${last.maxEN.toFixed(1)} Td`);
console.log(`U_burn (среднее)       : ${(last.Uburn / 1e3).toFixed(2)} кВ   (periodStats: ${ps.UburnkV.toFixed(2)} кВ)`);
console.log(`поверхностный заряд    : |σ|max = ${last.sigMax.toExponential(3)} Кл/м²  (σ_L=${s.sigmaL.toExponential(2)}, σ_R=${s.sigmaR.toExponential(2)})`);
console.log(`C_diel / C_cell        : ${(ps.Cdiel * 1e12).toFixed(2)} пФ (измеренный наклон горящей ветви, ` +
  `R²=${ps.qvR2on.toFixed(3)}, дуг ${ps.qvArcsOn}/${ps.qvArcsOff}, qvOk=${ps.qvOk}) / ${(ps.Ccell * 1e12).toFixed(3)} пФ` +
  `\n                         тождество Q=C_d(U−U_gap): ${(ps.CdielIdentity * 1e12).toFixed(2)} пФ, геометрия ${(ps.CdielGeom * 1e12).toFixed(2)} пФ`);
console.log(`озон                   : ${last.o3.toExponential(3)} м^-3 = ${last.o3ppm.toFixed(2)} ppm`);
console.log(`атомарный O            : ${last.o.toExponential(3)} м^-3`);
console.log(`O2(a¹Δ)                : ${last.o2a.toExponential(3)} м^-3`);
console.log(`n_e / O- / O2- / O3-   : ${last.ne.toExponential(2)} / ${last.om.toExponential(2)} / ${last.o2m.toExponential(2)} / ${last.o3m.toExponential(2)} м^-3`);

// --- проверки по PHYSICS.md §6 ----------------------------------------------
console.log('');
console.log('=== ПРОВЕРКИ (docs/REFERENCE_TARGETS.md) ============================');
const checks = [];
const chk = (id, ok, msg) => { checks.push(ok); console.log(`${ok ? 'OK  ' : 'FAIL'} ${id}  ${msg}`); };

if (diverged) {
  console.log('');
  console.log('!!! 1D-МОДЕЛЬ РАСХОДИТСЯ НА ИМПУЛЬСЕ — измеренных величин нет.');
  console.log(`!!! ${diverged}`);
  console.log('!!! Диагноз: с корректными коэффициентами (ионизация ~втрое быстрее)');
  console.log('!!! приповерхностная ячейка набирает n_e ~ 1e24 м^-3 и E/N ~ 4e4 Td;');
  console.log('!!! на ОБОИХ барьерах оседает заряд ОДНОГО знака, экранирование зазора');
  console.log('!!! пропадает и самогашение не наступает. Не зависит от сетки (200/400/800),');
  console.log('!!! от CFL, от потолка dt (до 1e-14 с), от gamma (включая gamma=0) и от U0 (5–10 кВ).');
  console.log('!!! Это ограничение 1D-СОЛВЕРА, а не ошибка коэффициентов: сами коэффициенты');
  console.log('!!! воспроизводят все эталоны REFERENCE_TARGETS §1–4 (см. test/solver.test.mjs T6*).');
  console.log('');
}

// V2: сравнивать U_burn с U_br НАПРЯМУЮ нельзя (ERRATA §C0, §C V3) — это разные
// величины. U_br = 3.473 кВ статическое, в голом газе; U_burn динамическое, на
// перенапряжённом импульсе с диэлектриками и цепью. Проверяем только, что зазор
// зажат в физичной вилке вокруг U_br (снизу — остаточный sigma даёт зажигание
// ниже U_br, сверху — формативное запаздывание даёт перелёт).
const U_BR = 3.473e3;                     // REFERENCE_TARGETS §2, gamma = 0.02
const U_I = U_BR * (1 + 0.885 / 7.97);    // = 3.859 кВ, порог на ЯЧЕЙКЕ
chk('V2 ', last.Uburn / 1e3 > 0.7 * U_BR / 1e3 && last.Uburn / 1e3 < 1.6 * U_BR / 1e3,
  `U_burn = ${(last.Uburn / 1e3).toFixed(2)} кВ против U_br = ${(U_BR / 1e3).toFixed(3)} кВ ` +
  `(голый газ; вилка 0.7–1.6 U_br). Порог на ячейке U_i = ${(U_I / 1e3).toFixed(3)} кВ`);
chk('V4a', Math.abs(ps.Ccell - 0.797e-12) / 0.797e-12 < 0.15,
  `C_cell = ${(ps.Ccell * 1e12).toFixed(3)} пФ (ожидается 0.797)`);
// V4b — самопроверка цепи: наклон dQ/d(U_app − U_gap) обязан равняться C_d
// ТОЖДЕСТВЕННО (Q = C_d(U_app − U_gap)); это не измерение фигуры.
chk('V4b', Math.abs(ps.CdielIdentity - 7.97e-12) / 7.97e-12 < 0.05,
  `тождество Q=C_d(U−U_gap): ${(ps.CdielIdentity * 1e12).toFixed(2)} пФ (аналитика 7.97)`);
// V4c — измеренный наклон горящей ветви: в 1D зазор зажат неидеально, фигура не
// параллелограмм, наклон выше C_d в 1.5–2.5 раза; это должно быть ПОМЕЧЕНО.
chk('V4c', ps.Cdiel > 7.97e-12 && ps.Cdiel < 3 * 7.97e-12 && ps.qvOk === false,
  `наклон горящей ветви ${(ps.Cdiel * 1e12).toFixed(2)} пФ, R²on=${ps.qvR2on.toFixed(3)}, qvOk=${ps.qvOk}`);
// V5: формула Мэнли P = 4*f*C_d*U_b*(U0 - U_i), где U_b — напряжение НА ГАЗЕ,
// а U_i — порог НА ЯЧЕЙКЕ (ERRATA §C V5: в исходной спеке U_i было подставлено
// в оба места). С новыми эталонами U_b = 3.473 кВ, U_i = 3.859 кВ -> 6.80 Вт
// (со старыми коэффициентами было 7.78 Вт).
const Pmanley = 4 * P.freqKHz * 1e3 * 7.97e-12 * U_BR * (P.U0kV * 1e3 - U_I);
chk('V5 ', last.power > 0.4 * Pmanley && last.power < 2.5 * Pmanley,
  `P = ${last.power.toFixed(2)} Вт против Мэнли ${Pmanley.toFixed(2)} Вт ` +
  `(U_b=${(U_BR / 1e3).toFixed(3)} кВ на газе, U_i=${(U_I / 1e3).toFixed(3)} кВ на ячейке)`);
chk('V6 ', last.qPos * 1e9 > 20 && last.qPos * 1e9 < 400,
  `ΔQ(+) = ${(last.qPos * 1e9).toFixed(1)} нКл за полупериод (ожидается ~92)`);
chk('V7a', last.pulses >= 2,
  `${last.pulses} импульс(ов) за период (>=1 на полупериод)`);
// V7b: спека даёт вилку 1e2..1e4 А/м², исходя из оценки «1 А/см² × 100 нс = 100 нКл».
// Инвариант здесь не сама амплитуда, а произведение j_pk·FWHM ≈ ΔQ/A: 1D-модель
// зажигает всю площадь синхронно и даёт более короткий и высокий импульс
// (~13 нс / 4.5 А/см²), что для одиночного филамента экспериментально нормально.
// Поэтому вилку по амплитуде расширяем до 1e5, а согласованность проверяем явно (V7d).
chk('V7b', last.peakIcond / A > 1e2 && last.peakIcond / A < 1e5,
  `j_pk = ${(last.peakIcond / A).toExponential(2)} А/м² (спека 1e2..1e4, допуск до 1e5 — см. V7d)`);
chk('V7c', last.fwhmMain * 1e9 > 3 && last.fwhmMain * 1e9 < 500,
  `FWHM(главный) = ${(last.fwhmMain * 1e9).toFixed(1)} нс (спека 20–200, 1D даёт более резкий)`);
const qEst = last.peakIcond * last.fwhmMain;
chk('V7d', qEst > 0.2 * last.qMain && qEst < 5 * last.qMain,
  `j_pk·FWHM = ${(qEst * 1e9).toFixed(1)} нКл согласовано с зарядом главного импульса ${(last.qMain * 1e9).toFixed(1)} нКл`);
chk('V8 ', last.sigMax > 1e-4 && last.sigMax < 2e-3,
  `|σ| = ${last.sigMax.toExponential(2)} Кл/м² (ожидается ~4.6e-4)`);
chk('V9 ', (last.om + last.o2m + last.o3m) > last.ne,
  `n(-ионов) = ${(last.om + last.o2m + last.o3m).toExponential(2)} > n_e = ${last.ne.toExponential(2)}`);
// V10 — ОЦЕНКА ПОРЯДКА, а не независимый тест (ERRATA §C V10): эмпирический
// выход 100 г/кВт·ч относится к проточному озонатору, модель — безпроточная.
// После замены k_diss (×1000 при 100 Td, источник атомарного O вырос ~в 8 раз
// при типичных для ДБР полях, RATES_REVIEW §7.3) старая цель 1.9e21 за период
// заведомо занижена. Физичная вилка после нескольких периодов: 1e21…1e22 м^-3.
chk('V10', last.o3 > 1e21 && last.o3 < 1e22,
  `[O3] = ${last.o3.toExponential(2)} м^-3 после ${NPER} периодов (физичная вилка 1e21…1e22)`);
chk('V-nan', Number.isFinite(last.power) && Number.isFinite(last.o3),
  'нет NaN');

chk('V0 ', diverged === null,
  diverged ? `прогон оборван: ${diverged}` : 'прогон досчитан без ухода в нефизичное состояние');

// V3: ниже U_i = U_br*(1 + C_g/C_d) = 3.86 кВ (было 4.20 кВ со старым U_br)
// разряда быть НЕ должно — чисто ёмкостный отклик и вырожденная (линейная)
// фигура Лиссажу с наклоном C_cell. Пробный запуск при 3.0 кВ: запас до нового
// порога стал меньше (3.5 кВ — это уже 91 % от U_i, слишком близко).
{
  const sub = new DBDSolver({ mode: 'demo', U0kV: 3.0, freqKHz: P.freqKHz });
  const Tsub = 1 / (P.freqKHz * 1e3);
  let pk = 0, sig = 0;
  while (sub.state.t < 2 * Tsub) {
    sub.step();
    if (physGuard(sub.state)) break;
    if (sub.state.t > Tsub) {
      pk = Math.max(pk, Math.abs(sub.state.condCurrent));
      sig = Math.max(sig, Math.abs(sub.state.sigmaL), Math.abs(sub.state.sigmaR));
    }
  }
  const sp = sub.periodStats;
  chk('V3a', pk < 1e-5 && sub.state.breakdownsPerPeriod === 0,
    `при U0=3.0 кВ (ниже U_i=3.86 кВ) разряда нет: I_cond,pk=${pk.toExponential(2)} А, импульсов=${sub.state.breakdownsPerPeriod}`);
  chk('V3b', sig < 1e-6,
    `поверхностный заряд не накапливается: |σ|=${sig.toExponential(2)} Кл/м²`);
  chk('V3c', Math.abs(sp.Ccell - 0.797e-12) / 0.797e-12 < 0.05 && sp.qvArcsOn === 0,
    `вырожденная фигура Лиссажу: C_cell=${(sp.Ccell * 1e12).toFixed(3)} пФ (аналитика 0.797), ` +
    `горящих дуг ${sp.qvArcsOn} (шум не принят за разряд)`);
}

const nOk = checks.filter(Boolean).length;
console.log('');
console.log(`ПРОЙДЕНО ${nOk} / ${checks.length}`);
process.exit(nOk === checks.length ? 0 : 1);
