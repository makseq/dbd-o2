// tests/physics.test.mjs — тесты физического модуля 2D-модели ДБР в O2.
// Запуск:  node sim2d/tests/physics.test.mjs
// Критерии F1..F7 — из задания (ERRATA §C V1/V2, §D, §D6-РАЗБОР, §E1).
// F8..F11 — реперы новых коэффициентов из docs/RATES_REVIEW.md.
//
// ⚠ ЭТАЛОНЫ F1/F2 ПЕРЕСЧИТАНЫ. Старые цели (101.3 Td, U_br = 3.76 кВ) выведены
// ВНУТРИ неверных фитов сечений и ОТМЕНЕНЫ (ERRATA §C0). Новые — в
// docs/REFERENCE_TARGETS.md. Подгонять коэффициенты под старые значения нельзя.

import * as P from '../physics2d.mjs';

// ---------------------------------------------------------------------------
// Минимальный тест-харнесс (без зависимостей)
// ---------------------------------------------------------------------------
let passed = 0;
let failed = 0;
const failures = [];

function check(name, ok, detail = '') {
  if (ok) {
    passed++;
    console.log(`  ok   ${name}${detail ? '  — ' + detail : ''}`);
  } else {
    failed++;
    failures.push(`${name}${detail ? '  — ' + detail : ''}`);
    console.log(`  FAIL ${name}${detail ? '  — ' + detail : ''}`);
  }
}

function group(title) {
  console.log(`\n${title}`);
}

function near(a, b, relTol) {
  return Math.abs(a - b) <= relTol * Math.abs(b);
}

const fmt = (v, d = 4) => (Number.isFinite(v) ? v.toPrecision(d) : String(v));

const N = P.N_REF;          // 2.446e25 м^-3 при 1 атм/300 К
const T = 300;
const D_GAP = 1e-3;         // зазор 1 мм
const GAMMA = 0.02;

// ===========================================================================
// F1. alpha = eta при E/N = 128 +/- 5 Td
//     ЭТАЛОН ПЕРЕСЧИТАН: старая цель 101.3 Td отменена (ERRATA §C0).
//     128 Td — расчёт по сечениям Phelps и IST-Lisbon (Biagi 132) и совпадает
//     с экспериментальным пробойным полем O2 30-32 кВ/см = 122-130 Td.
// ===========================================================================
group('F1. Точка alpha = eta (НОВЫЙ эталон: 128 +/- 5 Td)');
{
  const x = P.alphaEtaCrossEN(N);
  check('F1 пересечение alpha=eta в пределах 128 +/- 5 Td',
    Math.abs(x - 128) <= 5, `x = ${fmt(x)} Td`);

  // невязка в самой точке
  const res = P.alphaN(x) - P.etaN(x, N);
  check('F1 невязка alpha-eta в найденной точке пренебрежима',
    Math.abs(res) < 1e-28, `alpha-eta = ${fmt(res)} м^2`);

  // то же в виде критического поля: 31.3 кВ/см против эксперимента 30-32 кВ/см
  const Ecrit = x * P.TD * N;
  check('F1 критическое поле 30..32 кВ/см (эксперимент по пробою O2)',
    Ecrit >= 3.0e6 && Ecrit <= 3.2e6, `${fmt(Ecrit / 1e5)} кВ/см`);

  // старая цель ДОЛЖНА быть отвергнута — иначе кто-то вернул воздушные данные
  check('F1 отменённая цель 101.3 Td более не воспроизводится',
    Math.abs(x - 101.3) > 10, `x = ${fmt(x)} Td vs отменённые 101.3`);
}

// ===========================================================================
// F2. Самоподдержание по ЭЛЕКТРООТРИЦАТЕЛЬНОЙ формуле
//     gamma*alpha/(alpha-eta)*[exp((alpha-eta)*d)-1] = 1, d = 1 мм
//
//     ЭТАЛОН ПЕРЕСЧИТАН ЧИСЛЕННО с новыми alpha(E/N) и eta(E/N).
//     Старая цель (153.4 Td, 3.76 кВ) ОТМЕНЕНА.
//     ⚠ U_br УМЕНЬШИЛСЯ, а не вырос, вопреки ожиданию ERRATA §C0: критическое
//     поле выросло на 26 %, но alpha выше него растёт КРУЧЕ (alpha/N вдвое-втрое
//     больше фита Morrow), поэтому нужный запас над критическим полем меньше.
//     Проверено независимо по опубликованной таблице Phelps: 140.4 Td / 3.43 кВ.
// ===========================================================================
group('F2. Критерий самоподдержания (НОВЫЙ эталон: 142.0 Td, U_br = 3.473 кВ, 5%)');
{
  const x = P.breakdownEN(D_GAP, GAMMA, N);
  const E = x * P.TD * N;               // В/м
  const U = E * D_GAP;                  // В (однородное поле в зазоре)

  check('F2 E/N самоподдержания = 141.97 Td (5%)', near(x, 141.97, 0.05), `x = ${fmt(x)} Td`);
  check('F2 U_br = 3.473 кВ (5%)', near(U, 3473, 0.05), `U = ${fmt(U / 1000)} кВ`);

  // невязка критерия строго в корне
  const res = P.selfSustainResidual(x, D_GAP, GAMMA, N);
  check('F2 невязка критерия в корне < 1e-6', Math.abs(res) < 1e-6, `res = ${fmt(res)}`);

  // критерий должен быть монотонно растущим около корня
  check('F2 критерий не выполнен ниже корня и выполнен выше',
    P.selfSustainResidual(x * 0.98, D_GAP, GAMMA, N) < 0 &&
    P.selfSustainResidual(x * 1.02, D_GAP, GAMMA, N) > 0);

  // упрощённая формула даёт ДРУГОЙ ответ — фиксируем расхождение явно
  const target = Math.log(1 + 1 / GAMMA);
  let lo = 100, hi = 1000;
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (lo + hi);
    if (P.alphaEff(m, N) * D_GAP < target) lo = m; else hi = m;
  }
  const xSimple = 0.5 * (lo + hi);
  check('F2 упрощённая формула (alpha-eta)d=ln(1+1/g) даёт БОЛЬШЕЕ поле (её не используем)',
    xSimple > x, `упрощ. = ${fmt(xSimple)} Td против ${fmt(x)} Td`);

  // ПАРАМЕТРИЧЕСКАЯ зависимость U_br(gamma) — ERRATA §C1: gamma НЕ калибруется,
  // а предъявляется как параметр. Таблица — docs/REFERENCE_TARGETS.md.
  const GAMMA_TARGETS = [
    [0.005, 148.56, 3634],
    [0.010, 145.41, 3557],
    [0.020, 141.97, 3473],
    [0.050, 136.88, 3349],
  ];
  for (const [g, xRef, uRef] of GAMMA_TARGETS) {
    const xg = P.breakdownEN(D_GAP, g, N);
    const ug = xg * P.TD * N * D_GAP;
    check(`F2 gamma=${g}: E/N = ${xRef} Td, U_br = ${(uRef / 1000).toFixed(3)} кВ (5%)`,
      near(xg, xRef, 0.05) && near(ug, uRef, 0.05),
      `${fmt(xg)} Td, ${fmt(ug / 1000)} кВ`);
  }
  // монотонность: больше gamma -> легче поджиг -> меньше U_br
  let mono = true;
  for (let i = 1; i < GAMMA_TARGETS.length; i++) {
    if (!(P.breakdownEN(D_GAP, GAMMA_TARGETS[i][0], N)
        < P.breakdownEN(D_GAP, GAMMA_TARGETS[i - 1][0], N))) mono = false;
  }
  check('F2 U_br монотонно убывает с ростом gamma', mono);

  // самоподдержание обязано лежать ВЫШЕ критического поля alpha = eta
  check('F2 поле самоподдержания > поля alpha=eta',
    x > P.alphaEtaCrossEN(N), `${fmt(x)} Td > ${fmt(P.alphaEtaCrossEN(N))} Td`);
}

// ===========================================================================
// F3. Все константы скорости конечны и неотрицательны на сетке 0.1..1000 Td
// ===========================================================================
group('F3. Конечность/неотрицательность на сетке E/N = 0.1..1000 Td');
{
  const grid = [];
  for (let i = 0; i <= 400; i++) {
    grid.push(0.1 * Math.pow(10000, i / 400));    // лог-сетка 0.1..1000 Td
  }
  let bad = 0, badWhere = '';
  let zeroAbove20 = 0, zeroWhere = '';
  for (const x of grid) {
    for (const r of P.REACTIONS) {
      const k = r.rate(x, N, T);
      if (!Number.isFinite(k) || k < 0) {
        bad++;
        if (!badWhere) badWhere = `${r.id} @ ${fmt(x)} Td = ${k}`;
      }
      // выше 20 Td открыты все каналы — нулей быть не должно
      if (x >= 20 && k === 0) {
        zeroAbove20++;
        if (!zeroWhere) zeroWhere = `${r.id} @ ${fmt(x)} Td`;
      }
    }
    // транспорт электронов
    for (const [nm, v] of [['muE_N', P.muE_N(x)], ['DE_N', P.DE_N(x)],
      ['meanEnergy', P.meanEnergy(x)], ['alphaN', P.alphaN(x)],
      ['etaN', P.etaN(x, N)], ['kIoniz', P.kIoniz(x)]]) {
      if (!Number.isFinite(v) || v < 0) {
        bad++;
        if (!badWhere) badWhere = `${nm} @ ${fmt(x)} Td = ${v}`;
      }
    }
  }
  check('F3 нет NaN/Inf/отрицательных значений', bad === 0, badWhere || `${grid.length} узлов`);
  check('F3 при E/N >= 20 Td все каналы строго положительны',
    zeroAbove20 === 0, zeroWhere || 'ok');

  // граничные и патологические аргументы
  const edge = [0, 1e-12, 1e-6, 0.1, 1000, 5000];
  let badEdge = '';
  for (const x of edge) {
    for (const r of P.REACTIONS) {
      const k = r.rate(x, N, T);
      if (!Number.isFinite(k) || k < 0) badEdge ||= `${r.id} @ ${x} = ${k}`;
    }
    if (!Number.isFinite(P.muE_N(x)) || !Number.isFinite(P.meanEnergy(x))) badEdge ||= `transport @ ${x}`;
  }
  check('F3 краевые аргументы (0, 1e-12, 5000 Td) не ломают фиты', badEdge === '', badEdge);
}

// ===========================================================================
// F4. Модуль поля: k(-x) == k(+x) (ERRATA §E1)
// ===========================================================================
group('F4. Модуль приведённого поля (ERRATA E1)');
{
  const grid = [0.1, 1, 5, 20, 50, 101.3, 150, 153.4, 300, 700, 1000];
  let mismatch = '';
  for (const x of grid) {
    for (const r of P.REACTIONS) {
      const kp = r.rate(x, N, T);
      const km = r.rate(-x, N, T);
      if (!(kp === km)) mismatch ||= `${r.id} @ ${x}: ${kp} vs ${km}`;
    }
    const pairs = [
      ['muE_N', P.muE_N(x), P.muE_N(-x)],
      ['DE_N', P.DE_N(x), P.DE_N(-x)],
      ['meanEnergy', P.meanEnergy(x), P.meanEnergy(-x)],
      ['alphaN', P.alphaN(x), P.alphaN(-x)],
      ['etaN', P.etaN(x, N), P.etaN(-x, N)],
      ['kIoniz', P.kIoniz(x), P.kIoniz(-x)],
      ['electronTemp', P.electronTemp(x), P.electronTemp(-x)],
    ];
    for (const [nm, a, b] of pairs) if (a !== b) mismatch ||= `${nm} @ ${x}: ${a} vs ${b}`;
  }
  check('F4 все коэффициенты чётны по E/N (побитово)', mismatch === '', mismatch || `${grid.length} точек`);

  // reducedField тоже обязан брать модуль
  check('F4 reducedField(-E) == reducedField(+E)',
    P.reducedField(-3e6, N) === P.reducedField(3e6, N),
    `${fmt(P.reducedField(3e6, N))} Td при 3e6 В/м`);

  // kD6 от знака поля не зависит
  check('F4 kD6 чётна по E', P.kD6(-3e6, 2.37e-4, T) === P.kD6(3e6, 2.37e-4, T));
}

// ===========================================================================
// F5. Размерности: m3/s в 1e-22..1e-13, m6/s в 1e-46..1e-40
//     (документированные исключения помечены полем dimNote)
// ===========================================================================
group('F5. Размерности констант скорости');
{
  const GENERIC = { 'm3/s': [1e-22, 1e-13], 'm6/s': [1e-46, 1e-40] };
  // «унитарное окно» — грубая проверка, что нигде не потерян фактор 1e-6/1e-12
  const UNIT_SANITY = { 'm3/s': [1e-31, 1e-11], 'm6/s': [1e-47, 1e-39] };

  const grid = [];
  for (let i = 0; i <= 200; i++) grid.push(0.1 * Math.pow(10000, i / 200));

  for (const r of P.REACTIONS) {
    const [lo, hi] = r.kRange;
    let kmin = Infinity, kmax = -Infinity;
    for (const x of grid) {
      const k = r.rate(x, N, T);
      if (k < kmin) kmin = k;
      if (k > kmax) kmax = k;
    }
    const inRange = kmin >= lo && kmax <= hi;
    check(`F5 ${r.id} k в объявленном диапазоне [${lo}, ${hi}] ${r.kUnits}`,
      inRange, `k = ${fmt(kmin, 3)} .. ${fmt(kmax, 3)}`);

    // объявленный верхний предел обязан лежать внутри родового окна,
    // если у реакции нет явного обоснования отклонения (dimNote)
    const g = GENERIC[r.kUnits];
    const u = UNIT_SANITY[r.kUnits];
    if (r.dimNote) {
      check(`F5 ${r.id} исключение обосновано и внутри унитарного окна`,
        hi <= u[1] && (lo === 0 || lo >= u[0]), r.dimNote);
    } else {
      check(`F5 ${r.id} максимум внутри родового окна ${r.kUnits}`,
        kmax <= g[1] && kmax > 0, `kmax = ${fmt(kmax, 3)}`);
    }
  }

  // отдельная ловушка на «см вместо м»: три реперных значения
  check('F5 репер: k_ii(1 атм) = 5.1e-12 м^3/с (а не 2e-6 см^3/с = 2e-12)',
    near(P.kIonIon(N), 5.1e-12, 0.03), `${fmt(P.kIonIon(N))} м^3/с`);
  check('F5 репер: kC1(300) = 2.4e-42 м^6/с', near(P.kC1(300), 2.4e-42, 1e-9));
  check('F5 репер: kN1(300) = 6.0e-46 м^6/с', near(P.kN1(300), 6.0e-46, 1e-9));

  // характерные времена — они мгновенно ловят потерю степени N (PHYSICS.md §6 п.11)
  const tauCluster = 1 / (P.kC1(300) * N * N);
  const tauO3 = 1 / (P.kN1(300) * N * N);
  const tauAtt3 = 1 / (P.kAtt3(1) * N * N);
  check('F5 время O2+ -> O4+ ~ 0.7 нс', tauCluster > 3e-10 && tauCluster < 2e-9,
    `${fmt(tauCluster * 1e9, 3)} нс`);
  check('F5 время O -> O3 ~ 2.8 мкс', tauO3 > 1e-6 && tauO3 < 6e-6,
    `${fmt(tauO3 * 1e6, 3)} мкс`);
  // ЭТАЛОН ПЕРЕСЧИТАН: k3(1 Td) = 2.935e-42 м^6/с -> tau = 0.57 нс.
  // Допуск отражает ЧЕСТНУЮ неопределённость k3 (+/-40 %), а не подгонку.
  check('F5 время трёхтельного прилипания при 1 Td = 0.57 нс (+/-40 % по k3)',
    tauAtt3 > 0.41e-9 && tauAtt3 < 0.95e-9, `${fmt(tauAtt3 * 1e9, 3)} нс`);
}

// ===========================================================================
// F6. D6: обе конвенции массы (ERRATA §D6-РАЗБОР), 30 кВ/см, допуск 2%
// ===========================================================================
group('F6. D6 — обе конвенции массы (ERRATA D6-РАЗБОР, допуск 2%)');
{
  const E = 3.0e6;                 // 30 кВ/см
  const muO2m = 2.37e-4;           // ERRATA §D, 1 атм/300 К
  const vd = muO2m * E;
  check('F6 v_d(O2-) = 711 м/с при 30 кВ/см', near(vd, 711, 0.02), `${fmt(vd)} м/с`);

  const TeffN = P.wannierTemp(vd, T, 'neutral');
  const nuN = P.nuD6(E, muO2m, T, N, 'neutral');
  check('F6 neutral: T_eff = 948 K', near(TeffN, 948, 0.02), `${fmt(TeffN)} K`);
  check('F6 neutral: nu = 3.24e7 1/с', near(nuN, 3.24e7, 0.02), `${fmt(nuN)} 1/с`);

  const TeffR = P.wannierTemp(vd, T, 'reduced');
  const nuR = P.nuD6(E, muO2m, T, N, 'reduced');
  check('F6 reduced: T_eff = 624 K', near(TeffR, 624, 0.02), `${fmt(TeffR)} K`);
  check('F6 reduced: nu = 1.23e6 1/с', near(nuR, 1.23e6, 0.02), `${fmt(nuR)} 1/с`);

  check('F6 расхождение конвенций x26 при 30 кВ/см (зафиксированная неопределённость)',
    near(nuN / nuR, 26.3, 0.05), `x${fmt(nuN / nuR, 3)}`);

  // вторая строка таблицы ERRATA: 37.8 кВ/см -> 896 м/с
  const E2 = 3.78e6;
  const vd2 = muO2m * E2;
  check('F6 v_d = 896 м/с при 37.8 кВ/см', near(vd2, 896, 0.02), `${fmt(vd2)} м/с`);
  check('F6 neutral @37.8 кВ/см: T_eff = 1330 K',
    near(P.wannierTemp(vd2, T, 'neutral'), 1330, 0.02),
    `${fmt(P.wannierTemp(vd2, T, 'neutral'))} K`);
  check('F6 neutral @37.8 кВ/см: nu = 2.08e8 1/с',
    near(P.nuD6(E2, muO2m, T, N, 'neutral'), 2.08e8, 0.03),
    `${fmt(P.nuD6(E2, muO2m, T, N, 'neutral'))} 1/с`);
  check('F6 reduced @37.8 кВ/см: T_eff = 815 K',
    near(P.wannierTemp(vd2, T, 'reduced'), 815, 0.02),
    `${fmt(P.wannierTemp(vd2, T, 'reduced'))} K`);
  check('F6 reduced @37.8 кВ/см: nu = 1.14e7 1/с',
    near(P.nuD6(E2, muO2m, T, N, 'reduced'), 1.14e7, 0.03),
    `${fmt(P.nuD6(E2, muO2m, T, N, 'reduced'))} 1/с`);

  // дефолт обязан быть 'neutral'
  check('F6 дефолтная конвенция — neutral (ERRATA §D6)',
    P.TUNING_DEFAULTS.d6MassConvention === 'neutral' &&
    P.kD6(E, muO2m, T) === P.kD6(E, muO2m, T, 'neutral'));

  // переключение через TUNING работает
  P.setTuning({ d6MassConvention: 'reduced' });
  const switched = P.kD6(E, muO2m, T) === P.kD6(E, muO2m, T, 'reduced');
  P.resetTuning();
  check('F6 TUNING.d6MassConvention переключает ветку', switched);
}

// ===========================================================================
// F7. Сохранение заряда в КАЖДОЙ реакции — точно
// ===========================================================================
group('F7. Стехиометрия по заряду');
{
  let bad = '';
  for (const r of P.REACTIONS) {
    const q = P.chargeBalance(r);
    if (q !== 0) bad ||= `${r.id}: ${r.equation} -> dq = ${q}`;
  }
  check('F7 суммарный заряд сохраняется во всех реакциях', bad === '',
    bad || `${P.REACTIONS.length} реакций`);

  // и заодно: все участники реакций объявлены в SPECIES
  let unknown = '';
  for (const r of P.REACTIONS) {
    for (const s of [...r.reagents, ...r.products]) {
      if (!P.SPECIES_BY_ID[s]) unknown ||= `${r.id}: ${s}`;
    }
  }
  check('F7 все сорта из реакций объявлены в SPECIES', unknown === '', unknown);

  // сохранение числа ядер O (фон O2 участвует явно, поэтому баланс обязан сходиться)
  const nO = { e: 0, O2p: 2, O4p: 4, Om: 1, O2m: 2, O3m: 3, O: 1, O3: 3, O2a: 2, O2: 2 };
  let badO = '';
  for (const r of P.REACTIONS) {
    let d = 0;
    for (const s of r.reagents) d -= nO[s];
    for (const s of r.products) d += nO[s];
    if (d !== 0) badO ||= `${r.id}: ${r.equation} -> dN(O) = ${d}`;
  }
  check('F7 (доп.) баланс атомов кислорода в каждой реакции', badO === '', badO);
}

// ===========================================================================
// Дополнительно: набор частиц и подвижности (ERRATA §D)
// ===========================================================================
group('Дополнительно: SPECIES и подвижности (ERRATA D)');
{
  const need = ['e', 'O2p', 'O4p', 'Om', 'O2m', 'O3m', 'O', 'O3', 'O2a'];
  check('набор частиц полный (9 сортов + фон)',
    need.every((id) => P.SPECIES_BY_ID[id]) && P.SPECIES.length === 10,
    P.SPECIES.map((s) => s.id).join(','));

  for (const [id, expected] of [['O2p', 2.41e-4], ['Om', 3.51e-4], ['O2m', 2.37e-4]]) {
    const mu = P.ionMobility(id, N);
    check(`подвижность ${id} = ${expected} м^2/(В*с) при 1 атм/300 К (ERRATA §D)`,
      near(mu, expected, 0.01), `${fmt(mu)}`);
  }

  // масштабирование как 1/N
  const mu1 = P.ionMobility('O2p', N);
  const mu2 = P.ionMobility('O2p', 2 * N);
  check('подвижность масштабируется как 1/N', near(mu2, mu1 / 2, 1e-12),
    `${fmt(mu1)} -> ${fmt(mu2)}`);

  // диффузия ионов по Эйнштейну
  const Di = P.ionDiffusion('O2p', N, T);
  check('D(O2+) = mu*kT/e = 6.2e-6 м^2/с', near(Di, 6.23e-6, 0.02), `${fmt(Di)}`);

  // время пролёта частиц через зазор при НОВОМ поле самоподдержания (142 Td)
  const Ebr = 141.97 * P.TD * N;                 // 3.47e6 В/м
  const tIon = D_GAP / (P.ionMobility('O2p', N) * Ebr);
  check('время пролёта иона ~1.2 мкс', tIon > 0.7e-6 && tIon < 1.6e-6, `${fmt(tIon * 1e6, 3)} мкс`);
  // ЭТАЛОН ПЕРЕСЧИТАН: w_e(142 Td) = 2.31e5 м/с (было 1.51e5 на воздушной таблице)
  const tEl = D_GAP / P.driftElectron(141.97);
  check('время пролёта электрона ~4.3 нс (чистый O2, не воздух)',
    tEl > 3e-9 && tEl < 6e-9, `${fmt(tEl * 1e9, 3)} нс`);
}

// ===========================================================================
// F8. Диссоциация e + O2 -> 2O + e (RATES_REVIEW §2.2), допуск 10 %
//     Состав суммы: каналы 6.0 эВ + 8.4 эВ, полосы Герцберга 4.5 эВ НЕ входят.
// ===========================================================================
group('F8. k_diss (каналы 6.0 + 8.4 эВ), допуск 10 %');
{
  const REF = [[50, 4.729e-16], [100, 1.937e-15], [200, 5.877e-15], [400, 1.245e-14]];
  for (const [x, k] of REF) {
    check(`F8 k_diss(${x} Td) = ${k.toExponential(3)} м^3/с (10%)`,
      near(P.kDiss(x), k, 0.10), `${fmt(P.kDiss(x))} м^3/с`);
  }
  // старое значение обязано быть отвергнуто: занижение было ровно в ~1000 раз
  check('F8 старое значение ~1.8e-18 при 100 Td более не воспроизводится',
    P.kDiss(100) / 1.8e-18 > 100, `отношение ×${fmt(P.kDiss(100) / 1.8e-18, 3)}`);
  // монотонный рост в рабочем диапазоне
  let mono = true;
  for (let x = 20; x < 900; x *= 1.1) if (!(P.kDiss(x * 1.1) > P.kDiss(x))) mono = false;
  check('F8 k_diss монотонно растёт на 20..1000 Td', mono);
}

// ===========================================================================
// F9. Возбуждение O2(a1Dg) (RATES_REVIEW §3.1), допуск 15 %
// ===========================================================================
group('F9. k(O2 a1Dg), допуск 15 %');
{
  check('F9 k(O2a) при 100 Td = 7.0e-16 м^3/с (15%)',
    near(P.kExcO2a(100), 7.0e-16, 0.15), `${fmt(P.kExcO2a(100))} м^3/с`);
  check('F9 k(O2a) при 50 Td = 5.394e-16 м^3/с (15%)',
    near(P.kExcO2a(50), 5.394e-16, 0.15), `${fmt(P.kExcO2a(50))} м^3/с`);
  // насыщение 9.6e-16 при 400-600 Td и слабый спад выше
  check('F9 насыщение 9.6e-16 при 500 Td (15%)',
    near(P.kExcO2a(500), 9.6e-16, 0.15), `${fmt(P.kExcO2a(500))} м^3/с`);
  check('F9 выше насыщения — слабый спад (ФРЭЭ уходит за максимум сечения)',
    P.kExcO2a(1000) < P.kExcO2a(500), `${fmt(P.kExcO2a(1000))} < ${fmt(P.kExcO2a(500))}`);
  check('F9 старая оценка ~2e-17 при 100 Td более не воспроизводится',
    P.kExcO2a(100) / 2e-17 > 10, `отношение ×${fmt(P.kExcO2a(100) / 2e-17, 3)}`);
}

// ===========================================================================
// F10. Трёхтельное прилипание k3 (Chanin-Phelps-Biondi 1962), допуск 40 %
//      Допуск ЧЕСТНЫЙ: расчёт по сечению Phelps даёт 1.76e-42 в тепловом
//      пределе против измеренных 2.8e-42 (ERRATA §D4). Не ужимать.
// ===========================================================================
group('F10. k3 трёхтельное прилипание, допуск 40 % (реальная неопределённость)');
{
  const kTh = P.kAtt3(0);            // тепловой предел (аргумент клипуется вниз)
  check('F10 k3 тепловое = 2.8e-42 м^6/с (40%)', near(kTh, 2.8e-42, 0.40),
    `${fmt(kTh)} м^6/с`);
  check('F10 максимум ~5.0e-42 м^6/с при ~0.03 Td (CPB 1962)',
    near(P.kAtt3(0.03), 5.0e-42, 0.40), `${fmt(P.kAtt3(0.03))} м^6/с`);
  // Немонотонность CPB: от теплового 2.8e-42 вверх к максимуму 5.0e-42
  // при <eps> ~ 0.09 эВ (E/N ~ 0.03 Td), затем спад. При 1 Td кривая уже
  // опустилась почти до теплового уровня, но ЕЩЁ выше него — это нормально.
  check('F10 немонотонность: k3(тепл.) < k3(0.03 Td, максимум) > k3(1 Td) > k3(10 Td)',
    kTh < P.kAtt3(0.03) && P.kAtt3(0.03) > P.kAtt3(1) && P.kAtt3(1) > P.kAtt3(10),
    `${fmt(kTh)} < ${fmt(P.kAtt3(0.03))} > ${fmt(P.kAtt3(1))} > ${fmt(P.kAtt3(10))}`);
  check('F10 k3(100 Td) = 8.35e-44 м^6/с (40%)', near(P.kAtt3(100), 8.347e-44, 0.40),
    `${fmt(P.kAtt3(100))} м^6/с`);
  // спад ВЫШЕ 25 Td степенной, не экспоненциальный: k3(200)/k3(100) должно быть
  // ~0.62, а exp(-EN/60) дало бы 0.19 — это и был источник занижения в 3-5 раз
  const decay = P.kAtt3(200) / P.kAtt3(100);
  check('F10 спад степенной, а не exp(-EN/60): k3(200)/k3(100) ~ 0.62',
    decay > 0.5 && decay < 0.75, `${fmt(decay, 3)} (exp-фит дал бы 0.19)`);
  check('F10 старое тепловое 3.0e-43 более не воспроизводится',
    kTh / 3.0e-43 > 5, `отношение ×${fmt(kTh / 3.0e-43, 3)}`);
}

// ===========================================================================
// F11. Транспорт электронов в ЧИСТОМ O2 (RATES_REVIEW §7.4).
//      Ловит возврат «воздушной» таблицы: она даёт 1.15e5 м/с при 100 Td.
// ===========================================================================
group('F11. Дрейф и средняя энергия электронов в чистом O2');
{
  check('F11 w_e(100 Td) = 1.84e5 м/с (10%)', near(P.driftElectron(100), 1.84e5, 0.10),
    `${fmt(P.driftElectron(100))} м/с`);
  check('F11 воздушное значение 1.15e5 м/с при 100 Td ОТВЕРГНУТО',
    P.driftElectron(100) / 1.15e5 > 1.4, `отношение ×${fmt(P.driftElectron(100) / 1.15e5, 3)}`);
  check('F11 <eps>(100 Td) = 4.362 эВ (10%)', near(P.meanEnergy(100), 4.362, 0.10),
    `${fmt(P.meanEnergy(100))} эВ`);
  check('F11 <eps>(50 Td) = 3.424 эВ (10%)', near(P.meanEnergy(50), 3.424, 0.10),
    `${fmt(P.meanEnergy(50))} эВ`);
  check('F11 <eps>(200 Td) = 6.254 эВ (10%)', near(P.meanEnergy(200), 6.254, 0.10),
    `${fmt(P.meanEnergy(200))} эВ`);
  check('F11 mu*N(100 Td) = 1.830e24 1/(м*В*с) (10%)', near(P.muE_N(100), 1.830e24, 0.10),
    `${fmt(P.muE_N(100))}`);
  // D/mu при 100 Td: 3.41 эВ против 3.373 эВ в опубликованной таблице Phelps
  const Dmu = P.DE_N(100) / P.muE_N(100);
  check('F11 D/mu(100 Td) = 3.37 эВ по таблице Phelps (10%)', near(Dmu, 3.373, 0.10),
    `${fmt(Dmu)} эВ`);
  // w_e монотонно растёт по E/N (в O2 нет отрицательной дифференциальной подвижности)
  let mono = true;
  for (let x = 5; x < 900; x *= 1.2) if (!(P.driftElectron(x * 1.2) > P.driftElectron(x))) mono = false;
  check('F11 w_e монотонно растёт на 5..1000 Td', mono);
}

// ===========================================================================
// F12. Ионизация: замена фита Morrow (ERRATA §D3 п.3)
// ===========================================================================
group('F12. alpha/N — таблица вместо фита Morrow');
{
  const morrow = (x) => (x < 150 ? 6.619e-21 * Math.exp(-559.3 / x)
                                 : 2.0e-20 * Math.exp(-724.8 / x));
  const r200 = P.alphaN(200) / morrow(200);
  const r500 = P.alphaN(500) / morrow(500);
  check('F12 расхождение с фитом Morrow при 200 Td ~ ×2 (ERRATA §D3 п.3)',
    r200 > 1.7 && r200 < 2.4, `×${fmt(r200, 3)}`);
  check('F12 расхождение с фитом Morrow при 500 Td ~ ×2',
    r500 > 1.7 && r500 < 2.4, `×${fmt(r500, 3)}`);
  // сверка с опубликованной таблицей Phelps (Report 28): k_ion/w_e
  // при 150 Td = 1.04e-16 / 2.4135e5 = 4.31e-22 м^2 (расчёт консервативно ниже)
  const a150 = P.alphaN(150);
  check('F12 alpha/N(150 Td) согласуется с таблицей Phelps 4.31e-22 м^2 (15%)',
    near(a150, 4.31e-22, 0.15), `${fmt(a150)} м^2`);
  check('F12 alpha/N монотонно растёт', P.alphaN(50) < P.alphaN(100)
    && P.alphaN(100) < P.alphaN(200) && P.alphaN(200) < P.alphaN(500));
}

// ===========================================================================
// Дополнительно: масштабирующие множители TUNING
// ===========================================================================
group('Дополнительно: kDissScale / kO2aScale / kAtt3Scale');
{
  const x = 150;
  const kd0 = P.kDiss(x), ke0 = P.kExcO2a(x), ka0 = P.kAtt3(x);
  check('E4 при 150 Td = 3.856e-15 м^3/с (RATES_REVIEW §2.2)',
    near(kd0, 3.856e-15, 0.10), `${fmt(kd0)} м^3/с`);
  check('E5 при 150 Td = 7.981e-16 м^3/с (RATES_REVIEW §3.1)',
    near(ke0, 7.981e-16, 0.15), `${fmt(ke0)} м^3/с`);

  P.setTuning({ kDissScale: 2.5, kO2aScale: 0.4, kAtt3Scale: 1.4 });
  const ok = near(P.kDiss(x), 2.5 * kd0, 1e-12) && near(P.kExcO2a(x), 0.4 * ke0, 1e-12)
    && near(P.kAtt3(x), 1.4 * ka0, 1e-12);
  P.resetTuning();
  check('масштабы применяются линейно и сбрасываются', ok &&
    P.kDiss(x) === kd0 && P.kExcO2a(x) === ke0 && P.kAtt3(x) === ka0);

  // kAtt3Scale — рабочий инструмент прогона по неопределённости +/-40 %
  P.setTuning({ kAtt3Scale: 1.4 });
  const xHi = P.alphaEtaCrossEN(N);
  P.setTuning({ kAtt3Scale: 0.6 });
  const xLo = P.alphaEtaCrossEN(N);
  P.resetTuning();
  check('чувствительность alpha=eta к k3 +/-40 % мала (k3 при 128 Td даёт ~4% eta)',
    Math.abs(xHi - xLo) < 3, `${fmt(xLo)} .. ${fmt(xHi)} Td`);

  let threw = false;
  try { P.setTuning({ nonsense: 1 }); } catch { threw = true; }
  check('setTuning отвергает неизвестные ключи', threw);
}

// ===========================================================================
console.log(`\n=== ИТОГО: ${passed} passed, ${failed} failed ===`);
if (failed) {
  console.log('Провалено:');
  for (const f of failures) console.log('  - ' + f);
  process.exitCode = 1;
}
