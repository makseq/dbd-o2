// physics2d.mjs — химия и транспорт для 2D осесимметричной fluid-модели ДБР в чистом O2.
//
// Источник данных, в порядке приоритета:
//   1. docs/RATES_REVIEW.md — электронные константы и транспорт (таблицы Больцмана)
//   2. docs/ERRATA.md §C0/§D — исправления к PHYSICS.md
//   3. docs/PHYSICS.md — всё остальное (ионная/нейтральная химия, подвижности)
// Где PHYSICS.md противоречит ERRATA — действует ERRATA; где ERRATA противоречит
// RATES_REVIEW — действует RATES_REVIEW (см. ERRATA §C0, он это прямо предписывает).
// Пересчитанные валидационные эталоны: docs/REFERENCE_TARGETS.md.
// Все величины в СИ. Комментарии по-русски, идентификаторы английские.
//
// Соглашения, важные для солвера:
//   * ВСЕ функции от приведённого поля берут МОДУЛЬ аргумента (ERRATA E1):
//     x = |E|/(N*1e-21) Td. Отрицательный полупериод не даёт ни NaN, ни обращения
//     к неверному концу таблицы.
//   * REACTIONS[j].rate(...) возвращает СЫРУЮ константу скорости в единицах
//     reaction.kUnits ('m3/s' или 'm6/s'). Скорость реакции считается как
//     rate * произведение плотностей ВСЕХ реагентов из reaction.reagents
//     (третье тело O2 уже входит в reagents явно). Никаких «спрятанных» N.
//   * SPECIES[i].muN = mu*N [1/(м*В*с)] — инвариант по плотности; mu = muN/N.
//     Это и есть требуемое ERRATA §D масштабирование подвижности как 1/N.

// ===========================================================================
// ===  ЭЛЕКТРОННЫЕ КОНСТАНТЫ СКОРОСТИ — ТАБЛИЦЫ, А НЕ ФИТЫ                ===
// ===========================================================================
//
// Источник (единственный): docs/RATES_REVIEW.md — собственное решение
// двухчленного уравнения Больцмана (BOLOS, схема Hagelaar–Pitchford 2005) на
// сечениях LXCat/Phelps для чистого O2, T_gas = 300 K, стационарная ФРЭЭ.
// Валидировано против собственной опубликованной таблицы A.V. Phelps
// (JILA Information Center Report No. 28, 1985; Lawton & Phelps 1978) —
// согласие 1–4 % по транспорту и всем ключевым каналам.
// Разброс между наборами Phelps / IST-Lisbon / Biagi: ±5–7 %.
//
// Все табличные величины интерполируются ЛОГ-ЛОГ (ошибка < 1 % выше 25 Td,
// < 4 % при 8–25 Td) — RATES_REVIEW §7. Аналитические фиты сознательно НЕ
// используются: они дают ±30–50 % (RATES_REVIEW §7.1).
//
// ЧТО БЫЛО ЗАМЕНЕНО (ERRATA §C0/§D2/§D3):
//   E4 kDiss     аррениус по Te 4.23e-15*exp(-5.56/Te)  занижен в 7–10 раз
//                (а исходный фит PHYSICS.md — в 500…1e6 раз)
//   E5 kExcO2a   4.0e-16*exp(-2.14/Te)                  занижен в 1.3–2.5 раза
//   E3 kAtt3     2.0e-42*exp(-EN/30)+5e-45              форма неверна (не экспонента)
//   E2 kAtt2     8.0e-18*exp(-160/EN)/(1+EN/3000)       занижен в 8…1000 раз
//   E1 alphaN    фит Morrow 1985                        занижен ровно вдвое при 200–500 Td
//   ELECTRON_TABLE                                      была таблица ВОЗДУХА, не O2
//
// ⚠ СОСТАВ СУММЫ «ДИССОЦИАЦИЯ» — ЭТО НЕ ФОРМАЛЬНОСТЬ (ERRATA §D4).
// Неопределённость ОПРЕДЕЛЕНИЯ (фактор 3.3) в ~30 раз больше неопределённости
// самих данных (±6 %). Здесь принято:
//     k_diss = k(6.0 эВ) + k(8.4 эВ)
//   6.0 эВ = континуум Шумана–Рунге -> O(3P) + O(3P)   [IST-Lisbon: явная метка]
//   8.4 эВ = отталкивательный канал -> O(3P) + O(1D)   [IST-Lisbon: явная метка]
// НЕ включены:
//   4.5 эВ  полосы Герцберга (A3Su+, C3Du, c1Su-) — IST-Lisbon трактует их как
//           СВЯЗАННОЕ состояние (порог 4.5 эВ ниже предела диссоциации 5.12 эВ,
//           есть обратный сверхупругий процесс). Включение дало бы +21 % @100 Td.
//   9.97 эВ (<0.1 % до 300 Td), 14.7 эВ, диссоциативное прилипание (это E2).
// ⚠ Ловушка первоисточника: строка `SUM DISSOC.` в таблице самого Phelps — это
//   (прилипание + 8.4 + 9.97 + 14.7 эВ), она ИСКЛЮЧАЕТ канал 6.0 эВ и даёт
//   1.27e-15 @100 Td против нашей суммы 1.937e-15. Не сравнивать «в лоб».
// Каждый акт даёт 2 атома O — стехиометрия E4 это учитывает.
//
// Выше ~60 Td доминирует канал 8.4 эВ (64 % суммы при 100 Td, 78 % при 200 Td);
// именно он ПОЛНОСТЬЮ отсутствовал в аррениусовском фите — отсюда ×1000.
export const TUNING = {
  kDissScale: 1.0,          // множитель E4 (дефолт 1.0; неопределённость данных ±7 %,
                            //   неопределённость ОПРЕДЕЛЕНИЯ — фактор 3.3, см. выше)
  kO2aScale: 1.0,           // множитель E5 (дефолт 1.0; неопределённость ±5 %)
  kAtt3Scale: 1.0,          // множитель E3 (дефолт 1.0; ЧЕСТНАЯ неопределённость ±40 %:
                            //   расчёт по сечению Phelps даёт 1.76e-42 в тепловом
                            //   пределе против измеренных CPB-1962 2.8e-42)
  d6MassConvention: 'neutral', // 'neutral' (дефолт, ERRATA §D6-РАЗБОР) | 'reduced'
};

/** Снимок дефолтов TUNING — чтобы тесты и раннер могли восстановить состояние. */
export const TUNING_DEFAULTS = Object.freeze({ ...TUNING });

/** Частичная установка параметров тюнинга (валидирует имена ключей). */
export function setTuning(patch) {
  for (const k of Object.keys(patch)) {
    if (!(k in TUNING_DEFAULTS)) throw new Error(`setTuning: неизвестный ключ ${k}`);
    TUNING[k] = patch[k];
  }
  return TUNING;
}

/** Возврат к дефолтам. */
export function resetTuning() {
  Object.assign(TUNING, TUNING_DEFAULTS);
  return TUNING;
}
// ===========================================================================
// ===                КОНЕЦ БЛОКА ЭЛЕКТРОННЫХ КОНСТАНТ                     ===
// ===========================================================================

// ---------------------------------------------------------------------------
// Фундаментальные константы
// ---------------------------------------------------------------------------
export const QE = 1.602177e-19;        // элементарный заряд, Кл
export const KB = 1.380649e-23;        // постоянная Больцмана, Дж/К
export const EPS0 = 8.8541878128e-12;  // электрическая постоянная, Ф/м
export const ME = 9.10938e-31;         // масса электрона, кг
export const N_LOSCHMIDT = 2.687e25;   // число Лошмидта, м^-3 (нормировка K0)
export const TORR = 133.32236842;      // Па в 1 Торр
export const TD = 1e-21;               // 1 Td = 1e-21 В*м^2
export const M_O = 2.6567e-26;         // масса O,  кг
export const M_O2 = 5.3134e-26;        // масса O2, кг
export const M_O3 = 7.9701e-26;        // масса O3, кг
export const M_O4 = 1.06268e-25;       // масса O4, кг
export const EV_TO_K = 11604.5;        // 1 эВ в кельвинах

/** Плотность газа [м^-3] из давления [Торр] и температуры [К]: N = p/(kB*T). */
export function gasDensity(pressureTorr, tempK) {
  return (pressureTorr * TORR) / (KB * tempK);
}

/** Опорная плотность газа: 1 атм, 300 К. N = 2.446e25 м^-3. */
export const N_REF = gasDensity(760, 300);

/** Приведённое поле [Td] из поля [В/м] и плотности. ОБЯЗАТЕЛЬНО модуль (ERRATA E1). */
export function reducedField(E, N) {
  return Math.abs(E) / (N * TD);
}

// ---------------------------------------------------------------------------
// Электронная таблица LFA для ЧИСТОГО O2 — RATES_REVIEW.md §7.4 (+ D*N).
// столбцы: E/N [Td], mu*N [1/(м*В*с)], w_e [м/с], <eps> [эВ], D*N [1/(м*с)]
//
// ⚠ ПРЕДЫДУЩАЯ ТАБЛИЦА БЫЛА НЕ ДЛЯ O2 (ERRATA §D3 п.1, RATES_REVIEW §5.2).
// Судя по значениям (<eps> ~3.1 эВ и w_e = 1.15e5 м/с при 100 Td), она взята
// для ВОЗДУХА, где колебательное возбуждение N2 сильно охлаждает ФРЭЭ.
// В чистом O2 при 100 Td: <eps> = 4.362 эВ (было 3.10, занижено на 29 %),
// w_e = 1.830e5 м/с (было 1.15e5, занижено на 37 %).
// Это критично: eta = (k2*N + k3*N^2)/w_e — ошибка в w_e напрямую двигает
// точку alpha = eta, а через неё и напряжение горения.
//
// mu*N и <eps> — дословно из RATES_REVIEW §7.4.
// w_e — производная величина, w = (mu*N)*(E/N)*1e-21 (согласована с колонкой 1).
// D*N — из того же BOLOS-расчёта. Проверка: D/mu при 100 Td = 3.411 эВ против
// 3.373 эВ в опубликованной таблице Phelps (Report 28) — согласие 1.1 %.
// ⚠ Прежняя колонка D*N была построена по приближению D/mu = (2/3)<eps>
// (2.07 эВ при 100 Td) на заниженной <eps> — она занижена в 2.6 раза.
// ---------------------------------------------------------------------------
export const ELECTRON_TABLE = [
  [0.1, 1.958e25, 1.958e3, 0.130, 2.490e24],
  [0.3, 1.541e25, 4.623e3, 0.153, 2.480e24],
  [1, 1.184e25, 1.184e4, 0.202, 2.474e24],
  [3, 7.769e24, 2.331e4, 0.363, 2.523e24],
  [5, 5.397e24, 2.699e4, 0.741, 2.835e24],
  [10, 3.327e24, 3.327e4, 1.699, 3.933e24],
  [20, 2.646e24, 5.292e4, 2.524, 4.937e24],
  [30, 2.418e24, 7.254e4, 2.920, 5.329e24],
  [50, 2.174e24, 1.087e5, 3.424, 5.724e24],
  [75, 1.979e24, 1.484e5, 3.909, 6.021e24],
  [100, 1.830e24, 1.830e5, 4.362, 6.240e24],
  [150, 1.596e24, 2.394e5, 5.286, 6.583e24],
  [200, 1.419e24, 2.837e5, 6.254, 6.870e24],
  [300, 1.173e24, 3.519e5, 8.232, 7.348e24],
  [500, 9.139e23, 4.570e5, 11.931, 8.133e24],
  [700, 7.835e23, 5.484e5, 15.231, 8.785e24],
  [1000, 6.736e23, 6.736e5, 19.731, 9.635e24],
];

// Предвычисленные логарифмы — лог-лог интерполяция без аллокаций в горячем цикле.
const NT = ELECTRON_TABLE.length;
const LX = new Float64Array(NT);
const LMU = new Float64Array(NT);
const LEPS = new Float64Array(NT);
const LDE = new Float64Array(NT);
for (let i = 0; i < NT; i++) {
  LX[i] = Math.log(ELECTRON_TABLE[i][0]);
  LMU[i] = Math.log(ELECTRON_TABLE[i][1]);
  LEPS[i] = Math.log(ELECTRON_TABLE[i][3]);
  LDE[i] = Math.log(ELECTRON_TABLE[i][4]);
}
const LX0 = LX[0];
const LXN = LX[NT - 1];
const EN_MIN = 1e-6; // нижний клип E/N, защита от log(0)

function tableIndex(lx) {
  let lo = 0, hi = NT - 2;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (LX[mid] <= lx) lo = mid; else hi = mid - 1;
  }
  return lo;
}

function interpLog(lx, LY) {
  if (lx <= LX0) return Math.exp(LY[0]);      // clamp на концах
  if (lx >= LXN) return Math.exp(LY[NT - 1]);
  const i = tableIndex(lx);
  const w = (lx - LX[i]) / (LX[i + 1] - LX[i]);
  return Math.exp(LY[i] + w * (LY[i + 1] - LY[i]));
}

/** Модуль + клип аргумента E/N (ERRATA E1). */
function absEN(EN) {
  const x = EN < 0 ? -EN : EN;
  return x > EN_MIN ? x : EN_MIN;
}

/**
 * Фабрика лог-лог интерполятора по произвольной таблице k(E/N) (RATES_REVIEW §7).
 * Логарифмы предвычисляются один раз — в горячем цикле только бинпоиск и exp.
 * На концах — clamp (значения крайних узлов), как и в interpLog выше.
 * @param {number[]} X узлы E/N [Td], строго возрастающие
 * @param {number[]} Y значения k, строго положительные
 * @returns {(x:number)=>number}
 */
function makeLogLog(X, Y) {
  const n = X.length;
  const lx = new Float64Array(n), ly = new Float64Array(n);
  for (let i = 0; i < n; i++) { lx[i] = Math.log(X[i]); ly[i] = Math.log(Y[i]); }
  const lx0 = lx[0], lxn = lx[n - 1], y0 = Y[0], yn = Y[n - 1];
  return function (x) {
    const l = Math.log(x);
    if (l <= lx0) return y0;
    if (l >= lxn) return yn;
    let lo = 0, hi = n - 2;
    while (lo < hi) { const m = (lo + hi + 1) >> 1; if (lx[m] <= l) lo = m; else hi = m - 1; }
    const w = (l - lx[lo]) / (lx[lo + 1] - lx[lo]);
    return Math.exp(ly[lo] + w * (ly[lo + 1] - ly[lo]));
  };
}

/** Приведённая подвижность электронов mu_e*N [1/(м*В*с)] от E/N [Td]. */
export function muE_N(EN) {
  return interpLog(Math.log(absEN(EN)), LMU);
}

/** Приведённый коэффициент диффузии электронов D_e*N [1/(м*с)] от E/N [Td]. */
export function DE_N(EN) {
  return interpLog(Math.log(absEN(EN)), LDE);
}

/** Средняя энергия электронов <eps> [эВ] от E/N [Td]. */
export function meanEnergy(EN) {
  return interpLog(Math.log(absEN(EN)), LEPS);
}

/** Электронная температура [К]: Te = (2/3)*<eps>*11604.5. */
export function electronTemp(EN) {
  return (2 / 3) * meanEnergy(EN) * EV_TO_K;
}

/** Электронная температура [эВ]. */
export function electronTempEv(EN) {
  return (2 / 3) * meanEnergy(EN);
}

/** Дрейфовая скорость электронов [м/с] (модуль): w_e = (mu*N)*(E/N)*1e-21. */
export function driftElectron(EN) {
  const x = absEN(EN);
  return muE_N(x) * x * TD;
}

// ---------------------------------------------------------------------------
// Таунсендовские коэффициенты
// ---------------------------------------------------------------------------

// --- E1: приведённый коэффициент ионизации alpha/N [м^2] --------------------
// alpha/N = k_ioniz(E/N) / w_e(E/N), k_ioniz — канал IONIZATION 12.06 эВ
// того же решения уравнения Больцмана (сечения Phelps), w_e — из ELECTRON_TABLE.
//
// ⚠ ЗАМЕНЁН ФИТ MORROW (Phys. Rev. A 32 (1985) 1799): ERRATA §D3 п.3.
// Проверено численно (RATES_REVIEW §7.2 п.6 просил перепроверить):
//   отношение расчёт/Morrow = 3.10 (100 Td), 2.66 (125), 2.02 (200), 2.00 (500).
// То есть «вдвое при 200–500 Td» подтверждается ТОЧНО, а ниже 150 Td фит
// занижен ещё сильнее. С фитом Morrow и правильным eta точка alpha = eta
// уезжает на ~168 Td, что противоречит и расчёту (128 Td), и эксперименту
// (пробойное поле O2 30–32 кВ/см = 122–130 Td).
//
// Независимая проверка ПО ОПУБЛИКОВАННОЙ таблице самого Phelps (Report 28,
// строки IONIZATION 12.06 / 2-BODY ATTACH. / 3-BODY ATTACH. / W при 40…250 Td):
// alpha = eta при 126.6 Td против 127.8 Td по этой таблице — расхождение 1 %.
// Сам k_ioniz расчёта на 5–12 % ниже таблицы Phelps (консервативно).
const EN_ION = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 85, 90, 95, 100,
                110, 115, 120, 125, 130, 135, 140, 150, 160, 175, 190, 200, 225, 250,
                275, 300, 350, 400, 450, 500, 600, 700, 800, 850, 1000];
const A_ION = [2.7035e-32, 1.4795e-29, 6.3164e-28, 7.6897e-27, 4.6061e-26, 1.7792e-25,
               5.1386e-25, 1.2117e-24, 2.4677e-24, 4.5013e-24, 7.5471e-24, 1.1852e-23,
               1.7638e-23, 3.4562e-23, 4.6119e-23, 5.9988e-23, 7.6376e-23, 1.1726e-22,
               1.4205e-22, 1.6990e-22, 2.0092e-22, 2.3519e-22, 2.7282e-22, 3.1383e-22,
               4.0632e-22, 5.1290e-22, 6.9942e-22, 9.1754e-22, 1.0800e-21, 1.5428e-21,
               2.0789e-21, 2.6780e-21, 3.3293e-21, 4.7484e-21, 6.2660e-21, 7.8298e-21,
               9.4031e-21, 1.2489e-20, 1.5419e-20, 1.8156e-20, 1.9450e-20, 2.3048e-20];
const ALPHA_TAB = makeLogLog(EN_ION, A_ION);

/** Приведённый коэффициент ионизации alpha/N [м^2] (таблица + лог-лог). */
export function alphaN(EN) {
  const x = absEN(EN);
  if (x <= 10) return 0;          // ниже 10 Td alpha/N < 3e-35 — численный шум
  return ALPHA_TAB(x);
}

/**
 * Приведённый коэффициент прилипания eta/N [м^2], самосогласованно из E2+E3:
 * eta = (k2*N + k3*N^2)/w_e  =>  eta/N = (k2 + k3*N)/w_e.
 */
export function etaN(EN, N = N_REF) {
  const x = absEN(EN);
  const we = driftElectron(x);
  if (we <= 0) return 0;
  return (kAtt2(x) + kAtt3(x) * N) / we;
}

/** Эффективный коэффициент размножения (alpha - eta)*N [1/м]. */
export function alphaEff(EN, N = N_REF) {
  return (alphaN(EN) - etaN(EN, N)) * N;
}

// ---------------------------------------------------------------------------
// Константы электронных реакций E1..E8
// ---------------------------------------------------------------------------

/** E1: e + O2 -> 2e + O2+ [м^3/с] = (alpha/N)*w_e. */
export function kIoniz(EN) {
  return alphaN(EN) * driftElectron(EN);
}

// --- E2: e + O2 -> O^- + O, диссоциативное прилипание [м^3/с] ---------------
// ⚠ ЗАМЕНЁН фит 8.0e-18*exp(-160/EN)/(1+EN/3000) — занижен в 8…1000 раз
// (×27 при 100 Td, ERRATA §D3 п.2). Он ломал eta, alphaEff и весь баланс O^-.
// Сверка с опубликованной таблицей Phelps (строка 2-BODY ATTACH.):
// 4.201e-17 против 4.18e-17 при 100 Td — 0.5 %.
const EN_DA = [3, 4, 5, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60, 65, 75,
               90, 100, 125, 135, 150, 175, 200, 225, 250, 300, 350, 400, 500,
               600, 700, 850, 1000];
const K_DA = [6.188e-25, 8.957e-23, 1.799e-21, 5.241e-20, 1.454e-19, 5.812e-19,
              1.405e-18, 3.258e-18, 7.270e-18, 1.156e-17, 1.562e-17, 1.926e-17,
              2.249e-17, 2.781e-17, 3.196e-17, 3.370e-17, 3.668e-17, 4.018e-17,
              4.201e-17, 4.533e-17, 4.626e-17, 4.733e-17, 4.843e-17, 4.892e-17,
              4.898e-17, 4.876e-17, 4.780e-17, 4.654e-17, 4.521e-17, 4.271e-17,
              4.063e-17, 3.897e-17, 3.711e-17, 3.582e-17];
const ATT2_TAB = makeLogLog(EN_DA, K_DA);

/** E2: e + O2 -> O^- + O [м^3/с] (таблица + лог-лог, RATES_REVIEW §7). */
export function kAtt2(EN) {
  const x = absEN(EN);
  if (x <= 3) return 0;           // ниже 3 Td k < 6e-25, порог канала 4.2 эВ
  return ATT2_TAB(x);
}

// --- E3: e + 2 O2 -> O2^- + O2, трёхтельное прилипание [м^6/с] --------------
// ⚠ ЗАМЕНЁН фит A*exp(-EN/EN0)+B. Экспоненциальная форма НЕВЕРНА: реальная
// зависимость выше 25 Td степенная (~ (E/N)^-0.9), из-за чего при 150–500 Td
// прежний фит занижал в 3–5 раз, а в тепловом пределе — в 9 раз.
// Расчёт по сечению Phelps, низкополевой участок ПРИВЯЗАН к прямому измерению
// Chanin, Phelps & Biondi, Phys. Rev. 128, 219 (1962):
//   2.800e-42 м^6/с при 300 K (тепловые электроны), максимум 5.000e-42 при
//   <eps> = 0.09 эВ. Поправочный множитель плавно спадает 1.59 -> 1.00 на
//   отрезке 0.01…5 Td, выше 5 Td значения чисто расчётные.
// НЕМОНОТОННА: максимум при E/N ~ 0.03 Td, затем спад — это физика, не баг.
// ⚠ ЧЕСТНАЯ НЕОПРЕДЕЛЁННОСТЬ ±40 % (ERRATA §D4): расчёт по сечению Phelps даёт
// в тепловом пределе 1.76e-42 против измеренных 2.8e-42, формула Kossyi 1992 —
// 1.9e-42. Значение по набору Biagi (1.5e-41) ОТВЕРГНУТО как завышенное на
// порядок относительно прямого измерения. Масштаб — TUNING.kAtt3Scale.
const EN_A3 = [0.01, 0.03, 0.1, 0.3, 0.5, 1, 2, 3, 4, 5, 7, 8, 10, 12, 15, 20,
               25, 30, 35, 40, 50, 65, 75, 90, 100, 135, 150, 200, 225, 300,
               350, 400, 500, 600, 700, 850, 1000];
const K_A3 = [2.800e-42, 5.000e-42, 4.550e-42, 3.886e-42, 3.549e-42, 2.935e-42,
              2.251e-42, 1.806e-42, 1.418e-42, 1.090e-42, 6.697e-43, 5.460e-43,
              3.911e-43, 3.067e-43, 2.406e-43, 1.893e-43, 1.637e-43, 1.480e-43,
              1.377e-43, 1.295e-43, 1.173e-43, 1.043e-43, 9.739e-44, 8.851e-44,
              8.347e-44, 6.926e-44, 6.431e-44, 5.149e-44, 4.653e-44, 3.550e-44,
              3.040e-44, 2.643e-44, 2.082e-44, 1.705e-44, 1.432e-44, 1.160e-44,
              9.688e-45];
const ATT3_TAB = makeLogLog(EN_A3, K_A3);

/** E3: e + 2 O2 -> O2^- + O2 [м^6/с] (таблица + лог-лог). */
export function kAtt3(EN) {
  return TUNING.kAtt3Scale * ATT3_TAB(absEN(EN));
}

// --- E4: e + O2 -> e + O + O, диссоциация электронным ударом [м^3/с] --------
// СОСТАВ СУММЫ: k(6.0 эВ) + k(8.4 эВ). Полосы Герцберга 4.5 эВ НЕ включены.
// Полное обоснование и предупреждения — в шапке файла (блок TUNING).
// Опорные точки: 4.729e-16 (50 Td), 1.937e-15 (100), 5.877e-15 (200),
// 1.245e-14 (400). Каждый акт даёт 2 атома O (стехиометрия E4).
const EN_DISS = [3, 4, 5, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60, 65, 75,
                 90, 100, 125, 135, 150, 175, 200, 225, 250, 300, 350, 400,
                 500, 600, 700, 850, 1000];
const K_DISS = [5.800e-27, 5.496e-24, 3.304e-22, 3.241e-20, 1.307e-19, 8.799e-19,
                3.028e-18, 1.018e-17, 3.455e-17, 7.463e-17, 1.294e-16, 1.976e-16,
                2.783e-16, 4.729e-16, 7.067e-16, 8.366e-16, 1.119e-15, 1.594e-15,
                1.937e-15, 2.866e-15, 3.257e-15, 3.856e-15, 4.870e-15, 5.877e-15,
                6.858e-15, 7.800e-15, 9.543e-15, 1.109e-14, 1.245e-14, 1.470e-14,
                1.648e-14, 1.792e-14, 1.961e-14, 2.093e-14];
const DISS_TAB = makeLogLog(EN_DISS, K_DISS);

/** E4: e + O2 -> e + O + O [м^3/с] (таблица + лог-лог). */
export function kDiss(EN) {
  const x = absEN(EN);
  if (x <= 3) return 0;           // ниже 3 Td k < 6e-27 — пренебрежимо
  return TUNING.kDissScale * DISS_TAB(x);
}

// --- E5: e + O2 -> e + O2(a1Dg), порог 0.977 эВ [м^3/с] ---------------------
// Опорные точки: 5.394e-16 (50 Td), 6.982e-16 (100), 8.649e-16 (200);
// насыщение 9.612e-16 при 400–600 Td, дальше слабый спад (ФРЭЭ уходит выше
// максимума сечения). Прежний аррениусовский фит занижал в 1.3–2.5 раза при
// 25–200 Td, а исходный фит PHYSICS.md — в ~50 раз.
const EN_A1 = [1, 2, 3, 4, 5, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60, 65,
               75, 90, 100, 125, 135, 150, 175, 200, 225, 250, 300, 350, 400,
               500, 600, 700, 850, 1000];
const K_A1 = [2.249e-23, 6.320e-20, 1.635e-18, 9.058e-18, 2.494e-17, 7.460e-17,
              1.027e-16, 1.575e-16, 2.056e-16, 2.634e-16, 3.342e-16, 3.863e-16,
              4.276e-16, 4.615e-16, 4.908e-16, 5.394e-16, 5.795e-16, 5.973e-16,
              6.299e-16, 6.728e-16, 6.982e-16, 7.532e-16, 7.722e-16, 7.981e-16,
              8.349e-16, 8.649e-16, 8.891e-16, 9.085e-16, 9.358e-16, 9.516e-16,
              9.597e-16, 9.612e-16, 9.524e-16, 9.387e-16, 9.149e-16, 8.908e-16];
const A1_TAB = makeLogLog(EN_A1, K_A1);

/** E5: e + O2 -> e + O2(a1Dg) [м^3/с] (таблица + лог-лог). */
export function kExcO2a(EN) {
  const x = absEN(EN);
  if (x <= 0.5) return 0;
  return TUNING.kO2aScale * A1_TAB(x);
}

/** E6: e + O2+ -> O + O [м^3/с], Te в кельвинах. Kossyi 1992. */
export function kRecO2p(Te) {
  return 2.0e-13 * Math.pow(300 / Te, 0.7);
}

/** E7: e + O4+ -> O2 + O2 [м^3/с], Te в кельвинах. Kossyi 1992. */
export function kRecO4p(Te) {
  return 1.4e-12 * Math.pow(300 / Te, 0.5);
}

/** E8: e + O3 -> O- + O2 [м^3/с]. Kossyi 1992. */
export const K_E8 = 1.0e-15;

// ---------------------------------------------------------------------------
// Ион-ионная рекомбинация (ERRATA §D: трёхтельный вклад ОБЯЗАТЕЛЕН)
// ---------------------------------------------------------------------------
/** k_ii = 2e-13 + 2e-37*N [м^3/с]; при N = 2.45e25 даёт 5.1e-12. */
export function kIonIon(N) {
  return 2.0e-13 + 2.0e-37 * N;
}

// ---------------------------------------------------------------------------
// Конверсия / кластеризация ионов C1..C6 (Kossyi 1992)
// ---------------------------------------------------------------------------
export function kC1(T) { return 2.4e-42 * Math.pow(300 / T, 3.2); }  // [м^6/с]
export function kC2(T) { return 3.3e-12 * Math.pow(300 / T, 4) * Math.exp(-5030 / T); } // [м^3/с]
export const K_C3 = 1.1e-42;  // O- + 2O2 -> O3- + O2 [м^6/с]
export const K_C4 = 5.3e-16;  // O- + O3 -> O3- + O
export const K_C5 = 4.0e-16;  // O2- + O3 -> O3- + O2
export const K_C6 = 2.5e-16;  // O3- + O -> O2- + O2

// ---------------------------------------------------------------------------
// Отлипание D1..D7 (Kossyi 1992)
// ---------------------------------------------------------------------------
export const K_D1 = 3.0e-16;  // O- + O2(a) -> O3 + e
export const K_D2 = 5.0e-16;  // O- + O -> O2 + e
export const K_D3 = 2.0e-16;  // O2- + O2(a) -> 2 O2 + e
export const K_D4 = 1.5e-16;  // O2- + O -> O3 + e
export const K_D5 = 5.0e-21;  // O- + O2 -> O3 + e
export const K_D7 = 1.0e-17;  // O3- + O -> 2 O2 + e

/**
 * Эффективная температура иона по Ванье для D6.
 *
 * ERRATA §D6-РАЗБОР: два независимых ревью разошлись в ×18–26.
 *   'neutral' (ДЕФОЛТ): T_eff = T_g + m_g*v_d^2/(3*kB), полная масса нейтрала.
 *      Обоснование методологическое: фит k = 2.7e-16*sqrt(Teff/300)*exp(-5590/Teff)
 *      (Kossyi 1992) откалиброван именно в этой конвенции; подставлять в чужой
 *      предэкспонент чужое определение T_eff нельзя.
 *   'reduced': mu_r = m_g/2 для пары O2- / O2 (приведённая масса).
 * Обе ветви обязаны существовать — прогон чувствительности идёт по обеим.
 *
 * @param {number} vd   дрейфовая скорость иона [м/с]
 * @param {number} T    температура газа [К]
 * @param {string} conv 'neutral' | 'reduced'
 */
export function wannierTemp(vd, T, conv = TUNING.d6MassConvention) {
  const m = conv === 'reduced' ? M_O2 / 2 : M_O2;
  return T + (m * vd * vd) / (3 * KB);
}

/**
 * D6: O2- + O2 -> 2 O2 + e, полевое (столкновительное) отлипание [м^3/с].
 * @param {number} E     модуль поля [В/м] (модуль берётся внутри)
 * @param {number} muO2m подвижность O2- при текущей плотности [м^2/(В*с)]
 * @param {number} T     температура газа [К]
 * @param {string} conv  конвенция массы
 */
export function kD6(E, muO2m, T, conv = TUNING.d6MassConvention) {
  const vd = muO2m * Math.abs(E);
  const Teff = wannierTemp(vd, T, conv);
  return 2.7e-16 * Math.sqrt(Teff / 300) * Math.exp(-5590 / Teff);
}

/** Частота отлипания D6 [1/с] при плотности N: nu = k*N. */
export function nuD6(E, muO2m, T, N, conv = TUNING.d6MassConvention) {
  return kD6(E, muO2m, T, conv) * N;
}

// ---------------------------------------------------------------------------
// Нейтральная химия N1..N7
// ---------------------------------------------------------------------------
export function kN1(T) { return 6.0e-46 * Math.pow(300 / T, 2.4); } // O + O2 + M -> O3 + M [м^6/с]
export function kN2(T) { return 3.0e-46 * (300 / T); }              // O + O + M -> O2 + M [м^6/с]
export function kN3(T) { return 8.0e-18 * Math.exp(-2060 / T); }    // O + O3 -> 2 O2
export function kN4(T) { return 5.2e-17 * Math.exp(-2840 / T); }    // O3 + O2(a) -> O + 2 O2
export const K_N5 = 2.2e-24;                                        // O2(a) + O2 -> 2 O2
export const K_N6 = 7.0e-22;                                        // O2(a) + O -> O2 + O
export function kN7(T) { return 7.3e-16 * Math.exp(-11400 / T); }   // O3 + O2 -> O + 2 O2

// ---------------------------------------------------------------------------
// SPECIES — полный набор (ERRATA §D: без O4+/O3-/O2(a) модель неверна)
//   muN = mu*N [1/(м*В*с)]  => mu(N) = muN/N   (масштабирование как 1/N)
//   K0  — приведённая подвижность при N0 = 2.687e25 (Ellis 1976/1978)
//   mu1atm — контрольное значение при 1 атм/300 К, ERRATA §D
// ---------------------------------------------------------------------------
function ionMuN(K0) { return K0 * N_LOSCHMIDT; }

export const SPECIES = [
  {
    id: 'e', name: 'electron', formula: 'e', mass: ME, z: -1, charged: true,
    fieldDependent: true, muN: null, DN: null, n0: 1e13,
    note: 'LFA: mu_e(E/N) из таблицы, D_e = (D_e*N)/N',
  },
  { id: 'O2p', name: 'molecular oxygen ion', formula: 'O2+', mass: M_O2, z: +1, charged: true,
    K0: 2.20e-4, muN: ionMuN(2.20e-4), mu1atm: 2.41e-4, n0: 1e13 },
  { id: 'O4p', name: 'oxygen cluster ion', formula: 'O4+', mass: M_O4, z: +1, charged: true,
    K0: 2.16e-4, muN: ionMuN(2.16e-4), mu1atm: 2.37e-4, n0: 0 },
  { id: 'Om', name: 'atomic oxygen anion', formula: 'O-', mass: M_O, z: -1, charged: true,
    K0: 3.20e-4, muN: ionMuN(3.20e-4), mu1atm: 3.51e-4, n0: 0 },
  { id: 'O2m', name: 'molecular oxygen anion', formula: 'O2-', mass: M_O2, z: -1, charged: true,
    K0: 2.16e-4, muN: ionMuN(2.16e-4), mu1atm: 2.37e-4, n0: 0 },
  { id: 'O3m', name: 'ozonide ion', formula: 'O3-', mass: M_O3, z: -1, charged: true,
    K0: 2.40e-4, muN: ionMuN(2.40e-4), mu1atm: 2.63e-4, n0: 0 },
  { id: 'O', name: 'atomic oxygen', formula: 'O', mass: M_O, z: 0, charged: false,
    muN: 0, D300: 2.0e-5, gammaWall: 2e-3, n0: 1e15 },
  { id: 'O3', name: 'ozone', formula: 'O3', mass: M_O3, z: 0, charged: false,
    muN: 0, D300: 1.5e-5, gammaWall: 1e-5, n0: 0 },
  { id: 'O2a', name: 'singlet delta oxygen', formula: 'O2(a1Dg)', mass: M_O2, z: 0, charged: false,
    muN: 0, D300: 2.0e-5, gammaWall: 3e-4, n0: 0 },
  { id: 'O2', name: 'background oxygen', formula: 'O2', mass: M_O2, z: 0, charged: false,
    background: true, muN: 0, n0: N_REF },
];

export const SPECIES_BY_ID = Object.fromEntries(SPECIES.map((s) => [s.id, s]));

/** Порядок транспортируемых сортов (без фонового O2) — контракт state.n. */
export const SPECIES_IDS = ['e', 'O2p', 'O4p', 'Om', 'O2m', 'O3m', 'O', 'O3', 'O2a'];

/** Заряд сорта в единицах e (для проверок сохранения заряда). */
export function speciesCharge(id) {
  const s = SPECIES_BY_ID[id];
  if (!s) throw new Error(`speciesCharge: неизвестный сорт ${id}`);
  return s.z;
}

/** Подвижность иона [м^2/(В*с)] при плотности N: mu = muN/N. */
export function ionMobility(id, N) {
  const s = SPECIES_BY_ID[id];
  if (!s || s.muN == null) throw new Error(`ionMobility: нет muN для ${id}`);
  return s.muN / N;
}

/** Диффузия иона по Эйнштейну при T: D = mu*kB*T/e. */
export function ionDiffusion(id, N, T) {
  return (ionMobility(id, N) * KB * T) / QE;
}

/** Диффузия нейтрала: D(T,p) = D300 * (T/300)^1.75 * (N_REF/N). */
export function neutralDiffusion(id, N, T) {
  const s = SPECIES_BY_ID[id];
  if (!s || s.D300 == null) throw new Error(`neutralDiffusion: нет D300 для ${id}`);
  return s.D300 * Math.pow(T / 300, 1.75) * (N_REF / N);
}

/** Тепловая скорость (средняя максвелловская) [м/с]. */
export function vThermal(massKg, tempK) {
  return Math.sqrt((8 * KB * tempK) / (Math.PI * massKg));
}

/** Тепловая скорость электронов из <eps> LFA [м/с]. */
export function vThermalElectron(EN) {
  return Math.sqrt((8 * QE * electronTempEv(EN)) / (Math.PI * ME));
}

// ---------------------------------------------------------------------------
// REACTIONS — декларативный список.
//
//   rate(EN, N, T) -> СЫРАЯ константа скорости в единицах kUnits.
//   Скорость = rate * произведение n[s] по ВСЕМ элементам reagents (с кратностью).
//   Фон O2 входит в reagents явно, его плотность = N.
//   kRange — ожидаемый диапазон значений k, используется тестом размерностей F5.
//   dimNote — пояснение, если диапазон выходит за родовое «унитарное окно».
// ---------------------------------------------------------------------------
export const REACTIONS = [
  // --- электронный удар ---
  { id: 'E1', name: 'ionization', equation: 'e + O2 -> 2e + O2+', type: 'electron_impact',
    threshold: 12.06, reagents: ['e', 'O2'], products: ['e', 'e', 'O2p'],
    kUnits: 'm3/s', kRange: [0, 1e-13], rate: (EN) => kIoniz(EN),
    source: 'Boltzmann/Phelps (RATES_REVIEW §7.2 п.6); фит Morrow 1985 отвергнут (×2..×3 занижен)' },

  { id: 'E2', name: 'dissociative attachment', equation: 'e + O2 -> O- + O', type: 'electron_impact',
    threshold: 4.2, reagents: ['e', 'O2'], products: ['Om', 'O'],
    kUnits: 'm3/s', kRange: [0, 1e-13], rate: (EN) => kAtt2(EN),
    source: 'таблица Boltzmann/Phelps, RATES_REVIEW §7 (прежний фит занижен в 8..1000 раз)' },

  { id: 'E3', name: '3-body attachment', equation: 'e + O2 + O2 -> O2- + O2', type: 'electron_impact',
    reagents: ['e', 'O2', 'O2'], products: ['O2m', 'O2'],
    kUnits: 'm6/s', kRange: [5e-45, 1e-41], rate: (EN) => kAtt3(EN),
    dimNote: 'немонотонна: максимум 5.0e-42 при ~0.03 Td (CPB 1962), спад ~ (E/N)^-0.9; ±40 %',
    source: 'Chanin-Phelps-Biondi 1962 (тепловой якорь) + Boltzmann/Phelps, RATES_REVIEW §7' },

  { id: 'E4', name: 'dissociation', equation: 'e + O2 -> e + O + O', type: 'electron_impact',
    threshold: 6.0, reagents: ['e', 'O2'], products: ['e', 'O', 'O'],
    kUnits: 'm3/s', kRange: [0, 1e-13], rate: (EN) => kDiss(EN),
    source: 'таблица Boltzmann/Phelps, RATES_REVIEW §7; СУММА = каналы 6.0 эВ + 8.4 эВ (4.5 эВ НЕ включён)' },

  { id: 'E5', name: 'singlet excitation', equation: 'e + O2 -> e + O2(a)', type: 'electron_impact',
    threshold: 0.98, reagents: ['e', 'O2'], products: ['e', 'O2a'],
    kUnits: 'm3/s', kRange: [0, 1e-13], rate: (EN) => kExcO2a(EN),
    source: 'таблица Boltzmann/Phelps, RATES_REVIEW §7 (насыщение 9.6e-16 при 400-600 Td)' },

  { id: 'E6', name: 'dissociative recombination', equation: 'e + O2+ -> O + O', type: 'recombination',
    reagents: ['e', 'O2p'], products: ['O', 'O'],
    kUnits: 'm3/s', kRange: [1e-15, 1e-12],
    dimNote: 'диссоциативная рекомбинация быстрее родового окна: 2e-13 м^3/с = 2e-7 см^3/с',
    rate: (EN) => kRecO2p(electronTemp(EN)), source: 'Kossyi 1992' },

  { id: 'E7', name: 'cluster recombination', equation: 'e + O4+ -> O2 + O2', type: 'recombination',
    reagents: ['e', 'O4p'], products: ['O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-14, 1e-11],
    dimNote: 'кластерная рекомбинация 1.4e-12 м^3/с = 1.4e-6 см^3/с (Kossyi), выше родового окна',
    rate: (EN) => kRecO4p(electronTemp(EN)), source: 'Kossyi 1992' },

  { id: 'E8', name: 'ozone attachment', equation: 'e + O3 -> O- + O2', type: 'electron_impact',
    reagents: ['e', 'O3'], products: ['Om', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_E8, source: 'Kossyi 1992' },

  // --- конверсия ионов ---
  { id: 'C1', name: 'cluster formation', equation: 'O2+ + O2 + O2 -> O4+ + O2', type: 'conversion',
    reagents: ['O2p', 'O2', 'O2'], products: ['O4p', 'O2'],
    kUnits: 'm6/s', kRange: [1e-46, 1e-40], rate: (EN, N, T) => kC1(T), source: 'Kossyi 1992' },

  { id: 'C2', name: 'cluster breakup', equation: 'O4+ + O2 -> O2+ + O2 + O2', type: 'conversion',
    reagents: ['O4p', 'O2'], products: ['O2p', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: (EN, N, T) => kC2(T), source: 'Kossyi 1992' },

  { id: 'C3', name: 'ozonide formation', equation: 'O- + O2 + O2 -> O3- + O2', type: 'conversion',
    reagents: ['Om', 'O2', 'O2'], products: ['O3m', 'O2'],
    kUnits: 'm6/s', kRange: [1e-46, 1e-40], rate: () => K_C3, source: 'Kossyi 1992' },

  { id: 'C4', name: 'charge transfer', equation: 'O- + O3 -> O3- + O', type: 'conversion',
    reagents: ['Om', 'O3'], products: ['O3m', 'O'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_C4, source: 'Kossyi 1992' },

  { id: 'C5', name: 'charge transfer', equation: 'O2- + O3 -> O3- + O2', type: 'conversion',
    reagents: ['O2m', 'O3'], products: ['O3m', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_C5, source: 'Kossyi 1992' },

  { id: 'C6', name: 'charge transfer', equation: 'O3- + O -> O2- + O2', type: 'conversion',
    reagents: ['O3m', 'O'], products: ['O2m', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_C6, source: 'Kossyi 1992' },

  // --- отлипание ---
  { id: 'D1', name: 'detachment (singlet)', equation: 'O- + O2(a) -> O3 + e', type: 'detachment',
    reagents: ['Om', 'O2a'], products: ['O3', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D1, source: 'Kossyi 1992' },

  { id: 'D2', name: 'associative detachment', equation: 'O- + O -> O2 + e', type: 'detachment',
    reagents: ['Om', 'O'], products: ['O2', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D2, source: 'Kossyi 1992' },

  { id: 'D3', name: 'detachment (singlet)', equation: 'O2- + O2(a) -> O2 + O2 + e', type: 'detachment',
    reagents: ['O2m', 'O2a'], products: ['O2', 'O2', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D3, source: 'Kossyi 1992' },

  { id: 'D4', name: 'associative detachment', equation: 'O2- + O -> O3 + e', type: 'detachment',
    reagents: ['O2m', 'O'], products: ['O3', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D4, source: 'Kossyi 1992' },

  { id: 'D5', name: 'slow detachment', equation: 'O- + O2 -> O3 + e', type: 'detachment',
    reagents: ['Om', 'O2'], products: ['O3', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D5, source: 'Kossyi 1992' },

  { id: 'D6', name: 'field detachment', equation: 'O2- + O2 -> O2 + O2 + e', type: 'detachment',
    reagents: ['O2m', 'O2'], products: ['O2', 'O2', 'e'],
    kUnits: 'm3/s', kRange: [1e-30, 1e-14],
    dimNote: 'при E->0 экспонента exp(-5590/300) даёт 2.2e-24; диапазон широкий по построению',
    rate: (EN, N, T) => kD6(absEN(EN) * TD * N, ionMuN(2.16e-4) / N, T),
    source: 'Kossyi 1992 + Ванье, конвенция массы — TUNING.d6MassConvention' },

  { id: 'D7', name: 'associative detachment', equation: 'O3- + O -> O2 + O2 + e', type: 'detachment',
    reagents: ['O3m', 'O'], products: ['O2', 'O2', 'e'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: () => K_D7, source: 'Kossyi 1992' },

  // --- ион-ионная рекомбинация (k_ii с трёхтельным вкладом, ERRATA §D) ---
  { id: 'R1', name: 'ion-ion recombination', equation: 'O- + O2+ -> O + O2', type: 'ion_ion',
    reagents: ['Om', 'O2p'], products: ['O', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11],
    dimNote: 'k_ii = 2e-13 + 2e-37*N = 5.1e-12 при 1 атм — выше родового окна',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  { id: 'R2', name: 'ion-ion recombination', equation: 'O2- + O2+ -> O2 + O2', type: 'ion_ion',
    reagents: ['O2m', 'O2p'], products: ['O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11], dimNote: 'см. R1',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  { id: 'R3', name: 'ion-ion recombination', equation: 'O3- + O2+ -> O3 + O2', type: 'ion_ion',
    reagents: ['O3m', 'O2p'], products: ['O3', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11], dimNote: 'см. R1',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  { id: 'R4', name: 'ion-ion recombination', equation: 'O- + O4+ -> O + O2 + O2', type: 'ion_ion',
    reagents: ['Om', 'O4p'], products: ['O', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11], dimNote: 'см. R1',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  { id: 'R5', name: 'ion-ion recombination', equation: 'O2- + O4+ -> O2 + O2 + O2', type: 'ion_ion',
    reagents: ['O2m', 'O4p'], products: ['O2', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11], dimNote: 'см. R1',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  { id: 'R6', name: 'ion-ion recombination', equation: 'O3- + O4+ -> O3 + O2 + O2', type: 'ion_ion',
    reagents: ['O3m', 'O4p'], products: ['O3', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-11], dimNote: 'см. R1',
    rate: (EN, N) => kIonIon(N), source: 'Kossyi 1992' },

  // --- нейтральная химия ---
  { id: 'N1', name: 'ozone formation', equation: 'O + O2 + O2 -> O3 + O2', type: 'neutral_3body',
    reagents: ['O', 'O2', 'O2'], products: ['O3', 'O2'],
    kUnits: 'm6/s', kRange: [1e-46, 1e-40], rate: (EN, N, T) => kN1(T), source: 'JPL/Atkinson' },

  { id: 'N2', name: 'atom recombination', equation: 'O + O + O2 -> O2 + O2', type: 'neutral_3body',
    reagents: ['O', 'O', 'O2'], products: ['O2', 'O2'],
    kUnits: 'm6/s', kRange: [1e-46, 1e-40], rate: (EN, N, T) => kN2(T), source: 'Kossyi 1992' },

  { id: 'N3', name: 'ozone destruction', equation: 'O + O3 -> O2 + O2', type: 'neutral',
    reagents: ['O', 'O3'], products: ['O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: (EN, N, T) => kN3(T), source: 'JPL' },

  { id: 'N4', name: 'singlet-ozone', equation: 'O3 + O2(a) -> O + O2 + O2', type: 'neutral',
    reagents: ['O3', 'O2a'], products: ['O', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-22, 1e-13], rate: (EN, N, T) => kN4(T), source: 'Kossyi 1992' },

  { id: 'N5', name: 'singlet quenching', equation: 'O2(a) + O2 -> O2 + O2', type: 'quenching',
    reagents: ['O2a', 'O2'], products: ['O2', 'O2'],
    kUnits: 'm3/s', kRange: [1e-25, 1e-13],
    dimNote: 'тушение O2(a) очень медленное: 2.2e-24 м^3/с (tau ~ 19 мс) — ниже родового окна',
    rate: () => K_N5, source: 'Kossyi 1992' },

  { id: 'N6', name: 'singlet quenching', equation: 'O2(a) + O -> O2 + O', type: 'quenching',
    reagents: ['O2a', 'O'], products: ['O2', 'O'],
    kUnits: 'm3/s', kRange: [1e-25, 1e-13],
    dimNote: '7.0e-22 м^3/с — на нижней границе родового окна',
    rate: () => K_N6, source: 'Kossyi 1992' },

  { id: 'N7', name: 'thermal ozone decay', equation: 'O3 + O2 -> O + O2 + O2', type: 'neutral',
    reagents: ['O3', 'O2'], products: ['O', 'O2', 'O2'],
    kUnits: 'm3/s', kRange: [0, 1e-13],
    dimNote: 'при 300 К пренебрежимо мал (~1e-32), включён для полноты',
    rate: (EN, N, T) => kN7(T), source: 'Kossyi 1992' },
];

export const REACTIONS_BY_ID = Object.fromEntries(REACTIONS.map((r) => [r.id, r]));

/**
 * Стехиометрический вектор реакции: сколько частиц сорта s рождается (>0) /
 * гибнет (<0) за один акт. Фоновый O2 включён (солвер его игнорирует).
 */
export function stoichiometry(reaction) {
  const nu = Object.create(null);
  for (const s of reaction.reagents) nu[s] = (nu[s] || 0) - 1;
  for (const s of reaction.products) nu[s] = (nu[s] || 0) + 1;
  return nu;
}

/** Изменение суммарного заряда в реакции (должно быть строго 0). */
export function chargeBalance(reaction) {
  let q = 0;
  const nu = stoichiometry(reaction);
  for (const id of Object.keys(nu)) q += nu[id] * speciesCharge(id);
  return q;
}

// ---------------------------------------------------------------------------
// Критерии пробоя (валидация ERRATA §C V1/V2)
// ---------------------------------------------------------------------------

/**
 * Критерий самоподдержания для ЭЛЕКТРООТРИЦАТЕЛЬНОГО газа (ERRATA §C V2):
 *     f(x) = gamma * alpha/(alpha - eta) * [exp((alpha - eta)*d) - 1] - 1
 * Упрощённая форма (alpha-eta)d = ln(1+1/gamma) НЕ используется.
 */
export function selfSustainResidual(EN, d, gamma, N = N_REF) {
  const a = alphaN(EN) * N;
  const e = etaN(EN, N) * N;
  const ae = a - e;
  if (ae <= 0) return -1;               // размножения нет — критерий не выполнен
  const arg = Math.min(ae * d, 700);    // защита от переполнения exp
  return gamma * (a / ae) * Math.expm1(arg) - 1;
}

/** Корень критерия самоподдержания по E/N [Td] (бисекция). */
export function breakdownEN(d = 1e-3, gamma = 0.02, N = N_REF) {
  // Скобка начинается ниже критического поля alpha = eta (~128 Td при 1 атм):
  // там ae <= 0 и невязка тождественно -1, так что корень заведомо внутри.
  let lo = 30, hi = 1000;
  if (selfSustainResidual(lo, d, gamma, N) > 0) return lo;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (selfSustainResidual(mid, d, gamma, N) < 0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

/** Точка alpha = eta по E/N [Td] (ERRATA §C V1, «поле смены знака», не пробой). */
export function alphaEtaCrossEN(N = N_REF) {
  let lo = 30, hi = 300;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (alphaN(mid) - etaN(mid, N) < 0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// ---------------------------------------------------------------------------
// Численные утилиты, общие для солвера и тестов
// ---------------------------------------------------------------------------

/** Функция Бернулли B(x) = x/(exp(x)-1) — ядро схемы Шарфеттера–Гуммеля. */
export function bernoulli(x) {
  if (x > 40) return 0;
  if (x < -40) return -x;
  const ax = x < 0 ? -x : x;
  if (ax < 1e-10) return 1 - 0.5 * x;
  return x / Math.expm1(x);
}

/** phi1(z) = (1-exp(-z))/z — для экспоненциального интегратора источников. */
export function phi1(z) {
  const az = z < 0 ? -z : z;
  if (az < 1e-8) return 1 - 0.5 * z;
  return -Math.expm1(-z) / z;
}
