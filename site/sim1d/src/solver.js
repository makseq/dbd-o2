// solver.js — 1D fluid-солвер ДБР в чистом O2 (drift-diffusion + Пуассон + LFA).
// Схема: см. docs/NUMERICS.md
//   * конечно-объёмная неравномерная сетка по ВСЕМУ стеку (диэлектрик | газ | диэлектрик),
//     диэлектрики — ячейки с eps_r и нулевым транспортом;
//   * потоки Шарфеттера–Гуммеля, электроны — неявно (трёхдиагональная M-матрица),
//     ионы — явно с ограничением CFL;
//   * полунеявный (заряд-сохраняющий) Пуассон: eps_eff = eps + dt*sigma_cond,
//     затем коррекция потоков Gamma = Gamma* + z*mu*n_up*(E^{n+1} - E^n);
//   * поверхностный заряд на границах газ/диэлектрик через исключение потенциала поверхности;
//   * источники — экспоненциальный интегратор (строгая положительность);
//   * ток по Сато–Морроу с весовым (лапласовским) полем.

import {
  QE, KB, EPS0, TD,
  gasDensity, muE_N, DE_N, meanEnergy,
  kIoniz, kAtt2, kAtt3, kDiss, kExc, kRecO2p, kRecO4p, K_E8,
  kIonIon, kC1, kC2, K_C3, K_C4, K_C5, K_C6,
  K_D1, K_D2, K_D3, K_D4, K_D5, K_D7, kD6,
  kN1, kN2, kN3, kN4, K_N5, K_N6,
  SPECIES_BY_ID, vThermal, vThermalElectron,
  bernoulli, sgFlux, phi1, thomas,
} from './physics.js';

// Индексы заряженных сортов в горячем цикле
const S_E = 0, S_O2P = 1, S_O4P = 2, S_OM = 3, S_O2M = 4, S_O3M = 5;
const NCH = 6;
const CH_IDS = ['e', 'O2p', 'O4p', 'Om', 'O2m', 'O3m'];
const CH_Z = [-1, 1, 1, -1, -1, -1];

const ND = 10;              // ячеек в каждом диэлектрике (профиль там линеен => хватает)
// ERRATA A2: пол плотности — ИСКЛЮЧИТЕЛЬНО защита от round-off, один и тот же для
// всех сортов. Значение 1e10 (ранее применявшееся к e/O2+/O4+) эквивалентно скрытому
// источнику ~5e17 м^-3 с^-1 и смещало напряжение горения на ~6%. Физическая затравка
// задаётся ТОЛЬКО членом sBg и каналами отлипания D1..D7.
const N_FLOOR_MIN = 1.0;    // абсолютный пол для всех сортов, м^-3
const QV_CAP = 4096;        // отсчётов Q–V на период
const QV_GHOSTS = 8;        // сколько ЗАВЕРШЁННЫХ периодов Q–V хранить для UI
const QV_GHOST_N = 512;     // отсчётов в сохранённой огибающей периода
const HIST_CAP = 4096;
const RT_N = 2048;          // узлов в таблицах коэффициентов по E/N

const DEFAULTS = {
  gapMM: 1.0,
  dielMM1: 0.5,
  dielMM2: 0.5,
  epsR: 9.0,
  areaCM2: 1.0,
  U0kV: 10.0,
  freqKHz: 10.0,
  pressureTorr: 760.0,
  tempK: 300.0,
  gamma: 0.02,
  seedDensity: 1e13,
  nCells: null,        // null => берётся из пресета режима (demo 200 / default 300 / accurate 800)
  ballastOhm: 0.0,
  mode: 'default',
  // дополнительные (не входят в обязательный контракт, но нужны спеке)
  sBg: 1e13,          // фоновый источник пар e + O2+, м^-3 с^-1
  betaGrid: 2.4,      // сила tanh-сгущения к поверхностям
  chemistry: true,    // выключается в тестах транспорта/сохранения заряда
  surfLeakTau: 0,     // 0 = идеальный диэлектрик; иначе постоянная утечки sigma, с
  // --- пристеночное ГУ (см. _wallCoeffs) ---
  wallBC: 'hagelaar', // 'hagelaar' (по умолчанию) | 'legacy' (старая форма, для сравнения)
  // reflE — коэффициент отражения ЭЛЕКТРОНОВ от диэлектрика, физический диапазон 0..0.2.
  // ВЫБОР МОДЕЛИ, НЕ ПОДГОНКА: r=0 (идеальное поглощение) — идеализация, при которой
  // в пристеночной ячейке остаётся непогасающая «плазма одной ячейки» и атомарный
  // кислород выходит в 140 раз плотнее самого газа ([O]/N = 1.4e2 — заведомый артефакт).
  // При r >= 0.05 артефакт исчезает ([O]/N = 1.0e-2) и восстанавливается нормальная
  // сеточная сходимость озона: 200 ячеек 5208 ppm -> 300 561 ppm -> 400(β=3.2) 309 ppm,
  // последнее уже в эталонном диапазоне 40..400 ppm (1e21..1e22 м^-3).
  // Устойчивость от reflE НЕ зависит вообще; ёмкости не меняются в 4-м знаке.
  // Развёртка: test/wall-bc-sweep.mjs, разбор: docs/DIVERGENCE_ANALYSIS.md §9.
  reflE: 0.05,
  reflI: 0.0,         // коэффициент отражения ИОНОВ от диэлектрика
};

const MODE_PRESETS = {
  // eCFL — множитель дрейфового CFL электронов. Транспорт электронов НЕЯВНЫЙ,
  // поэтому это ограничение чисто по точности (размытие фронта лавины), а не по
  // устойчивости; именно оно определяет скорость счёта (см. NUMERICS §5.1).
  // demo поднят с 200 до 300 ячеек: на 200 химия непригодна (озон 5208 ppm против
  // эталонных 40..400), на 300 — 561 ppm при том же порядке скорости счёта.
  demo: { nCells: 300, cfl: 0.8, dETol: 0.05, eCFL: 8 },
  default: { nCells: 300, cfl: 0.4, dETol: 0.03, eCFL: 4 },
  accurate: { nCells: 800, cfl: 0.2, dETol: 0.015, eCFL: 1.5 },
};

const now = (typeof performance !== 'undefined' && performance.now)
  ? () => performance.now()
  : () => Date.now();

export class DBDSolver {
  constructor(params = {}) {
    this.params = Object.assign({}, DEFAULTS, params);
    this.targetSimTime = Infinity;   // advance() останавливается по достижении
    this._buildAll();
  }

  // -------------------------------------------------------------------------
  // Параметры и сетка
  // -------------------------------------------------------------------------

  setParams(obj) {
    const geomKeys = ['gapMM', 'dielMM1', 'dielMM2', 'nCells', 'mode', 'betaGrid'];
    let rebuild = false, epsChanged = false;
    for (const k of Object.keys(obj)) {
      if (geomKeys.includes(k) && obj[k] !== this.params[k]) rebuild = true;
      if (k === 'epsR' && obj[k] !== this.params[k]) epsChanged = true;
      this.params[k] = obj[k];
    }
    if (rebuild) this._buildAll();
    else {
      this._deriveConstants();
      // eps_r входит не только в диагностические ёмкости, но и в МАТРИЦУ Пуассона
      // (epsCell) и в весовое поле Сато–Морроу (ELw). Без этого слайдер ε_r менял
      // только подписи, а решалась задача со старым ε_r (ошибка U_gap до 35%).
      if (epsChanged) this._refreshEps();
    }
    return this.params;
  }

  /** Перезалить eps по ячейкам и пересобрать весовое поле (горячая смена ε_r). */
  _refreshEps() {
    if (!this.epsCell) return;
    const p = this.params;
    for (let i = 0; i < this.Nt; i++) {
      this.epsCell[i] = this.gasMask[i] ? EPS0 : EPS0 * p.epsR;
    }
    this._buildWeightField();
  }

  _buildAll() {
    this._deriveConstants();
    this._buildGrid();
    this._allocate();
    this.reset();
  }

  _deriveConstants() {
    const p = this.params;
    // Конфигурация БЕЗ одного барьера схемой не поддерживается: при d=0 граница
    // газ/диэлектрик совпадает с металлом, матрица Пуассона получает вклад sigFace,
    // а _fieldFromPhi трактует эту грань как чистый металл — заряд копится на
    // металлическом электроде и решение расходится (U_gap ~ 1e9 В) молча, без NaN.
    // Лучше явный отказ, чем тихий мусор.
    if (!(p.dielMM1 > 0) || !(p.dielMM2 > 0)) {
      throw new Error('DBDSolver: dielMM1 и dielMM2 должны быть > 0 ' +
        '(односторонний ДБР схемой не поддерживается; используйте тонкий барьер, напр. 0.01 мм)');
    }
    if (!(p.gapMM > 0)) throw new Error('DBDSolver: gapMM должен быть > 0');
    const preset = MODE_PRESETS[p.mode] || MODE_PRESETS.default;
    this.cfl = preset.cfl;
    this.dETol = preset.dETol;
    this.eCFL = preset.eCFL;
    if (p.nCells == null) p.nCells = preset.nCells;

    this.Lg = p.gapMM * 1e-3;
    this.d1 = p.dielMM1 * 1e-3;
    this.d2 = p.dielMM2 * 1e-3;
    this.A = p.areaCM2 * 1e-4;
    this.U0 = p.U0kV * 1e3;
    this.freq = p.freqKHz * 1e3;
    this.omega = 2 * Math.PI * this.freq;
    this.period = 1 / this.freq;
    this.T = p.tempK;
    this.N = gasDensity(p.pressureTorr, p.tempK);

    // Ёмкости (аналитика, спека transport.circuit)
    this.Cd = (this.d1 + this.d2) > 0
      ? (EPS0 * p.epsR * this.A) / (this.d1 + this.d2)
      : Infinity;
    this.Cg = (EPS0 * this.A) / this.Lg;
    this.dEff = this.Lg + (this.d1 + this.d2) / p.epsR;   // эффективная толщина весового поля
    this.Ccell = (EPS0 * this.A) / this.dEff;

    // Подвижности ионов при текущих N, T (Ellis: mu = K0*N0/N)
    this.muIon = new Float64Array(NCH);
    this.DIon = new Float64Array(NCH);
    this.vthIon = new Float64Array(NCH);
    for (let s = 1; s < NCH; s++) {
      const sp = SPECIES_BY_ID[CH_IDS[s]];
      this.muIon[s] = sp.muN / this.N;
      this.DIon[s] = (this.muIon[s] * KB * this.T) / QE;
      this.vthIon[s] = vThermal(sp.mass, this.T);
    }
    // Нейтралы: тепловые скорости для пристеночной гибели
    this.vthO = vThermal(SPECIES_BY_ID.O.mass, this.T);
    this.vthO3 = vThermal(SPECIES_BY_ID.O3.mass, this.T);
    this.vthO2a = vThermal(SPECIES_BY_ID.O2a.mass, this.T);

    // Константы, не зависящие от поля
    this.kII = kIonIon(this.N);
    this.kC1N2 = kC1(this.T) * this.N * this.N;      // 1/с
    this.kC2N = kC2(this.T) * this.N;                // 1/с
    this.kC3N2 = K_C3 * this.N * this.N;             // 1/с
    this.kD5N = K_D5 * this.N;                       // 1/с
    this.kN1N2 = kN1(this.T) * this.N * this.N;      // 1/с (на атом O)
    this.kN2N = kN2(this.T) * this.N;                // м^3/с
    this.kN3v = kN3(this.T);
    this.kN4v = kN4(this.T);
    this.kN5N = K_N5 * this.N;                       // 1/с

    this._buildRateTables();
  }

  /**
   * Предвычисленные таблицы коэффициентов по E/N (лог-равномерная сетка).
   * В LFA все они — функции ТОЛЬКО от E/N при фиксированных N и T, поэтому
   * табулирование точно эквивалентно прямому счёту, но убирает ~12 вызовов
   * exp/pow на ячейку на шаг (это ~60% времени горячего цикла).
   */
  _buildRateTables() {
    const M = RT_N;
    this.rtLogMin = Math.log(1e-4);
    this.rtInvStep = (M - 1) / (Math.log(5000) - this.rtLogMin);
    const a = () => new Float64Array(M);
    const kiz = a(), ka2 = a(), ka3 = a(), kds = a(), kex = a();
    const kr2 = a(), kr4 = a(), kd6 = a(), mue = a(), de = a(), eps = a();
    const N = this.N, T = this.T, muO2m = SPECIES_BY_ID.O2m.muN / N;
    for (let j = 0; j < M; j++) {
      const EN = Math.exp(this.rtLogMin + j / this.rtInvStep);
      const em = meanEnergy(EN);
      const Te = (2 / 3) * em * 11604.5;
      kiz[j] = kIoniz(EN);
      ka2[j] = kAtt2(EN) * N;          // эффективная 2-тельная, 1/с на электрон/N
      ka3[j] = kAtt3(EN) * N * N;
      kds[j] = kDiss(EN) * N;
      kex[j] = kExc(EN) * N;
      kr2[j] = kRecO2p(Te);
      kr4[j] = kRecO4p(Te);
      kd6[j] = kD6(EN * TD * N, muO2m, T) * N;
      mue[j] = muE_N(EN) / N;
      de[j] = DE_N(EN) / N;
      eps[j] = em;
    }
    this.rt = { kiz, ka2, ka3, kds, kex, kr2, kr4, kd6, mue, de, eps };
  }

  /** Линейная интерполяция таблицы по log(E/N) с clamp на концах. */
  _rtIdx(EN) {
    const u = (Math.log(EN > 1e-4 ? EN : 1e-4) - this.rtLogMin) * this.rtInvStep;
    const j = u < 0 ? 0 : (u > RT_N - 2 ? RT_N - 2 : u | 0);
    this._rtW = u - j;
    return j;
  }

  _buildGrid() {
    const p = this.params;
    const Ng = Math.max(20, p.nCells | 0);
    const nd1 = this.d1 > 0 ? ND : 0;
    const nd2 = this.d2 > 0 ? ND : 0;
    const Nt = nd1 + Ng + nd2;
    this.Ng = Ng; this.nd1 = nd1; this.nd2 = nd2; this.Nt = Nt;
    this.i0 = nd1;                 // первая газовая ячейка
    this.i1 = nd1 + Ng - 1;        // последняя газовая ячейка
    this.fL = nd1;                 // левая грань газ/диэлектрик
    this.fR = nd1 + Ng;            // правая грань газ/диэлектрик

    const xf = new Float64Array(Nt + 1);
    // левый барьер — равномерно
    for (let k = 0; k <= nd1; k++) xf[k] = (this.d1 * k) / (nd1 || 1);
    if (nd1 === 0) xf[0] = 0;
    // газ — симметричное tanh-сгущение к обеим поверхностям
    // betaGrid = 0 — это РАВНОМЕРНАЯ сетка; без отдельной ветки tanh(0) = 0 в
    // знаменателе давал NaN во ВСЕЙ геометрии молча (шаг падал до 1e-14, и прогон
    // «просто не шёл»). Порог 1e-6 — там, где tanh(beta*s/2)/tanh(beta/2) уже
    // неотличимо от s в double.
    const beta = p.betaGrid;
    const th = Math.tanh(beta / 2);
    const uni = !(Math.abs(beta) > 1e-6) || !(Math.abs(th) > 0);
    for (let k = 0; k <= Ng; k++) {
      const s = (2 * k) / Ng - 1;
      const u = uni ? s : Math.tanh((beta * s) / 2) / th;
      xf[nd1 + k] = this.d1 + 0.5 * this.Lg * (1 + u);
    }
    // правый барьер
    for (let k = 1; k <= nd2; k++) xf[nd1 + Ng + k] = this.d1 + this.Lg + (this.d2 * k) / (nd2 || 1);

    const x = new Float64Array(Nt);
    const dx = new Float64Array(Nt);
    const gasMask = new Uint8Array(Nt);
    const epsCell = new Float64Array(Nt);
    for (let i = 0; i < Nt; i++) {
      dx[i] = xf[i + 1] - xf[i];
      x[i] = 0.5 * (xf[i] + xf[i + 1]);
      const gas = i >= this.i0 && i <= this.i1;
      gasMask[i] = gas ? 1 : 0;
      epsCell[i] = gas ? EPS0 : EPS0 * p.epsR;
    }
    const h = new Float64Array(Nt + 1);   // расстояние между центрами (на гранях)
    h[0] = dx[0] / 2;
    for (let f = 1; f < Nt; f++) h[f] = x[f] - x[f - 1];
    h[Nt] = dx[Nt - 1] / 2;

    this.xf = xf; this.x = x; this.dx = dx; this.gasMask = gasMask; this.epsCell = epsCell; this.h = h;

    // «ширина», приписанная грани, для квадратуры интеграла Сато–Морроу по газу
    const hw = new Float64Array(Nt + 1);
    for (let f = this.fL; f <= this.fR; f++) {
      hw[f] = (f === this.fL) ? dx[this.i0] / 2
        : (f === this.fR) ? dx[this.i1] / 2
          : h[f];
    }
    this.hw = hw;

    // минимальная ширина ячейки, примыкающей к грани — для ПОЛОКАЛЬНОГО CFL
    // (глобальный min(dx)/max(v) занижает шаг в разы: max|E| и min(dx) в
    //  разных местах зазора не совпадают)
    const hmin = new Float64Array(Nt + 1);
    for (let f = 0; f <= Nt; f++) {
      const a = f > 0 ? dx[f - 1] : dx[0];
      const b = f < Nt ? dx[f] : dx[Nt - 1];
      hmin[f] = a < b ? a : b;
    }
    this.hmin = hmin;

    // Контроль гладкости сетки (спека §1.2: |dx_{i+1}/dx_i - 1| < 0.08)
    let maxJump = 0;
    for (let i = this.i0; i < this.i1; i++) {
      const r = dx[i + 1] / dx[i];
      maxJump = Math.max(maxJump, Math.abs(r - 1));
    }
    this.gridSmoothness = maxJump;
    this.dxMin = Math.min(...Array.from(dx.slice(this.i0, this.i1 + 1)));
    this.dxMax = Math.max(...Array.from(dx.slice(this.i0, this.i1 + 1)));
  }

  _allocate() {
    const Nt = this.Nt, Nf = Nt + 1;
    // пристеночные коэффициенты ГУ (см. _wallCoeffs): [сорт][0]=левая, [1]=правая
    this.kwL = new Float64Array(NCH);
    this.kwR = new Float64Array(NCH);
    this.linL = new Float64Array(NCH);
    this.linR = new Float64Array(NCH);
    this.n = new Array(NCH);
    this.nNew = new Array(NCH);
    for (let s = 0; s < NCH; s++) {
      this.n[s] = new Float64Array(Nt);
      this.nNew[s] = new Float64Array(Nt);
    }
    this.nO = new Float64Array(Nt);
    this.nO3 = new Float64Array(Nt);
    this.nO2a = new Float64Array(Nt);
    // аккумуляторы медленной химии
    this.accO = new Float64Array(Nt);
    this.accO3 = new Float64Array(Nt);
    this.accO2a = new Float64Array(Nt);

    this.phi = new Float64Array(Nt);
    this.phiA = new Float64Array(Nt);
    this.phiB = new Float64Array(Nt);
    this.phiL = new Float64Array(Nt);      // весовое (лапласовское) поле, U=1
    this.rho = new Float64Array(Nt);
    this.Ecell = new Float64Array(Nt);
    this.ENcell = new Float64Array(Nt);
    this.ionizRate = new Float64Array(Nt);

    this.Ef = new Float64Array(Nf);        // поле на гранях (в газе — со стороны газа)
    this.EfOld = new Float64Array(Nf);
    this.EfA = new Float64Array(Nf);
    this.EfB = new Float64Array(Nf);
    this.ELw = new Float64Array(Nf);       // весовое поле Сато–Морроу [1/м]
    this.Gam = new Array(NCH);
    this.sigS = new Array(NCH);
    for (let s = 0; s < NCH; s++) {
      this.Gam[s] = new Float64Array(Nf);
      this.sigS[s] = new Float64Array(Nf);
    }
    this.sigCond = new Float64Array(Nf);
    this.Jstar = new Float64Array(Nf);
    this.sigFace = new Float64Array(Nf);   // поверхностный заряд, привязанный к грани
    this.muEf = new Float64Array(Nf);
    this.DEf = new Float64Array(Nf);

    this.Ma = new Float64Array(Nt);
    this.Mb = new Float64Array(Nt);
    this.Mc = new Float64Array(Nt);
    this.Mr = new Float64Array(Nt);
    this.cw = new Float64Array(Nt);
    this.dw = new Float64Array(Nt);
    this.Gface = new Float64Array(Nf);
    this.wLf = new Float64Array(Nf);
    this.wRf = new Float64Array(Nf);
    this.gLm = new Float64Array(Nf);   // модифицированные (полунеявные) полупроводимости
    this.gRm = new Float64Array(Nf);
    this._SigTmp = new Float64Array(Nf);

    // Теневые копии для ОТБРАКОВКИ шага (NUMERICS §5.2)
    this.shN = new Array(NCH);
    for (let s = 0; s < NCH; s++) this.shN[s] = new Float64Array(Nt);
    this.shO = new Float64Array(Nt);
    this.shO3 = new Float64Array(Nt);
    this.shO2a = new Float64Array(Nt);
    this.shAccO = new Float64Array(Nt);
    this.shAccO3 = new Float64Array(Nt);
    this.shAccO2a = new Float64Array(Nt);
    this.shEf = new Float64Array(Nf);
    this.shSigFace = new Float64Array(Nf);

    // Q–V выборка за период
    this.qvU = new Float64Array(QV_CAP);
    this.qvQ = new Float64Array(QV_CAP);
    this.qvI = new Float64Array(QV_CAP);
    this.qvUg = new Float64Array(QV_CAP);
    // огибающие Q–V ЗАВЕРШЁННЫХ периодов (UI рисует «призраки» из них, а не из
    // history: кольцо history вмещает всего ~2 периода)
    this.qvHistory = {
      cap: QV_GHOSTS, n: 0, head: 0, len: QV_GHOST_N,
      u: [], q: [], t0: new Float64Array(QV_GHOSTS),
    };
    for (let k = 0; k < QV_GHOSTS; k++) {
      this.qvHistory.u.push(new Float64Array(QV_GHOST_N));
      this.qvHistory.q.push(new Float64Array(QV_GHOST_N));
    }

    this.history = {
      t: new Float64Array(HIST_CAP),
      Uapp: new Float64Array(HIST_CAP),
      Ugap: new Float64Array(HIST_CAP),
      current: new Float64Array(HIST_CAP),
      // расширение контракта: чистый кондукционный (разрядный) ток. Без него UI
      // вынужден вычитать C_cell*dU/dt численно и ловит шум дифференцирования.
      condCurrent: new Float64Array(HIST_CAP),
      charge: new Float64Array(HIST_CAP),
      len: 0, head: 0, capacity: HIST_CAP,
    };
  }

  // -------------------------------------------------------------------------
  reset() {
    const p = this.params;
    const Nt = this.Nt;
    for (let s = 0; s < NCH; s++) { this.n[s].fill(0); this.nNew[s].fill(0); }
    this.nO.fill(0); this.nO3.fill(0); this.nO2a.fill(0);
    this.accO.fill(0); this.accO3.fill(0); this.accO2a.fill(0);
    for (let i = this.i0; i <= this.i1; i++) {
      this.n[S_E][i] = p.seedDensity;
      this.n[S_O2P][i] = p.seedDensity;
      this.nO[i] = SPECIES_BY_ID.O.n0;
    }
    this.sigL = 0; this.sigR = 0;
    this.sigFace.fill(0);
    this.t = 0;
    this.dt = 1e-12;
    this.Uel = 0; this.UelPrev = 0;
    this.charge = 0; this.chargeC = 0;   // Кэхэн-компенсация
    this.current = 0; this.condCurrent = 0; this.dispCurrent = 0; this._Iprev = 0;
    this._eSignificant = false;
    this.power = 0;
    this.peakCurrent = 0; this.maxEN = 0;
    this.steps = 0; this.rejects = 0;
    this.rejectReason = { neg: 0, field: 0, poisson: 0 };
    this.tAcc = 0;
    this.qvLen = 0;
    this.qvHistory.n = 0; this.qvHistory.head = 0;
    this.periodIndex = 0; this.periodT0 = 0;
    this.pulseCount = 0; this.pulseOn = false;
    this.breakdownsPerPeriod = 0;
    this.clipCharge = 0;
    this.history.len = 0; this.history.head = 0;
    this.lastHistT = -1;
    this.periodStats = {
      energyPerPeriodJ: 0, powerW: 0,
      Cdiel: this.Cd, CdielRaw: this.Cd, CdielIdentity: this.Cd,
      CdielGeom: this.Cd, CcellGeom: this.Ccell, Ccell: this.Ccell, UburnkV: 0,
      qvArcsOn: 0, qvArcsOff: 0, qvR2on: 0, qvR2off: 0, qvOk: false, qvFallback: true,
    };
    this.phi.fill(0); this.Ef.fill(0); this.EfOld.fill(0);
    this._buildWeightField();
    this._poissonSolve(0, true);        // стартовое поле
    this.EfOld.set(this.Ef);
    this._updateDiagnostics();
  }

  /** Напряжение источника */
  Usrc(t) { return this.U0 * Math.sin(this.omega * t); }

  // -------------------------------------------------------------------------
  // Пуассон
  // -------------------------------------------------------------------------

  /**
   * Сборка трёхдиагональной матрицы.
   * dtImp > 0 => полунеявная поправка eps_eff = eps + dt*sigma_cond.
   * Возвращает коэффициенты в Ma/Mb/Mc, проводимости граней в Gface, веса wLf/wRf.
   */
  _assembleMatrix(dtImp) {
    const Nt = this.Nt, dx = this.dx, eps = this.epsCell, gas = this.gasMask;
    const G = this.Gface, wL = this.wLf, wR = this.wRf, sc = this.sigCond;
    for (let f = 0; f <= Nt; f++) {
      if (f === 0) {
        G[f] = (2 * eps[0]) / dx[0];           // металл: бесконечная проводимость слева
        wL[f] = 0; wR[f] = 1;
      } else if (f === Nt) {
        G[f] = (2 * eps[Nt - 1]) / dx[Nt - 1];
        wL[f] = 1; wR[f] = 0;
      } else {
        const iL = f - 1, iR = f;
        const addL = gas[iL] ? dtImp * sc[f] : 0;
        const addR = gas[iR] ? dtImp * sc[f] : 0;
        const gL = (2 * (eps[iL] + addL)) / dx[iL];
        const gR = (2 * (eps[iR] + addR)) / dx[iR];
        G[f] = (gL * gR) / (gL + gR);
        wL[f] = gL / (gL + gR);
        wR[f] = gR / (gL + gR);
        this.gLm[f] = gL; this.gRm[f] = gR;
      }
    }
    const a = this.Ma, b = this.Mb, c = this.Mc;
    for (let i = 0; i < Nt; i++) {
      a[i] = i === 0 ? 0 : -G[i];
      c[i] = i === Nt - 1 ? 0 : -G[i + 1];
      b[i] = G[i] + G[i + 1];
    }
  }

  /**
   * Правая часть: заряд, поверхностный заряд, полунеявный поток, напряжение электрода.
   *
   * Полунеявный член F_f = dt*(J*_f - sigma_f*E^n_f) входит по-разному:
   *  - на внутренних газовых гранях это обычная дивергенция: rhs_i += F_left - F_right;
   *  - на границах газ/диэлектрик он входит В УРАВНЕНИЕ ПОВЕРХНОСТНОГО ЗАРЯДА
   *    (sigma^{n+1} = sigma^n -+ dt*J), поэтому складывается с sigma и получает
   *    тот же вес w, что и sigma. Знак противоположен на левой и правой поверхности,
   *    т.к. dSigma_L/dt = -J[fL], а dSigma_R/dt = +J[fR].
   */
  _assembleRhs(Uel, useCharge, dtImp) {
    const Nt = this.Nt, dx = this.dx, r = this.Mr;
    const G = this.Gface, wL = this.wLf, wR = this.wRf, sf = this.sigFace;
    const sc = this.sigCond, Js = this.Jstar, Eo = this.EfOld;
    const fL = this.fL, fR = this.fR;
    const Sig = this._SigTmp;
    if (useCharge) {
      Sig.set(sf);
      if (dtImp > 0) {
        Sig[fL] = sf[fL] - dtImp * (Js[fL] - sc[fL] * Eo[fL]);
        Sig[fR] = sf[fR] + dtImp * (Js[fR] - sc[fR] * Eo[fR]);
      }
    }
    for (let i = 0; i < Nt; i++) {
      let v = 0;
      const fl = i, fr = i + 1;
      if (useCharge) {
        v = this.rho[i] * dx[i] + wR[fl] * Sig[fl] + wL[fr] * Sig[fr];
        if (dtImp > 0 && this.gasMask[i]) {
          // дивергентный вклад потока — по всем граням газовой ячейки, включая
          // пристеночные (поток на стенку и уносит заряд из ячейки, и заряжает
          // поверхность — это два РАЗНЫХ вклада, оба обязательны)
          v += dtImp * (Js[fl] - sc[fl] * Eo[fl]);
          v -= dtImp * (Js[fr] - sc[fr] * Eo[fr]);
        }
      }
      r[i] = v;
    }
    r[0] += G[0] * Uel;
    // правый металл заземлён => вклад 0
  }

  /** Полный расчёт поля для заданного U (без полунеявности) — старт и тесты. */
  _poissonSolve(Uel, useCharge) {
    this._computeRho();
    this.sigCond.fill(0);
    this._assembleMatrix(0);
    this._assembleRhs(Uel, useCharge, 0);
    thomas(this.Ma, this.Mb, this.Mc, this.Mr, this.phi, this.Nt, this.cw, this.dw);
    this.Uel = Uel;
    this._fieldFromPhi(this.phi, this._SigTmp, this.Ef, Uel);
  }

  /** Весовое (лапласовское) поле: rho=0, sigma=0, U=1. Даёт E_L для Сато–Морроу. */
  _buildWeightField() {
    const sc = this.sigCond; sc.fill(0);
    this._assembleMatrix(0);
    this._assembleRhs(1, false, 0);
    thomas(this.Ma, this.Mb, this.Mc, this.Mr, this.phiL, this.Nt, this.cw, this.dw);
    this._fieldFromPhi(this.phiL, null, this.ELw, 1);
  }

  /**
   * Поле на гранях из потенциала.
   * На гранях газ/диэлектрик возвращается поле СО СТОРОНЫ ГАЗА (через потенциал поверхности).
   * useSigma=false — для линейных откликов (phiB, phiL).
   */
  _fieldFromPhi(phi, sigArr, out, Uleft) {
    const Nt = this.Nt, dx = this.dx, h = this.h;
    for (let f = 0; f <= Nt; f++) {
      if (f === 0) {
        out[f] = -(phi[0] - Uleft) / (dx[0] / 2);
      } else if (f === Nt) {
        out[f] = -(0 - phi[Nt - 1]) / (dx[Nt - 1] / 2);
      } else {
        const iL = f - 1, iR = f;
        if (this.gasMask[iL] === this.gasMask[iR]) {
          out[f] = -(phi[iR] - phi[iL]) / h[f];
        } else {
          // ВАЖНО: используем ТЕ ЖЕ модифицированные полупроводимости и тот же
          // эффективный поверхностный заряд, что и матрица — иначе поле у стенки
          // не согласовано с полунеявной схемой и Пуассон перестаёт выполняться.
          const gL = this.gLm[f], gR = this.gRm[f];
          const s = sigArr ? sigArr[f] : 0;
          const phis = (s + gL * phi[iL] + gR * phi[iR]) / (gL + gR);
          out[f] = this.gasMask[iR]
            ? -(phi[iR] - phis) / (dx[iR] / 2)
            : -(phis - phi[iL]) / (dx[iL] / 2);
        }
      }
    }
  }

  _computeRho() {
    const rho = this.rho;
    rho.fill(0);
    for (let s = 0; s < NCH; s++) {
      const ns = this.n[s], q = QE * CH_Z[s];
      for (let i = this.i0; i <= this.i1; i++) rho[i] += q * ns[i];
    }
  }

  // -------------------------------------------------------------------------
  // Один шаг
  // -------------------------------------------------------------------------

  /**
   * Один принятый шаг. Внутри — попытка + ОТБРАКОВКА (NUMERICS §5.2):
   * попытка откатывается, если (a) появились заметно отрицательные плотности,
   * (b) поле изменилось больше чем на 3*dE_tol, (c) невязка Пуассона > 1e-8.
   * После отката dt *= 0.4 и попытка повторяется (максимум 6 раз).
   */
  step() {
    let dtTry = this.dt;
    this._snapshot();
    let r = null;
    for (let attempt = 0; ; attempt++) {
      r = this._attempt(dtTry);
      if (r.ok || attempt >= 6 || r.dt <= 1.01e-14) break;
      this.rejects++;
      this.rejectReason[r.why] = (this.rejectReason[r.why] || 0) + 1;
      this._restore();
      dtTry = Math.max(1e-14, r.dt * 0.4);
    }
    const dt = r.dt;

    this.t = this._tNew;
    this.steps++;
    this._updateDiagnostics();
    this._record();
    this._periodBookkeeping();
    this.dt = this._nextDt(dt, r.maxIonV, r.maxWe, r.maxNu);
    return dt;
  }

  _snapshot() {
    for (let s = 0; s < NCH; s++) this.shN[s].set(this.n[s]);
    this.shO.set(this.nO); this.shO3.set(this.nO3); this.shO2a.set(this.nO2a);
    this.shAccO.set(this.accO); this.shAccO3.set(this.accO3); this.shAccO2a.set(this.accO2a);
    this.shEf.set(this.Ef); this.shSigFace.set(this.sigFace);
    this._sh = {
      sigL: this.sigL, sigR: this.sigR, charge: this.charge, chargeC: this.chargeC,
      Iprev: this._Iprev, Uel: this.Uel, UelPrev: this.UelPrev,
      clipCharge: this.clipCharge, tAcc: this.tAcc,
    };
  }

  _restore() {
    for (let s = 0; s < NCH; s++) this.n[s].set(this.shN[s]);
    this.nO.set(this.shO); this.nO3.set(this.shO3); this.nO2a.set(this.shO2a);
    this.accO.set(this.shAccO); this.accO3.set(this.shAccO3); this.accO2a.set(this.shAccO2a);
    this.Ef.set(this.shEf); this.sigFace.set(this.shSigFace);
    const s = this._sh;
    this.sigL = s.sigL; this.sigR = s.sigR; this.charge = s.charge; this.chargeC = s.chargeC;
    this._Iprev = s.Iprev; this.Uel = s.Uel; this.UelPrev = s.UelPrev;
    this.clipCharge = s.clipCharge; this.tAcc = s.tAcc;
  }

  /** Попытка шага. Возвращает {ok, why, dt, maxIonV, maxWe, maxNu}. */
  _attempt(dtIn) {
    const p = this.params;
    let dt = dtIn;
    const i0 = this.i0, i1 = this.i1, fL = this.fL, fR = this.fR;
    const dx = this.dx, h = this.h, N = this.N;

    // --- 1. Поле старого слоя, коэффициенты LFA на гранях -------------------
    this.EfOld.set(this.Ef);
    const Eo = this.EfOld;
    let maxIonV = 0, maxWe = 0, maxNuNet = 0;
    // rho СТАРОГО слоя — именно она входит в полунеявный Пуассон (§4.2),
    // перенос учитывается отдельно через div(J*).
    this._computeRho();

    const rt = this.rt, invNTD = 1 / (N * TD);
    let eCflLoc = Infinity, ionCflLoc = Infinity;
    for (let f = fL; f <= fR; f++) {
      const aE = Math.abs(Eo[f]);
      const j = this._rtIdx(aE * invNTD), w = this._rtW;
      const muE = rt.mue[j] + w * (rt.mue[j + 1] - rt.mue[j]);
      this.muEf[f] = muE;
      this.DEf[f] = rt.de[j] + w * (rt.de[j + 1] - rt.de[j]);
      const we = muE * aE;
      if (we > maxWe) maxWe = we;
      const hm = this.hmin[f];
      if (we > 0) { const c = hm / we; if (c < eCflLoc) eCflLoc = c; }
    }
    let maxMuIon = 0;
    for (let s = 1; s < NCH; s++) if (this.muIon[s] > maxMuIon) maxMuIon = this.muIon[s];

    // --- 2. Явные SG-потоки для ионов, неявный транспорт электронов ---------
    // 2a0. пристеночные коэффициенты ГУ (Хагелаар) — нужны и для линеаризации ниже
    this._wallCoeffs();

    // 2a. проводимость граней и апвинд-плотности (для полунеявного Пуассона)
    const sc = this.sigCond;
    sc.fill(0);
    for (let s = 0; s < NCH; s++) this.sigS[s].fill(0);

    for (let f = fL; f <= fR; f++) {
      const E = Eo[f];
      for (let s = 0; s < NCH; s++) {
        const z = CH_Z[s];
        const mu = s === S_E ? this.muEf[f] : this.muIon[s];
        let nup;
        if (f === fL) {
          // Пристеночная грань: линеаризация ГУ Хагелаара по полю.
          // d(Gamma·n^)/dE = lin * (v_d·n^)/E = -lin*z*mu (слева n^ = -z), а код
          // ниже применяет поправку G[f] += z*mu*nup*(E^{n+1}-E^n) с G[fL] = -Gamma·n^,
          // поэтому nup = lin*n. При обрезанном нулём потоке lin = 0 (поток не
          // зависит от поля), при открытом lin = 1/(1+r) — В ТОМ ЧИСЛЕ когда дрейф
          // направлен ОТ стенки: там поток тоже полевой, только с обратным знаком.
          nup = this.linL[s] * this.n[s][i0];
        } else if (f === fR) {
          nup = this.linR[s] * this.n[s][i1];
        } else {
          nup = (z * E >= 0) ? this.n[s][f - 1] : this.n[s][f];
        }
        const g = QE * mu * nup;
        this.sigS[s][f] = g;
        sc[f] += g;
      }
      if (f > fL && f < fR) {
        const v = maxMuIon * Math.abs(E);
        if (v > maxIonV) maxIonV = v;
        if (v > 0) { const c = this.hmin[f] / v; if (c < ionCflLoc) ionCflLoc = c; }
      }
    }

    // 2b. ограничение шага ДО транспорта (стабильность явной части — ионы)
    this._ionCfl = ionCflLoc; this._eCflT = eCflLoc;
    dt = this._limitDt(dt, maxIonV, maxWe);

    // 2c. неявный транспорт электронов (SG, M-матрица) при поле E^n
    this._electronTransportImplicit(dt);

    // 2d. явные SG-потоки ионов и пристеночные потоки
    this._ionFluxes(dt);

    // 2e. суммарный ток на гранях J* (для полунеявной правой части)
    const Js = this.Jstar;
    Js.fill(0);
    for (let s = 0; s < NCH; s++) {
      const G = this.Gam[s], q = QE * CH_Z[s];
      for (let f = fL; f <= fR; f++) Js[f] += q * G[f];
    }

    // --- 3. Полунеявный Пуассон + цепь -------------------------------------
    this._assembleMatrix(dt);
    this._assembleRhs(0, true, dt);
    thomas(this.Ma, this.Mb, this.Mc, this.Mr, this.phiA, this.Nt, this.cw, this.dw);
    this._assembleRhs(1, false, dt);
    thomas(this.Ma, this.Mb, this.Mc, this.Mr, this.phiB, this.Nt, this.cw, this.dw);

    this._fieldFromPhi(this.phiA, this._SigTmp, this.EfA, 0);
    this._fieldFromPhi(this.phiB, null, this.EfB, 1);

    const tNew = this.t + dt;
    const Usrc = this.Usrc(tNew);
    let Uel;
    // ток линеен по U_el: I = IA + Uel*kB  (Сато–Морроу + смещение)
    let IA = 0, kB = 0;
    for (let f = fL; f <= fR; f++) {
      const w = this.ELw[f] * this.hw[f];
      IA += (Js[f] + sc[f] * (this.EfA[f] - Eo[f])) * w;
      kB += sc[f] * this.EfB[f] * w;
    }
    IA *= this.A; kB *= this.A;
    const IAtot = IA - (this.Ccell * this.Uel) / dt;
    const kBtot = kB + this.Ccell / dt;
    if (p.ballastOhm > 0) {
      Uel = (Usrc - p.ballastOhm * IAtot) / (1 + p.ballastOhm * kBtot);
    } else {
      Uel = Usrc;
    }
    this.UelPrev = this.Uel;
    this.Uel = Uel;

    const phi = this.phi;
    for (let i = 0; i < this.Nt; i++) phi[i] = this.phiA[i] + Uel * this.phiB[i];
    // по линейности: поле — та же суперпозиция (согласовано с полунеявной матрицей)
    for (let f = 0; f <= this.Nt; f++) this.Ef[f] = this.EfA[f] + Uel * this.EfB[f];

    // --- 4. Коррекция потоков по новому полю (заряд-сохраняющая) -----------
    const Ef = this.Ef;
    for (let s = 0; s < NCH; s++) {
      const G = this.Gam[s], sg = this.sigS[s], z = CH_Z[s];
      for (let f = fL; f <= fR; f++) {
        G[f] += (z * sg[f] / QE) * (Ef[f] - Eo[f]);
      }
    }

    // --- 5. Консервативное обновление плотностей ---------------------------
    for (let s = 0; s < NCH; s++) {
      const G = this.Gam[s], ns = this.n[s];
      if (s === S_E) {
        // электроны уже перенесены неявно при E^n; применяем только поправку потока
        const sg = this.sigS[s], z = CH_Z[s];
        for (let i = i0; i <= i1; i++) {
          const dGl = (z * sg[i] / QE) * (Ef[i] - Eo[i]);
          const dGr = (z * sg[i + 1] / QE) * (Ef[i + 1] - Eo[i + 1]);
          ns[i] -= (dt * (dGr - dGl)) / dx[i];
        }
      } else {
        for (let i = i0; i <= i1; i++) {
          ns[i] -= (dt * (G[i + 1] - G[i])) / dx[i];
        }
      }
    }

    // --- 6. Поверхностный заряд -------------------------------------------
    let JL = 0, JR = 0;
    for (let s = 0; s < NCH; s++) {
      const q = QE * CH_Z[s];
      JL += q * this.Gam[s][fL];
      JR += q * this.Gam[s][fR];
    }
    this.sigL += -JL * dt;
    this.sigR += JR * dt;
    if (p.surfLeakTau > 0) {
      const dec = Math.exp(-dt / p.surfLeakTau);
      this.sigL *= dec; this.sigR *= dec;
    }
    this.sigFace[fL] = this.sigL;
    this.sigFace[fR] = this.sigR;

    // --- 7. Источники (reaction extent, строго сохраняющие заряд) -----------
    if (p.chemistry) maxNuNet = this._sources(dt);

    // --- 8. Критерий отбраковки (a): заметно отрицательные плотности --------
    let negBad = false;
    for (let s = 0; s < NCH; s++) {
      const ns = this.n[s];
      let mn = Infinity, mx = 0;
      for (let i = i0; i <= i1; i++) { const v = ns[i]; if (v < mn) mn = v; if (v > mx) mx = v; }
      if (mn < -1e-6 * mx) { negBad = true; break; }
    }
    if (negBad) return { ok: false, why: 'neg', dt, maxIonV, maxWe, maxNu: maxNuNet };

    // --- 8b. Пол положительности (round-off) -------------------------------
    let clip = 0;
    for (let s = 0; s < NCH; s++) {
      const ns = this.n[s];
      for (let i = i0; i <= i1; i++) {
        if (!(ns[i] > N_FLOOR_MIN)) { clip += (N_FLOOR_MIN - ns[i]) * dx[i] * CH_Z[s]; ns[i] = N_FLOOR_MIN; }
      }
    }
    this.clipCharge += QE * clip;
    for (let i = i0; i <= i1; i++) {
      if (!(this.nO[i] > 0)) this.nO[i] = 0;
      if (!(this.nO3[i] > 0)) this.nO3[i] = 0;
      if (!(this.nO2a[i] > 0)) this.nO2a[i] = 0;
    }

    // --- 9. Ток, заряд, диагностика ---------------------------------------
    let Icond = 0;
    for (let f = fL; f <= fR; f++) {
      let J = 0;
      for (let s = 0; s < NCH; s++) J += QE * CH_Z[s] * this.Gam[s][f];
      Icond += J * this.ELw[f] * this.hw[f];
    }
    Icond *= this.A;
    const Idisp = (this.Ccell * (this.Uel - this.UelPrev)) / dt;
    this.condCurrent = Icond;
    this.dispCurrent = Idisp;
    this.current = Icond + Idisp;

    // накопление заряда (трапеции + компенсация Кэхэна)
    const dQ = 0.5 * (this.current + this._Iprev) * dt;
    this._Iprev = this.current;
    const y = dQ - this.chargeC;
    const tt = this.charge + y;
    this.chargeC = (tt - this.charge) - y;
    this.charge = tt;
    this._tNew = tNew;

    // --- 10. Критерии отбраковки (b) скачок поля, (c) невязка Пуассона ------
    let maxE = 1, maxDE = 0;
    for (let f = fL; f <= fR; f++) {
      const ae = Math.abs(this.Ef[f]);
      if (ae > maxE) maxE = ae;
      const d = Math.abs(this.Ef[f] - Eo[f]);
      if (d > maxDE) maxDE = d;
    }
    this._maxDErel = maxDE / maxE;
    if (this._maxDErel > 3 * this.dETol) {
      return { ok: false, why: 'field', dt, maxIonV, maxWe, maxNu: maxNuNet };
    }
    // источники теперь сохраняют заряд точно, поэтому rho^{n+1} обязана быть
    // согласована с phi^{n+1} машинно; рост невязки = признак ухода схемы
    this._resid = this.poissonResidual();
    if (!(this._resid < 1e-8)) {
      return { ok: false, why: 'poisson', dt, maxIonV, maxWe, maxNu: maxNuNet };
    }
    return { ok: true, why: '', dt, maxIonV, maxWe, maxNu: maxNuNet };
  }

  // -------------------------------------------------------------------------
  // Пристеночное ГУ
  // -------------------------------------------------------------------------

  /**
   * ГУ ХАГЕЛААРА для потока заряженных частиц на диэлектрик.
   *
   * Hagelaar G J M, de Hoog F J, Kroesen G M W, "Boundary conditions in fluid
   * models of gas discharges", Phys. Rev. E 62 (2000) 1452, ур. (18)–(21).
   *
   * Реализуемая формула (поток НАРУЖУ, т.е. вдоль внешней нормали n^ стенки):
   *
   *     Gamma_s·n^ = (1-r)/(1+r) * (1/4)*v_th,s*n_s
   *                + (2a-1)/(1+r) * |mu_s*E_n|*n_s
   *                - (2/(1+r)) * D_s*(dn_s/dn)            <- НЕ реализован, см. ниже
   *     Gamma_s·n^ = max(0, ...)                          <- п.2 задания
   *
   *   r  — коэффициент отражения частиц сорта s от поверхности (params.reflE / reflI),
   *   a  = 1, если дрейф направлен К стенке (v_d·n^ > 0), и a = 0, если ОТ стенки;
   *   v_th,s = sqrt(8 k T_s / (pi m_s)) — средняя тепловая скорость (physics.vThermal,
   *            для электронов vThermalElectron(E/N) по LFA-температуре).
   *
   * Так как (2a-1)*|mu_s*E_n| — это ровно ЗНАКОВАЯ проекция дрейфовой скорости на
   * внешнюю нормаль, реализация записана в эквивалентной бесветвевой форме
   *
   *     kw_s = max(0,  (1-r)/(1+r) * (1/4)*v_th,s  +  (v_d·n^)/(1+r) ),
   *     Gamma_s·n^ = kw_s * n_s(пристеночная ячейка)                        [м^-2 с^-1]
   *
   * где v_d = z_s*mu_s*E (вектор вдоль +z), n^ = -z слева и +z справа.
   *
   * РАЗМЕРНОСТЬ: [v_th] = [mu*E] = м/с; kw*n = (м/с)(м^-3) = м^-2 с^-1 = поток. ОК.
   * Коэффициенты (1-r)/(1+r) и 1/(1+r) безразмерны, при r=0 равны 1.
   *
   * ЧТО ЭТО МЕНЯЕТ ПО СУЩЕСТВУ. Старое ГУ было
   *     kw = max(0, v_d·n^) + (1/4)*v_th
   * то есть при отталкивающем поле (v_d·n^ < 0) дрейфовые «ворота» просто
   * закрывались, и на стенку шёл ПОЛНЫЙ тепловой поток, НЕ ЗАВИСЯЩИЙ ОТ ПОЛЯ
   * (docs/DIVERGENCE_ANALYSIS.md §8.4: Gamma_e,dr == 0 на всём срыве). У Хагелаара
   * дрейфовый член при a=0 ВЫЧИТАЕТСЯ: поток на стенку падает вплоть до нуля при
   * |mu_e*E_n| >= (1-r)/(1+r)*(1/4)*v_th,e. Это ЛИНЕЙНОЕ потоко-ограничение, а не
   * больцмановский множитель exp(-e*Phi/kT) (тот пробовали, §8.4 — стало хуже).
   *
   * Член -2/(1+r)*D*dn/dn СОЗНАТЕЛЬНО НЕ реализован (в формуле выше он оставлен
   * для полноты ссылки). Причина: у Хагелаара он появляется потому, что там ГУ
   * замыкает дрейф-диффузионный поток В САМОЙ граничной точке, а здесь схема
   * конечно-объёмная и диффузия между пристеночной ячейкой и её соседом уже
   * посчитана внутренней гранью Шарфеттера–Гуммеля. Добавление градиента ещё и
   * на стенке дало бы двойной счёт диффузии в первой ячейке. Пристеночная ячейка
   * в этой схеме играет роль «граничной точки» Хагелаара, а её плотность n_s —
   * роль n на стенке.
   *
   * Заполняет kwL/kwR (>=0, м/с) и linL/linR — множитель при n для ЛИНЕАРИЗАЦИИ
   * пристеночного потока по полю в полунеявном Пуассоне (d Gamma/d E = lin*z*mu;
   * при закрытом (обрезанном нулём) потоке производная нулевая).
   */
  _wallCoeffs() {
    const fL = this.fL, fR = this.fR, Eo = this.EfOld, p = this.params;
    const legacy = p.wallBC === 'legacy';
    const clampR = (v) => (v > 0 ? (v < 0.95 ? v : 0.95) : 0);
    const rE = clampR(p.reflE), rI = clampR(p.reflI);
    const thE = legacy ? 0.25 : 0.25 * (1 - rE) / (1 + rE);
    const feE = legacy ? 1 : 1 / (1 + rE);
    const thI = legacy ? 0.25 : 0.25 * (1 - rI) / (1 + rI);
    const feI = legacy ? 1 : 1 / (1 + rI);

    for (let s = 0; s < NCH; s++) {
      const z = CH_Z[s];
      const isE = (s === S_E);
      const muL = isE ? this.muEf[fL] : this.muIon[s];
      const muR = isE ? this.muEf[fR] : this.muIon[s];
      const vthL = isE ? vThermalElectron(Math.abs(Eo[fL]) / (this.N * TD)) : this.vthIon[s];
      const vthR = isE ? vThermalElectron(Math.abs(Eo[fR]) / (this.N * TD)) : this.vthIon[s];
      const th = isE ? thE : thI, fe = isE ? feE : feI;
      // проекции дрейфовой скорости на ВНЕШНЮЮ нормаль: слева n^ = -z, справа n^ = +z
      const vnL = -(z * muL * Eo[fL]);
      const vnR = +(z * muR * Eo[fR]);
      let kL, kR;
      if (legacy) {
        kL = th * vthL + (vnL > 0 ? vnL : 0);
        kR = th * vthR + (vnR > 0 ? vnR : 0);
        this.linL[s] = vnL > 0 ? 1 : 0;
        this.linR[s] = vnR > 0 ? 1 : 0;
      } else {
        kL = th * vthL + fe * vnL;
        kR = th * vthR + fe * vnR;
        // п.2: поток не может быть направлен из стенки в газ (эмиссия — отдельные члены)
        if (!(kL > 0)) { kL = 0; this.linL[s] = 0; } else this.linL[s] = fe;
        if (!(kR > 0)) { kR = 0; this.linR[s] = 0; } else this.linR[s] = fe;
      }
      this.kwL[s] = kL;
      this.kwR[s] = kR;
    }
  }

  // -------------------------------------------------------------------------
  // Транспорт
  // -------------------------------------------------------------------------

  /** Неявный SG-транспорт электронов (снимает жёсткость дрейфа+диффузии электронов). */
  _electronTransportImplicit(dt) {
    const i0 = this.i0, i1 = this.i1, fL = this.fL, fR = this.fR;
    const dx = this.dx, h = this.h, ne = this.n[S_E];
    const a = this.Ma, b = this.Mb, c = this.Mc, r = this.Mr;
    const Eo = this.EfOld;
    const Ng = this.Ng;
    // alpha/beta коэффициентов SG на гранях: Gamma_f = alpha*n_left - beta*n_right
    const al = this._alFace || (this._alFace = new Float64Array(this.Nt + 1));
    const be = this._beFace || (this._beFace = new Float64Array(this.Nt + 1));
    for (let f = fL + 1; f < fR; f++) {
      const mu = this.muEf[f], D = this.DEf[f], E = Eo[f], hh = h[f];
      const X = (-1 * mu * E * hh) / D;         // z = -1
      const bx = bernoulli(X);
      al[f] = (D / hh) * (bx + X);
      be[f] = (D / hh) * bx;
    }
    // пристеночные коэффициенты — ГУ Хагелаара, посчитано в _wallCoeffs()
    const kwL = this.kwL[S_E], kwR = this.kwR[S_E];
    this._kwL = kwL; this._kwR = kwR;

    // вторичная эмиссия — по потокам положительных ионов на стенку (явно, по E^n)
    const semL = this.params.gamma * this._posWallFlux(fL, i0, -1);
    const semR = this.params.gamma * this._posWallFlux(fR, i1, +1);
    this._semL = semL; this._semR = semR;

    for (let k = 0; k < Ng; k++) {
      const i = i0 + k;
      const fl = i, fr = i + 1;
      const s = dt / dx[i];
      let bb = 1, aa = 0, ccx = 0, rr = ne[i];
      if (fl === fL) { bb += s * kwL; rr += s * semL; }
      else { aa = -s * al[fl]; bb += s * be[fl]; }
      if (fr === fR) { bb += s * kwR; rr += s * semR; }
      else { ccx = -s * be[fr]; bb += s * al[fr]; }
      a[k] = aa; b[k] = bb; c[k] = ccx; r[k] = rr;
    }
    const out = this._eTmp || (this._eTmp = new Float64Array(this.Nt));
    thomas(a, b, c, r, out, Ng, this.cw, this.dw);
    // потоки по НОВЫМ плотностям (согласованность: div(Gamma) = -(n^{n+1}-n^n)/dt)
    const G = this.Gam[S_E];
    G.fill(0);
    for (let k = 0; k < Ng; k++) ne[i0 + k] = out[k] > 0 ? out[k] : 0;
    for (let f = fL + 1; f < fR; f++) G[f] = al[f] * ne[f - 1] - be[f] * ne[f];
    G[fL] = -(kwL * ne[i0] - semL);
    G[fR] = kwR * ne[i1] - semR;
  }

  /**
   * Суммарный поток положительных ионов на стенку (модуль, >=0). dir=-1 левая, +1 правая.
   * Использует ТЕ ЖЕ коэффициенты Хагелаара, что и сам перенос ионов (_wallCoeffs):
   * на аноде, где поле отталкивает ионы, поток обнуляется, а вместе с ним и
   * вторичная эмиссия — иначе эмиссия шла бы с поверхности, на которую ионы не летят.
   */
  _posWallFlux(f, iw, dir) {
    const kw = dir < 0 ? this.kwL : this.kwR;
    let tot = 0;
    for (let s = 1; s < NCH; s++) {
      if (CH_Z[s] <= 0) continue;
      tot += kw[s] * this.n[s][iw];
    }
    return tot;
  }

  /**
   * Явные SG-потоки ионов + пристеночные потоки.
   *
   * Пристеночный отток ЯВНЫМ он быть не может: скорость на стенке — дрейф ПЛЮС
   * тепловая 0.25*v_th (~110–160 м/с), а пристеночная ячейка самая узкая (tanh-
   * сгущение), поэтому число Куранта dt*w/dx доходило до 2.2 (demo) и ~9 (accurate)
   * — плотность уходила в минус, обрезалась полом, а поверхность получала БОЛЬШЕ
   * заряда, чем в ячейке было. Здесь используется ТОЧНОЕ среднее по шагу решение
   * задачи опустошения ячейки  n' = -(w/dx) n:
   *
   *      w_eff = (dx/dt) * (1 - exp(-w*dt/dx))
   *
   * при малом Куранте w_eff -> w (явная схема), при большом — насыщается на dx/dt,
   * то есть за шаг уносится не более того, что в ячейке есть. Схема остаётся строго
   * консервативной (поверхность получает ровно то, что теряет ячейка) и положительной,
   * а шаг по времени ограничивать не нужно.
   */
  _ionFluxes(dt) {
    const fL = this.fL, fR = this.fR, i0 = this.i0, i1 = this.i1;
    const h = this.h, Eo = this.EfOld, dx = this.dx;
    const cL = dx[i0] / dt, cR = dx[i1] / dt;
    let cflL = 0, cflR = 0;
    for (let s = 1; s < NCH; s++) {
      const G = this.Gam[s], ns = this.n[s], mu = this.muIon[s], D = this.DIon[s], z = CH_Z[s];
      G.fill(0);
      for (let f = fL + 1; f < fR; f++) {
        G[f] = sgFlux(ns[f - 1], ns[f], mu, D, z, Eo[f], h[f]);
      }
      // левая стенка (коэффициенты ГУ Хагелаара из _wallCoeffs, >= 0)
      const wallL = this.kwL[s];
      if (wallL / cL > cflL) cflL = wallL / cL;
      G[fL] = -cL * (-Math.expm1(-wallL / cL)) * ns[i0];
      // правая стенка
      const wallR = this.kwR[s];
      if (wallR / cR > cflR) cflR = wallR / cR;
      G[fR] = cR * (-Math.expm1(-wallR / cR)) * ns[i1];
    }
    this.wallCfl = Math.max(cflL, cflR);   // диагностика (без ограничения шага)
  }

  // -------------------------------------------------------------------------
  // Источники
  // -------------------------------------------------------------------------

  /**
   * Источники по формулировке REACTION EXTENT (ERRATA A5) — строго сохраняют заряд.
   *
   * Прежняя схема (независимое экспоненциальное обновление каждого сорта) заряд НЕ
   * сохраняла: электрон получал ионизацию как точный экспоненциальный интеграл
   * n0*(e^z-1), а партнёрский O2+ — как трапецию k*eAvg*dt, что при z=1.4 больше на
   * 15%. Плюс sBg рождал электрон без иона вовсе. Итог — +1.7e-9 Кл фиктивного
   * заряда за период, оседавшего на диэлектриках (асимметрия sigL/sigR ~2:1).
   *
   * Схема здесь двухстадийная:
   *   стадия 1 — экспоненциальный предиктор даёт СРЕДНЮЮ ПО ШАГУ плотность
   *              nbar_s = (1/dt)∫n_s dt для каждого сорта (для линейного канала
   *              L*nbar*dt тождественно равно точному числу актов реакции);
   *   стадия 2 — для каждой реакции считается ОДИН неотрицательный ΔR_j
   *              (с ограничением по запасу реагентов) и применяется ОДНОЙ И ТОЙ ЖЕ
   *              величиной ко всем реагентам и продуктам по стехиометрии.
   * Поскольку каждая реакция здесь либо нейтральна по заряду, либо переносит его
   * между двумя сортами, сумма z_s*Δn_s равна нулю машинно точно.
   *
   * Возвращает max чистой скорости РОСТА (nu_iz − nu_att) для контроллера шага.
   */
  _sources(dt) {
    const i0 = this.i0, i1 = this.i1, N = this.N, T = this.T, p = this.params;
    const ne = this.n[S_E], nO2p = this.n[S_O2P], nO4p = this.n[S_O4P];
    const nOm = this.n[S_OM], nO2m = this.n[S_O2M], nO3m = this.n[S_O3M];
    const nO = this.nO, nO3 = this.nO3, nO2a = this.nO2a;
    const Ef = this.Ef;
    let maxNu = 0;
    const kII = this.kII;

    // Средняя по шагу плотность для dn/dt = P - L*n (точно):
    //   ∫n dt = n0*g + P*h2,  g = (1-e^{-z})/L,  h2 = (dt-g)/L,  z = L*dt.
    // Обе величины неотрицательны при любом знаке L, поэтому nbar >= 0 всегда.
    const nbarOf = (n0v, P, L) => {
      let zz = L * dt;
      if (zz > 700) zz = 700;
      if (zz < -40) zz = -40;             // защита от переполнения при взрывном росте
      const Lz = zz / dt;
      let g, h2;
      if (Math.abs(zz) > 1e-6) {
        const em = Math.expm1(-zz);
        g = -em / Lz;
        h2 = (dt - g) / Lz;
      } else {
        g = dt * (1 - 0.5 * zz + zz * zz / 6);
        h2 = dt * dt * (0.5 - zz / 6 + zz * zz / 24);
      }
      const r = (n0v * g + P * h2) / dt;
      return r > 0 ? r : 0;
    };
    // коэффициент ограничения расхода сорта: не больше, чем есть + приход
    const lim = (avail, cons) => (cons > avail && cons > 0 ? avail / cons : 1);

    for (let i = i0; i <= i1; i++) {
      const E = 0.5 * (Ef[i] + Ef[i + 1]);
      const aE = Math.abs(E);
      const EN = aE / (N * TD);
      this.ENcell[i] = EN;
      this.Ecell[i] = E;

      const j = this._rtIdx(EN), w = this._rtW, rt = this.rt;
      const kizN = rt.kiz[j] + w * (rt.kiz[j + 1] - rt.kiz[j]);   // м^3/с
      const ka2N = rt.ka2[j] + w * (rt.ka2[j + 1] - rt.ka2[j]);   // 1/с
      const ka3N = rt.ka3[j] + w * (rt.ka3[j + 1] - rt.ka3[j]);   // 1/с
      const kdsN = rt.kds[j] + w * (rt.kds[j + 1] - rt.kds[j]);   // 1/с
      const kexN = rt.kex[j] + w * (rt.kex[j + 1] - rt.kex[j]);   // 1/с
      const kr2 = rt.kr2[j] + w * (rt.kr2[j + 1] - rt.kr2[j]);
      const kr4 = rt.kr4[j] + w * (rt.kr4[j + 1] - rt.kr4[j]);
      const kd6N = rt.kd6[j] + w * (rt.kd6[j + 1] - rt.kd6[j]);   // 1/с
      const kizNN = kizN * N;                                     // 1/с

      const e0 = ne[i], p2 = nO2p[i], p4 = nO4p[i];
      const m1 = nOm[i], m2 = nO2m[i], m3 = nO3m[i];
      const a0 = nO[i], o3 = nO3[i], sa = nO2a[i];

      // ---- псевдо-первопорядковые коэффициенты (нейтральные партнёры заморожены) --
      const cD1 = K_D1 * sa, cD2 = K_D2 * a0, cD5 = this.kD5N;      // O-  -> e
      const cD3 = K_D3 * sa, cD4 = K_D4 * a0, cD6 = kd6N;           // O2- -> e
      const cD7 = K_D7 * a0;                                        // O3- -> e
      const cC3 = this.kC3N2, cC4 = K_C4 * o3;                      // O-  -> O3-
      const cC5 = K_C5 * o3;                                        // O2- -> O3-
      const cC6 = K_C6 * a0;                                        // O3- -> O2-
      const cE8 = K_E8 * o3;                                        // e + O3 -> O-
      const mTot = m1 + m2 + m3;

      // ---- стадия 1: средние по шагу плотности ------------------------------
      const Pe = p.sBg + m1 * (cD1 + cD2 + cD5) + m2 * (cD3 + cD4 + cD6) + m3 * cD7;
      const Le = ka2N + ka3N + cE8 + kr2 * p2 + kr4 * p4 - kizNN;
      const nuGrow = -Le;
      if (nuGrow > maxNu) maxNu = nuGrow;
      const be = nbarOf(e0, Pe, Le);
      this.ionizRate[i] = kizNN * be;

      const bp2 = nbarOf(p2, kizNN * be + p.sBg + this.kC2N * p4,
        kr2 * be + this.kC1N2 + kII * mTot);
      const bp4 = nbarOf(p4, this.kC1N2 * bp2, kr4 * be + this.kC2N + kII * mTot);
      const bpT = bp2 + bp4;
      const bm1 = nbarOf(m1, (ka2N + cE8) * be, cD1 + cD2 + cD5 + cC3 + cC4 + kII * bpT);
      const bm2 = nbarOf(m2, ka3N * be + cC6 * m3, cD3 + cD4 + cD6 + cC5 + kII * bpT);
      const bm3 = nbarOf(m3, (cC3 + cC4) * bm1 + cC5 * bm2, cD7 + cC6 + kII * bpT);

      // ---- стадия 2: экстенты реакций (все >= 0) ----------------------------
      let Riz = kizNN * be * dt;            // e + O2 -> 2e + O2+
      let Rbg = p.sBg * dt;                 // фон: -> e + O2+ (ПАРА, а не один электрон)
      let Ra2 = ka2N * be * dt;             // e + O2 -> O- + O
      let Ra3 = ka3N * be * dt;             // e + 2O2 -> O2- + O2
      let Re8 = cE8 * be * dt;              // e + O3 -> O- + O2
      let Rr2 = kr2 * be * bp2 * dt;        // e + O2+ -> нейтралы
      let Rr4 = kr4 * be * bp4 * dt;        // e + O4+ -> нейтралы
      let RC1 = this.kC1N2 * bp2 * dt;      // O2+ -> O4+
      let RC2 = this.kC2N * bp4 * dt;       // O4+ -> O2+
      let RD1 = cD1 * bm1 * dt, RD2 = cD2 * bm1 * dt, RD5 = cD5 * bm1 * dt;
      let RD3 = cD3 * bm2 * dt, RD4 = cD4 * bm2 * dt, RD6 = cD6 * bm2 * dt;
      let RD7 = cD7 * bm3 * dt;
      let RC3 = cC3 * bm1 * dt, RC4 = cC4 * bm1 * dt;
      let RC5 = cC5 * bm2 * dt, RC6 = cC6 * bm3 * dt;
      const kIIdt = kII * dt;
      // ион-ионная рекомбинация: ОДИН экстент на пару (раньше положительные и
      // отрицательные считали её независимо и теряли заряд)
      let Ri21 = kIIdt * bp2 * bm1, Ri22 = kIIdt * bp2 * bm2, Ri23 = kIIdt * bp2 * bm3;
      let Ri41 = kIIdt * bp4 * bm1, Ri42 = kIIdt * bp4 * bm2, Ri43 = kIIdt * bp4 * bm3;

      // ---- ограничение по запасу реагентов (применяется к ЭКСТЕНТАМ) --------
      const Ge = Riz + Rbg + RD1 + RD2 + RD5 + RD3 + RD4 + RD6 + RD7;
      const Ce = Ra2 + Ra3 + Re8 + Rr2 + Rr4;
      const Gp2 = Riz + Rbg + RC2, Cp2 = Rr2 + RC1 + Ri21 + Ri22 + Ri23;
      const Gp4 = RC1, Cp4 = Rr4 + RC2 + Ri41 + Ri42 + Ri43;
      const Gm1 = Ra2 + Re8, Cm1 = RD1 + RD2 + RD5 + RC3 + RC4 + Ri21 + Ri41;
      const Gm2 = Ra3 + RC6, Cm2 = RD3 + RD4 + RD6 + RC5 + Ri22 + Ri42;
      const Gm3 = RC3 + RC4 + RC5, Cm3 = RD7 + RC6 + Ri23 + Ri43;
      const fe = lim(e0 + Ge, Ce), fp2 = lim(p2 + Gp2, Cp2), fp4 = lim(p4 + Gp4, Cp4);
      const fm1 = lim(m1 + Gm1, Cm1), fm2 = lim(m2 + Gm2, Cm2), fm3 = lim(m3 + Gm3, Cm3);
      if (fe < 1) { Ra2 *= fe; Ra3 *= fe; Re8 *= fe; Rr2 *= fe; Rr4 *= fe; }
      if (fp2 < 1) { Rr2 *= fp2; RC1 *= fp2; Ri21 *= fp2; Ri22 *= fp2; Ri23 *= fp2; }
      if (fp4 < 1) { Rr4 *= fp4; RC2 *= fp4; Ri41 *= fp4; Ri42 *= fp4; Ri43 *= fp4; }
      if (fm1 < 1) { RD1 *= fm1; RD2 *= fm1; RD5 *= fm1; RC3 *= fm1; RC4 *= fm1; Ri21 *= fm1; Ri41 *= fm1; }
      if (fm2 < 1) { RD3 *= fm2; RD4 *= fm2; RD6 *= fm2; RC5 *= fm2; Ri22 *= fm2; Ri42 *= fm2; }
      if (fm3 < 1) { RD7 *= fm3; RC6 *= fm3; Ri23 *= fm3; Ri43 *= fm3; }

      // ---- применение: n = n0 + приход - расход (одни и те же ΔR) -----------
      ne[i] = e0 + (Riz + Rbg + RD1 + RD2 + RD5 + RD3 + RD4 + RD6 + RD7)
        - (Ra2 + Ra3 + Re8 + Rr2 + Rr4);
      nO2p[i] = p2 + (Riz + Rbg + RC2) - (Rr2 + RC1 + Ri21 + Ri22 + Ri23);
      nO4p[i] = p4 + RC1 - (Rr4 + RC2 + Ri41 + Ri42 + Ri43);
      nOm[i] = m1 + (Ra2 + Re8) - (RD1 + RD2 + RD5 + RC3 + RC4 + Ri21 + Ri41);
      nO2m[i] = m2 + (Ra3 + RC6) - (RD3 + RD4 + RD6 + RC5 + Ri22 + Ri42);
      nO3m[i] = m3 + (RC3 + RC4 + RC5) - (RD7 + RC6 + Ri23 + Ri43);

      // --- накопление источников нейтралов (медленная химия, §4.4) ---
      const Rds = kdsN * be * dt;           // e + O2 -> 2O + e
      this.accO[i] += 2 * Rds + Ra2 + 2 * Rr2 + (Ri21 + Ri41) + RC4;
      this.accO3[i] += this.kN1N2 * a0 * dt + RD1 + RD5 + RD4 + (Ri23 + Ri43);
      this.accO2a[i] += kexN * be * dt;
    }

    // медленное расщепление для нейтралов
    this.tAcc += dt;
    if (this.tAcc >= 1e-7) this._slowChemistry(this.tAcc);
    return maxNu;
  }

  /** Медленная химия нейтралов O, O3, O2(a) с накопленным производством. */
  _slowChemistry(dts) {
    const i0 = this.i0, i1 = this.i1;
    const nO = this.nO, nO3 = this.nO3, nO2a = this.nO2a;
    const nOm = this.n[S_OM], nO2m = this.n[S_O2M], nO3m = this.n[S_O3M], ne = this.n[S_E];
    const inv = 1 / dts;
    // пристеночная гибель — как объёмная в пристеночных ячейках
    const wO = (0.25 * SPECIES_BY_ID.O.gammaWall * this.vthO);
    const wO3 = (0.25 * SPECIES_BY_ID.O3.gammaWall * this.vthO3);
    const wO2a = (0.25 * SPECIES_BY_ID.O2a.gammaWall * this.vthO2a);

    for (let i = i0; i <= i1; i++) {
      const a0 = nO[i], o3 = nO3[i], sa = nO2a[i];
      const wall = (i === i0 || i === i1) ? 1 / this.dx[i] : 0;

      // O
      {
        const P = this.accO[i] * inv + this.kN4v * o3 * sa;
        const L = this.kN1N2 + 2 * this.kN2N * a0 + this.kN3v * o3
          + K_C6 * nO3m[i] + K_D2 * nOm[i] + K_D4 * nO2m[i] + K_D7 * nO3m[i]
          + K_N6 * sa + wall * wO;
        const z = Math.min(30, L * dts);
        nO[i] = a0 * Math.exp(-z) + P * phi1(z) * dts;
      }
      // O3
      {
        const P = this.accO3[i] * inv;
        const L = K_E8 * ne[i] + K_C4 * nOm[i] + K_C5 * nO2m[i]
          + this.kN3v * a0 + this.kN4v * sa + wall * wO3;
        const z = Math.min(30, L * dts);
        nO3[i] = o3 * Math.exp(-z) + P * phi1(z) * dts;
      }
      // O2(a)
      {
        const P = this.accO2a[i] * inv;
        const L = this.kN5N + K_N6 * a0 + this.kN4v * o3 + wall * wO2a;
        const z = Math.min(30, L * dts);
        nO2a[i] = sa * Math.exp(-z) + P * phi1(z) * dts;
      }
      this.accO[i] = 0; this.accO3[i] = 0; this.accO2a[i] = 0;
    }
    this.tAcc = 0;
  }

  // -------------------------------------------------------------------------
  // Контроллер шага
  // -------------------------------------------------------------------------

  _limitDt(dt, maxIonV, maxWe) {
    // жёсткие ограничения устойчивости явной части (ионы) — применяются ДО транспорта
    let lim = 2e-8;
    if (isFinite(this._ionCfl)) lim = Math.min(lim, this.cfl * this._ionCfl);
    // электроны неявные, но точность фронта лавины требует умеренного CFL там,
    // где их плотность значима
    if (this._eSignificant && isFinite(this._eCflT)) lim = Math.min(lim, this.eCFL * this._eCflT);
    if (dt > lim) dt = lim;
    if (dt < 1e-14) dt = 1e-14;
    return dt;
  }

  _nextDt(dt, maxIonV, maxWe, maxNu) {
    const p = this.params;
    let next = dt * 1.25;
    // ионный CFL
    if (isFinite(this._ionCfl)) next = Math.min(next, this.cfl * this._ionCfl);
    // значимость электронов
    let neMax = 0;
    const ne = this.n[S_E];
    for (let i = this.i0; i <= this.i1; i++) if (ne[i] > neMax) neMax = ne[i];
    this._eSignificant = neMax > 1e14;
    if (this._eSignificant && isFinite(this._eCflT)) next = Math.min(next, this.eCFL * this._eCflT);
    // Точность лавины. Ограничитель применяется по ЧИСТОЙ СКОРОСТИ РОСТА (nu_iz - nu_att)
    // и БЕЗ гейта по плотности: раньше он был закрыт условием n_e > 1e14 и потому не
    // работал именно на фронте лавины, где z = nu*dt доходило до 1.4. Линейные каналы
    // потерь интегрируются точно, поэтому ограничивать нужно только рост.
    if (maxNu > 0) next = Math.min(next, 0.5 / maxNu);
    // изменение поля
    const maxDErel = this._maxDErel || 0;
    if (maxDErel > 0) next = Math.min(next, (dt * this.dETol) / maxDErel);
    // напряжение источника: не более 1/400 периода
    next = Math.min(next, this.period / 400);
    if (next > 2e-8) next = 2e-8;
    if (next < 1e-14) next = 1e-14;
    return next;
  }

  // -------------------------------------------------------------------------
  // Диагностика
  // -------------------------------------------------------------------------

  _updateDiagnostics() {
    const i0 = this.i0, i1 = this.i1, N = this.N;
    // поле и E/N по ячейкам (включая диэлектрики — для отрисовки)
    for (let i = 0; i < this.Nt; i++) {
      const E = 0.5 * (this.Ef[i] + this.Ef[i + 1]);
      this.Ecell[i] = E;
      this.ENcell[i] = this.gasMask[i] ? Math.abs(E) / (N * TD) : 0;
    }
    let mx = 0, o3sum = 0;
    for (let i = i0; i <= i1; i++) {
      if (this.ENcell[i] > mx) mx = this.ENcell[i];
      o3sum += this.nO3[i] * this.dx[i];
    }
    this.maxEN = mx;
    this.o3ppm = (o3sum / this.Lg / N) * 1e6;
    // напряжение на зазоре: потенциалы поверхностей
    this.Ugap = this._surfacePotential(this.fL) - this._surfacePotential(this.fR);
    const ai = Math.abs(this.condCurrent);
    if (ai > this.peakCurrent) this.peakCurrent = ai;
    this.power = this.Uel * this.current;
  }

  _surfacePotential(f) {
    if (f === 0) return this.Uel;
    if (f === this.Nt) return 0;
    const iL = f - 1, iR = f;
    const gL = (2 * this.epsCell[iL]) / this.dx[iL];
    const gR = (2 * this.epsCell[iR]) / this.dx[iR];
    return (this.sigFace[f] + gL * this.phi[iL] + gR * this.phi[iR]) / (gL + gR);
  }

  _record() {
    const H = this.history;
    const dtRec = this.period / 2000;
    if (this.t - this.lastHistT < dtRec) return;
    this.lastHistT = this.t;
    const k = H.head;
    H.t[k] = this.t; H.Uapp[k] = this.Uel; H.Ugap[k] = this.Ugap;
    H.current[k] = this.current; H.condCurrent[k] = this.condCurrent;
    H.charge[k] = this.charge;
    H.head = (k + 1) % H.capacity;
    if (H.len < H.capacity) H.len++;
  }

  _periodBookkeeping() {
    // Выборка Q–V. ВАЖНО: заполняем ВЕСЬ диапазон [qvLen, slot], а не только slot.
    // Шаг по времени спокойно перепрыгивает несколько слотов (при f=100 кГц — до 10),
    // и раньше пропущенные ячейки оставались с данными ПРОШЛОГО периода, но входили
    // и в qvLen, и в замкнутый интеграл ∮U dQ (ошибка мощности до 25%), и в маску
    // горения (C_cell ошибался в 5.6 раза).
    const frac = (this.t - this.periodT0) / this.period;
    const slot = Math.min(QV_CAP - 1, Math.floor(frac * QV_CAP));
    if (slot >= this.qvLen) {
      for (let k = this.qvLen; k <= slot; k++) {
        this.qvU[k] = this.Uel; this.qvQ[k] = this.charge;
        this.qvI[k] = this.condCurrent; this.qvUg[k] = this.Ugap;
      }
      this.qvLen = slot + 1;
    }
    // счётчик импульсов
    const thr = Math.max(1e-6, 0.05 * this.peakCurrent);
    const on = Math.abs(this.condCurrent) > thr;
    if (on && !this.pulseOn) this.pulseCount++;
    this.pulseOn = on;

    if (this.t - this.periodT0 >= this.period) {
      this._finishPeriod();
      this.periodT0 = this.t;
      this.periodIndex++;
      this.qvLen = 0;
      this.breakdownsPerPeriod = this.pulseCount;
      this.pulseCount = 0;
      this.peakCurrent *= 0.5;  // мягкий сброс пиковой метрики
    }
  }

  /** Метод Мэнли: энергия за период и наклоны параллелограмма Лиссажу. */
  _finishPeriod() {
    const n = this.qvLen;
    if (n < 64) return;
    // энергия = площадь фигуры Q–V = ∮ U dQ
    let E = 0, comp = 0;
    for (let k = 0; k < n; k++) {
      const k2 = (k + 1) % n;
      const term = 0.5 * (this.qvU[k] + this.qvU[k2]) * (this.qvQ[k2] - this.qvQ[k]);
      const y = term - comp; const t2 = E + y; comp = (t2 - E) - y; E = t2;
    }
    E = Math.abs(E);
    const sl = this._lissajousSlopes(n);
    let ugmax = 0;
    for (let k = 0; k < n; k++) ugmax = Math.max(ugmax, Math.abs(this.qvUg[k]));
    this.periodStats = {
      energyPerPeriodJ: E,
      powerW: E * this.freq,
      // Cdiel — ИЗМЕРЕННЫЙ наклон горящей ветви dQ/dU_app (то, что реально снимают
      // с фигуры Лиссажу). CdielIdentity — величина dQ/d(U_app-U_gap): это ТОЖДЕСТВО
      // последовательной цепи Q = C_d*(U_app - U_gap), оно равно C_d на ЛЮБОЙ дуге,
      // включая полностью тёмную, и потому годится только как самопроверка схемы.
      Cdiel: sl.CdielRaw, CdielRaw: sl.CdielRaw, CdielIdentity: sl.CdielIdentity,
      CdielGeom: this.Cd, CcellGeom: this.Ccell, Ccell: sl.Ccell,
      qvArcsOn: sl.arcsOn, qvArcsOff: sl.arcsOff, qvFallback: sl.fallback,
      qvR2on: sl.r2on, qvR2off: sl.r2off, qvOk: sl.ok,
      // напряжение горения = потолок, на котором зажат зазор (V2/V8)
      UburnkV: ugmax / 1e3,
    };
    this._storeQvGhost(n);
  }

  /** Сохранить огибающую Q–V завершённого периода (для «призраков» в UI). */
  _storeQvGhost(n) {
    const H = this.qvHistory, k = H.head;
    const u = H.u[k], q = H.q[k];
    for (let j = 0; j < QV_GHOST_N; j++) {
      const src = Math.min(n - 1, Math.floor((j * n) / QV_GHOST_N));
      u[j] = this.qvU[src]; q[j] = this.qvQ[src];
    }
    H.t0[k] = this.periodT0;
    H.head = (k + 1) % QV_GHOSTS;
    if (H.n < QV_GHOSTS) H.n++;
  }

  /**
   * Наклоны сторон параллелограмма Лиссажу (метод Мэнли).
   *
   * Две ловушки, обе проверены на реальных данных:
   *  1) Общая регрессия Q(U) по всем «тёмным» точкам НЕ работает: тёмных дуг
   *     за период две (плюс разрыв на пике U), они лежат на параллельных прямых
   *     с разными смещениями, и общая регрессия даёт наклон где-то посередине.
   *  2) Локальный (скользящий) наклон тоже НЕ работает на горящей ветви: разряд
   *     идёт пачкой коротких импульсов, между которыми наклон ~C_cell, а внутри
   *     импульса ~40 пФ. C_diel — это наклон ОГИБАЮЩЕЙ всей горящей полуволны.
   *
   * Поэтому: помечаем «горит/не горит», морфологически замыкаем паузы между
   * микроразрядами, режем период на СВЯЗНЫЕ дуги и делаем регрессию внутри
   * каждой дуги отдельно (внутри дуги смещение общее — регрессия корректна).
   */
  _lissajousSlopes(n) {
    const fail = {
      CdielRaw: this.Cd, CdielIdentity: this.Cd, Ccell: this.Ccell,
      arcsOn: 0, arcsOff: 0, r2on: 0, r2off: 0, ok: false,
    };
    let imax = 0, umax = 0;
    for (let k = 0; k < n; k++) {
      imax = Math.max(imax, Math.abs(this.qvI[k]));
      umax = Math.max(umax, Math.abs(this.qvU[k]));
    }
    if (imax <= 0 || umax <= 0) return fail;
    // Порог «горит» — относительный И абсолютный. Без абсолютного даже полностью
    // тёмный прогон (U0 ниже пробоя) размечался на «горящие» дуги по шуму
    // кондукционного тока, и с них снимался бессмысленный C_diel.
    // «Тёмная» дуга по определению та, где кондукционный ток пренебрежим по
    // сравнению с током смещения холодной ячейки C_cell*dU/dt — только тогда
    // dQ/dU_app = C_cell. Чисто относительный порог 0.02*max|I| на высокой частоте
    // (где пики в десятки ампер) пропускал в «тёмную» ветвь заметное остаточное
    // горение и завышал C_cell на 38%.
    const iScale = this.Ccell * this.U0 * this.omega;   // масштаб тока смещения
    const thr = 0.05 * iScale;
    const on = this._onMask && this._onMask.length >= n
      ? this._onMask : (this._onMask = new Uint8Array(QV_CAP));
    for (let k = 0; k < n; k++) on[k] = Math.abs(this.qvI[k]) > thr ? 1 : 0;
    // замыкание пауз между микроразрядами (радиус ~1/24 периода = 15 град.)
    const rad = Math.max(2, Math.round(n / 24));
    const dil = this._onDil && this._onDil.length >= n
      ? this._onDil : (this._onDil = new Uint8Array(QV_CAP));
    for (let k = 0; k < n; k++) {
      let v = 0;
      for (let j = -rad; j <= rad && !v; j++) v = on[((k + j) % n + n) % n];
      dil[k] = v;
    }
    // связные дуги (по кругу): стартуем с точки смены состояния
    let start = 0;
    for (let k = 0; k < n; k++) if (dil[k] !== dil[(k + n - 1) % n]) { start = k; break; }
    let sumOn = 0, wOn = 0, sumOff = 0, wOff = 0, sumOnC = 0, wOnC = 0;
    let arcsOn = 0, arcsOff = 0, r2on = 0, r2off = 0;
    // Тёмная ветвь принимается в фит только если это ДЛИННАЯ дуга с заметным размахом
    // напряжения: короткий «огрызок» у пика |U| (dQ/dU там не равен C_cell) давал
    // ошибку C_cell в 5.6 раза на 50 кГц.
    const minOff = Math.max(8, Math.round(n / 10));
    const minSpanOff = 0.25 * umax;
    let k0 = 0;
    while (k0 < n) {
      const st = dil[(start + k0) % n];
      let len = 1;
      while (k0 + len < n && dil[(start + k0 + len) % n] === st) len++;
      if (st && len >= Math.max(8, n >> 5)) {
        const f = this._fitArc(start + k0, len, n, false);
        const fc = this._fitArc(start + k0, len, n, true);
        // измеренный наклон принимаем КАК ЕСТЬ (это и есть результат измерения),
        // R² только протоколируется и определяет флаг качества qvOk
        if (f.slope > 0) { sumOn += f.slope * len; wOn += len; r2on += f.r2 * len; arcsOn++; }
        if (fc.slope > 0) { sumOnC += fc.slope * len; wOnC += len; }
      } else if (!st && len >= minOff) {
        const f = this._fitArc(start + k0, len, n, false);
        if (f.slope > 0 && f.span > minSpanOff && f.r2 > 0.98) {
          // Внутри тёмной дуги берём МИНИМАЛЬНЫЙ локальный наклон по скользящему
          // окну: любое остаточное горение может наклон только увеличить, поэтому
          // C_cell — это нижняя огибающая. Без этого на 50 кГц, где тёмные дуги
          // короткие и подпорчены хвостами микроразрядов, C_cell завышался на 38%.
          const wWin = Math.max(8, Math.round(n / 20));
          let best = f.slope, bestR2 = f.r2;
          for (let q = 0; q + wWin <= len; q += Math.max(1, wWin >> 1)) {
            const g = this._fitArc(start + k0 + q, wWin, n, false);
            if (g.slope > 0 && g.r2 > 0.995 && g.span > 0.05 * umax && g.slope < best) {
              best = g.slope; bestR2 = g.r2;
            }
          }
          sumOff += best * len; wOff += len; r2off += bestR2 * len; arcsOff++;
        }
      }
      k0 += len;
    }
    return {
      // измеренный наклон горящей ветви — это и есть «C_diel по Мэнли»
      CdielRaw: wOn > 0 ? sumOn / wOn : this.Cd,
      // тождество последовательной цепи Q = C_d*(U_app - U_gap): диагностика схемы
      CdielIdentity: wOnC > 0 ? sumOnC / wOnC : this.Cd,
      Ccell: wOff > 0 ? sumOff / wOff : this.Ccell,
      arcsOn, arcsOff,
      // true, если наклон(и) не удалось снять и подставлена геометрическая константа
      fallback: wOn === 0 || wOff === 0,
      r2on: wOn > 0 ? r2on / wOn : 0,
      r2off: wOff > 0 ? r2off / wOff : 0,
      // фигура считается пригодной для метода Мэнли, только если найдены обе горящие
      // и обе тёмные стороны и обе легли на прямую
      ok: arcsOn >= 2 && arcsOff >= 2 && wOn > 0 && wOff > 0
        && r2on / Math.max(1, wOn) > 0.98 && r2off / Math.max(1, wOff) > 0.98,
    };
  }

  /**
   * МНК-наклон на связной дуге [k0, k0+len) (индексы по модулю n).
   * gapCorrected=false -> dQ/dU_app; true -> dQ/d(U_app - U_gap).
   * Возвращает {slope, r2, span} — R² и размах по x нужны для оценки качества фигуры.
   */
  _fitArc(k0, len, n, gapCorrected) {
    let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0, xmin = Infinity, xmax = -Infinity;
    for (let j = 0; j < len; j++) {
      const k = (k0 + j) % n;
      const u = gapCorrected ? this.qvU[k] - this.qvUg[k] : this.qvU[k];
      const q = this.qvQ[k];
      sx += u; sy += q; sxx += u * u; sxy += u * q; syy += q * q;
      if (u < xmin) xmin = u; if (u > xmax) xmax = u;
    }
    const den = len * sxx - sx * sx;
    if (Math.abs(den) < 1e-30) return { slope: 0, r2: 0, span: 0 };
    const slope = (len * sxy - sx * sy) / den;
    const dy = len * syy - sy * sy;
    const cov = len * sxy - sx * sy;
    const r2 = dy > 1e-300 ? (cov * cov) / (den * dy) : 0;
    return { slope, r2, span: xmax - xmin };
  }

  // -------------------------------------------------------------------------
  // Публичный API
  // -------------------------------------------------------------------------

  advance(wallClockBudgetMs = 8) {
    const t0 = now();
    let steps = 0;
    const tStart = this.t;
    do {
      this.step();
      steps++;
      if (this.t >= this.targetSimTime) break;
    } while (now() - t0 < wallClockBudgetMs);
    return { steps, simTime: this.t - tStart };
  }

  get state() {
    return {
      t: this.t, dt: this.dt,
      x: this.x, xFaces: this.xf, gasMask: this.gasMask,
      n: {
        e: this.n[S_E], O2p: this.n[S_O2P], O4p: this.n[S_O4P],
        Om: this.n[S_OM], O2m: this.n[S_O2M], O3m: this.n[S_O3M],
        O: this.nO, O3: this.nO3, O2a: this.nO2a,
      },
      E: this.Ecell, phi: this.phi, rho: this.rho,
      ionizRate: this.ionizRate, EN: this.ENcell,
      sigmaL: this.sigL, sigmaR: this.sigR,
      Uapp: this.Uel, Ugap: this.Ugap,
      current: this.current, condCurrent: this.condCurrent,
      dispCurrent: this.dispCurrent, charge: this.charge, power: this.power,
      o3ppm: this.o3ppm, maxEN: this.maxEN,
      peakCurrent: this.peakCurrent, breakdownsPerPeriod: this.breakdownsPerPeriod,
      // численное здоровье (расширение контракта): доля отбракованных шагов,
      // компенсирующий заряд пола и невязка Пуассона последнего шага
      steps: this.steps, rejects: this.rejects,
      clipCharge: this.clipCharge * this.A,
      residual: this._resid || 0, wallCfl: this.wallCfl || 0,
    };
  }

  /** Полный заряд системы (объём + обе поверхности), Кл. Для теста сохранения. */
  totalCharge() {
    let q = 0;
    for (let s = 0; s < NCH; s++) {
      const ns = this.n[s], z = QE * CH_Z[s];
      for (let i = this.i0; i <= this.i1; i++) q += z * ns[i] * this.dx[i];
    }
    return (q + this.sigL + this.sigR) * this.A;
  }

  /** Невязка Пуассона max|div(eps grad phi) + rho| (диагностика C1). */
  poissonResidual() {
    this._computeRho();
    this.sigCond.fill(0);
    this._assembleMatrix(0);
    this._assembleRhs(this.Uel, true, 0);
    let mx = 0, scale = 0;
    for (let i = 0; i < this.Nt; i++) {
      const l = this.Mb[i] * this.phi[i]
        + (i > 0 ? this.Ma[i] * this.phi[i - 1] : 0)
        + (i < this.Nt - 1 ? this.Mc[i] * this.phi[i + 1] : 0);
      mx = Math.max(mx, Math.abs(l - this.Mr[i]));
      scale = Math.max(scale, Math.abs(this.Mr[i]));
    }
    return mx / (scale + 1e-30);
  }
}

export { sgFlux, bernoulli, thomas, phi1 };
export default DBDSolver;
