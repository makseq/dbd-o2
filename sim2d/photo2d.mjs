// photo2d.mjs — фотопроцессы 2D осесимметричной модели ДБР в чистом O2.
//
// Главный источник: docs/PHOTO_PROCESSES.md (§2, §3, §5, §7), ERRATA.md раздел F.
// Иерархия документов: ERRATA.md > PHOTO_PROCESSES.md > PHYSICS.md/NUMERICS_2D.md.
//
// Три независимых канала, каждый под своим флагом:
//   1) photoIonization — окно A (< 102.8 нм), три группы Гельмгольца (§2.6, §7.2);
//   2) photoEmission   — окна A и B, прямой view factor на ОБЕ поверхности (§3, §7.4);
//   3) photoDetachment — окно C (> 200 нм, газ прозрачен), прямое ядро 1/(4piR^2) (§5, §7.5).
//
// Всё в СИ. Только typed arrays. Никаких аллокаций в горячем цикле (update()):
// вся память (матрицы view factor, ядро прозрачного переноса, рабочие буферы)
// выделяется один раз в конструкторе / при первом build().
//
// ВНИМАНИЕ (ERRATA F3 / PHOTO_PROCESSES §7.1): уравнение Гельмгольца собирается
// ОТДЕЛЬНО, с eps = 1, только на газовой части сетки — НЕ сдвигом диагонали матрицы
// Пуассона (та собрана с eps(z), eps_r = 9 в барьерах). Здесь это выражено тем, что
// модуль НИКОГДА не трогает solvePoisson(), а вызывает только solveHelmholtz().
// Контракт solveHelmholtz описан в sim2d/PHOTO_API.md.

const EV = 1.602176634e-19;          // Дж/эВ
const QE = 1.602176634e-19;          // Кл

// ---------------------------------------------------------------------------
// Наборы весов трёхгруппового разложения окна A.
//
// lambda_j переносятся ЗАКОННО [B]: они описывают поглощение O2 и масштабируются
// с парциальным давлением O2 (воздух -> чистый O2 = сокращение всех длин в 5 раз).
// Веса A_j/lambda_j^2 переносятся НЕЗАКОННО [C]: они кодируют форму АЗОТНОГО
// эмиссионного спектра в окне 98–102.5 нм. Поэтому веса — параметр модели,
// а не константа (PHOTO_PROCESSES.md §2.6).
// ---------------------------------------------------------------------------
export const PHOTO_ION_WEIGHT_SETS = {
  air:   [0.07, 0.26, 0.67],   // Бурдон при p_O2 = 760 Торр — дефолт
  short: [0.02, 0.13, 0.85],   // гипотеза «спектр смещён к 98 нм»
  long:  [0.20, 0.40, 0.40],   // гипотеза «спектр у порога 102.7 нм»
};

// lambda_j [м^-1] — обратные длины поглощения: 238 / 90 / 14.8 мкм.
export const PHOTO_ION_LAMBDA = [4.202e3, 1.111e4, 6.757e4];

export const PHOTO_DEFAULTS = {
  // --- флаги ---
  photoIonization: true,
  photoEmission:   true,
  photoDetachment: true,

  // --- окно A: фотоионизация + жёсткая фотоэмиссия ---
  etaGamma:        1e-3,        // фотонов окна A на одну ударную ионизацию [1e-5..1e-1] [C]
  PhiIon:          0.30,        // квантовый выход ионизации на поглощённый фотон [0.05..0.9] [B]
  photoIonWeights: 'air',       // 'air' | 'short' | 'long' | [w1,w2,w3]
  photoIonLambda:  PHOTO_ION_LAMBDA,
  hvA:             12.3,        // эВ, средняя энергия фотона окна A

  // --- фотоэмиссия ---
  Eth:             9.0,         // эВ, порог выхода Al2O3 [5.0..10.5] [B/C]
  Yref:            5e-3,        // выход при Eref [1e-4..5e-2] [C]
  Eref:            12.5,        // эВ
  fVUV:            1e-4,        // доля мощности разряда выше Eth [1e-6..5e-3] [C]
  Tw:              0.1,         // доля испущенных фотонов, доходящих до стенки [0.05..0.3]
  gammaEff:        0.02,        // «эффективный» гамма из PHYSICS.md §4.1
  Ub:              3.78e3,      // В, рабочая точка филамента (для оценки F_ph/F_i)
  hvB:             9.51,        // эВ, O I 130.4 нм
  ellB:            41e-6,       // м, длина поглощения на 130.4 нм
  ratioWindow:     'A',         // окно, по которому оценивается F_ph/F_i ('A' | 'B')

  // --- окно C: фотоотлипание ---
  fDet:            1e-4,        // доля мощности в прозрачное излучение [1e-6..1e-2] [C]
  hvC:             1.6,         // эВ, O I 777.4 нм
  sigmaPd:         [6.0e-22, 2.4e-22, 1.5e-22],  // O-, O2-, O3- [A/A/B]
  Cself:           0.72,        // самоячейка прозрачного ядра, PHOTO_PROCESSES §7.5

  // --- численные ---
  photoIonMode:     'helmholtz', // 'helmholtz' | 'local' | 'off'
  photoDetachMode:  'kernel',    // 'kernel' | 'uniform' | 'off'
  photoDetachNear:  true,        // ближняя поправка к огрублённому ядру
  photoNearFactor:  1.0,         // радиус ближней зоны в диагоналях огрублённой ячейки
  photoCoarse:      [32, 48],    // огрублённая сетка [ncr, ncz] для окна C
  nPhotoSubstep:    3,           // обновлять фотополя раз в N шагов
  vfAzimuthNodes:   64,          // базовое число узлов азимутальной квадратуры view factor
  vfDropRel:        1e-7,        // порог отбрасывания элемента матрицы view factor
  GslabRel:         1.21e-3,     // м, G_slab для режима 'uniform'
  photoCheck:       true,        // проверка finite/неотрицательности каждый апдейт (PH8)
  lazyBuild:        true,        // строить тяжёлые матрицы при первом update(), а не в конструкторе
};

// ---------------------------------------------------------------------------
// Фаулеровский выход фотоэмиссии (PHOTO_PROCESSES.md §3.2).
//   Y_ph(hv) = Y_ref * max(0, (hv - E_th)/(E_ref - E_th))^2
// max(0, ...) обязателен: именно он выключает группу PE-B при E_th > 9.51 эВ
// (критерий VF7). Без отсечки порог перестаёт быть качественным переключателем.
// ---------------------------------------------------------------------------
export function fowlerYield(hv, Eth, Yref, Eref) {
  const den = Eref - Eth;
  if (!(den > 0)) return 0;
  const x = (hv - Eth) / den;
  return x > 0 ? Yref * x * x : 0;
}

// ---------------------------------------------------------------------------
// Ядро прозрачного переноса (окно C), азимут взят АНАЛИТИЧЕСКИ.
//
//   Psi(r,z) = INT q(r',z') r' dr' dz'
//              / ( 2 * sqrt( [(r+r')^2 + dz^2] * [(r-r')^2 + dz^2] ) )
//
// Вывод: INT_0^{2pi} dphi/(a - b cos phi) = 2pi/sqrt(a^2 - b^2),
//   a = r^2 + r'^2 + dz^2, b = 2 r r',  a^2 - b^2 = [(r-r')^2+dz^2]*[(r+r')^2+dz^2].
// Проверяется тестом PH4 прямой квадратурой по углу (сходимость < 1e-6).
//
// Возвращает множитель при (r' dr' dz'), то есть 1/(2*sqrt(P1*P2)).
// ---------------------------------------------------------------------------
export function transparentKernel(r, z, rp, zp) {
  const dz2 = (z - zp) * (z - zp);
  const dm = r - rp, dp = r + rp;
  const P1 = dm * dm + dz2;
  const P2 = dp * dp + dz2;
  return 1 / (2 * Math.sqrt(P1 * P2));
}

// ---------------------------------------------------------------------------
// Нормализация описания сетки. Принимаем минимум {nr, nz, rf, zf, JG0, JG1},
// остальное выводим сами — чтобы модуль не зависел от внутренних имён solver2d.
// ---------------------------------------------------------------------------
function normalizeGrid(g) {
  const nr = g.nr, nz = g.nz;
  const rf = g.rf, zf = g.zf;
  if (!rf || rf.length !== nr + 1) throw new Error('photo2d: grid.rf must have nr+1 entries');
  if (!zf || zf.length !== nz + 1) throw new Error('photo2d: grid.zf must have nz+1 entries');
  const rc = g.rc ?? new Float64Array(nr);
  const dr = g.dr ?? new Float64Array(nr);
  const Acell = g.Acell ?? new Float64Array(nr);
  if (!g.rc) for (let i = 0; i < nr; i++) rc[i] = 0.5 * (rf[i] + rf[i + 1]);
  if (!g.dr) for (let i = 0; i < nr; i++) dr[i] = rf[i + 1] - rf[i];
  if (!g.Acell) for (let i = 0; i < nr; i++) Acell[i] = Math.PI * (rf[i + 1] * rf[i + 1] - rf[i] * rf[i]);
  const zc = g.zc ?? new Float64Array(nz);
  const dz = g.dz ?? new Float64Array(nz);
  if (!g.zc) for (let j = 0; j < nz; j++) zc[j] = 0.5 * (zf[j] + zf[j + 1]);
  if (!g.dz) for (let j = 0; j < nz; j++) dz[j] = zf[j + 1] - zf[j];
  const JG0 = g.JG0 ?? 0, JG1 = g.JG1 ?? nz - 1;
  return { nr, nz, rf, rc, dr, Acell, zf, zc, dz, JG0, JG1 };
}

export class PhotoModule {
  /**
   * @param {object} grid   {nr, nz, rf, zf, JG0, JG1, [rc, dr, Acell, zc, dz]}
   * @param {object} params см. PHOTO_DEFAULTS
   * @param {object} solver объект с solveHelmholtz(q, kappa, out) (см. PHOTO_API.md).
   *                        Обязателен только при photoIonization = true.
   */
  constructor(grid, params = {}, solver = null) {
    const G = normalizeGrid(grid);
    this.g = G;
    this.solver = solver;
    const p = { ...PHOTO_DEFAULTS, ...params };
    this.p = p;

    const nr = G.nr, nz = G.nz, N = nr * nz;
    this.N = N;
    this.jg0 = G.JG0;
    this.jg1 = G.JG1;
    this.ngz = G.JG1 - G.JG0 + 1;

    // --- веса окна A, нормировка sum w_j = 1 (VF4: sum A_j/lambda_j^2 = 1) ---
    let w = p.photoIonWeights;
    if (typeof w === 'string') {
      if (!PHOTO_ION_WEIGHT_SETS[w]) throw new Error(`photo2d: unknown photoIonWeights preset '${w}'`);
      w = PHOTO_ION_WEIGHT_SETS[w];
    }
    if (!w || w.length !== 3) throw new Error('photo2d: photoIonWeights must be 3 numbers or a preset name');
    let ws = 0;
    for (let j = 0; j < 3; j++) {
      if (!(w[j] >= 0)) throw new Error('photo2d: photoIonWeights must be non-negative');
      ws += w[j];
    }
    if (!(ws > 0)) throw new Error('photo2d: photoIonWeights sum must be > 0');
    this.wA = Float64Array.from(w, (x) => x / ws);
    this.lamA = Float64Array.from(p.photoIonLambda);
    // A_j = w_j * lambda_j^2 — так функция Грина exp(-lam R)/(4 pi R) даёт
    // INT S_j dV = w_j * INT q dV, а сумма по j сохраняет источник тождественно.
    this.Aj = new Float64Array(3);
    for (let j = 0; j < 3; j++) this.Aj[j] = this.wA[j] * this.lamA[j] * this.lamA[j];
    // Алгебраическая проверка VF4.
    let chk = 0;
    for (let j = 0; j < 3; j++) chk += this.Aj[j] / (this.lamA[j] * this.lamA[j]);
    if (Math.abs(chk - 1) > 1e-3) throw new Error(`photo2d: sum A_j/lambda_j^2 = ${chk}, must be 1`);

    // --- выходы фотоэмиссии по группам ---
    // Группы: PE-A1/A2/A3 (окно A, hv = 12.3 эВ, ell = 238/90/14.8 мкм) и
    //         PE-B (окно B, hv = 9.51 эВ, ell = 41 мкм).
    // §7.4: усреднённая ell_eff = 50 мкм ТЕРЯЕТ длинный хвост, который и даёт ореол
    // rho_50 = sqrt(3)*h, определяющий место следующего зажигания. Поэтому для окна A
    // сохраняем полное трёхэкспоненциальное разложение (три подматрицы), а не одну ell_eff.
    const Ya = fowlerYield(p.hvA, p.Eth, p.Yref, p.Eref);
    const Yb = fowlerYield(p.hvB, p.Eth, p.Yref, p.Eref);
    this.YA = Ya;
    this.YB = Yb;
    this.emitGroups = [];
    for (let j = 0; j < 3; j++) {
      this.emitGroups.push({ name: `PE-A${j + 1}`, window: 'A', mu: this.lamA[j], Y: Ya, wsrc: this.wA[j] });
    }
    this.emitGroups.push({ name: 'PE-B', window: 'B', mu: 1 / p.ellB, Y: Yb, wsrc: 1 });

    // ---------------------------------------------------------------------
    // ЗАПРЕТ ДВОЙНОГО СЧЁТА (ERRATA F2, PHOTO_PROCESSES.md §3.7).
    //
    // gamma_eff = 0.02 в PHYSICS.md объявлен ЭФФЕКТИВНЫМ: в него УЖЕ свёрнута
    // фотоэмиссия. После введения явного члена Gamma_pe оставлять gamma_i = 0.02
    // НЕЛЬЗЯ — вторичные электроны будут посчитаны дважды. Разделяем явно:
    //
    //     gamma_i = gamma_eff / (1 + F_ph/F_i)|рабочая точка
    //
    // Тогда полный вторичный поток
    //     gamma_i*Gamma_i + F_ph = gamma_i*Gamma_i*(1 + F_ph/F_i) = gamma_eff*Gamma_i
    // то есть при номинале интегральная доза не меняется, а меняется её ВРЕМЕННАЯ
    // структура (фотонная часть приходит за ~3 пс, ионная за 0.5–1.4 мкс) — ровно
    // то, ради чего канал и вводится (§3.6). Это тождество проверяется тестом PH7.
    //
    // Структурная форма отношения (§3.5), самосогласованная с U_b, E_ph, gamma:
    //     F_ph/F_i = Y_ph * f_VUV * T_w * (e*U_b) / (gamma_eff * E_ph)
    // Контрольная константа ревью: 1.05e4 * Y_ph * f_VUV * T_w / beta_i.
    //
    // При photoEmission = false фотонного члена нет, и правильный коэффициент —
    // именно gamma_eff (в нём фотоэмиссия «зашита»); никакого деления не делаем.
    // ---------------------------------------------------------------------
    const hvRef = p.ratioWindow === 'B' ? p.hvB : p.hvA;
    const Yref_w = p.ratioWindow === 'B' ? Yb : Ya;
    // e*U_b/(gamma_eff*E_ph) с E_ph в эВ: e сокращается, остаётся U_b[В]/(gamma*hv[эВ]).
    this.structFactor = p.Ub / (p.gammaEff * hvRef);
    this.FphOverFi = p.photoEmission ? Yref_w * p.fVUV * p.Tw * this.structFactor : 0;
    this.FphOverFiCheck = p.photoEmission ? 1.05e4 * Yref_w * p.fVUV * p.Tw : 0; // контроль ревью
    this.gammaI = p.photoEmission ? p.gammaEff / (1 + this.FphOverFi) : p.gammaEff;

    // --- выходные массивы (всегда существуют, всегда нули при выключенных флагах) ---
    this.photoIonRate  = new Float64Array(N);   // S_pi, м^-3 с^-1
    this.psi           = new Float64Array(N);   // Psi, поток прозрачных фотонов, м^-2 с^-1
    this.nuPdOm        = new Float64Array(N);   // частота фотоотлипания O-,  с^-1
    this.nuPdO2m       = new Float64Array(N);   // O2-
    this.nuPdO3m       = new Float64Array(N);   // O3-
    this.photoDetachRate = new Float64Array(N); // заполняется солвером или detachRate()
    this.photoEmitFluxL = new Float64Array(nr); // Gamma_pe на нижней поверхности, м^-2 с^-1
    this.photoEmitFluxR = new Float64Array(nr); // на верхней

    // --- рабочие буферы ---
    this._qA  = new Float64Array(N);   // источник окна A (фот/м^3/с)
    this._qB  = new Float64Array(N);   // источник окна B
    this._qC  = new Float64Array(N);   // источник окна C
    this._Sj  = new Float64Array(N);   // решение одной группы Гельмгольца
    this._scr = new Float64Array(N);   // общий скретч

    this._built = false;
    this._stepCount = -1;
    this._lastUpdateStep = -1;

    // Объёмы ячеек (нужны и view factor'у, и ядру окна C).
    this.V = new Float64Array(N);
    for (let i = 0; i < nr; i++) {
      const a = G.Acell[i];
      for (let j = 0; j < nz; j++) this.V[i * nz + j] = a * G.dz[j];
    }

    if (!p.lazyBuild) this.build();
  }

  /** Полное число ячеек газа (для диагностики). */
  get gasCellCount() { return this.g.nr * this.ngz; }

  // -------------------------------------------------------------------------
  // Построение тяжёлых предвычисленных структур. Идемпотентно.
  // -------------------------------------------------------------------------
  build() {
    if (this._built) return;
    const p = this.p;
    if (p.photoEmission) this._buildViewFactor();
    if (p.photoDetachment && p.photoDetachMode === 'kernel') this._buildTransparentKernel();
    this._built = true;
  }

  // =========================================================================
  // ФОТОЭМИССИЯ: прямой view factor (PHOTO_PROCESSES.md §7.4)
  //
  //   Gamma_g(r_w) = INT_V q_g(r',z') * exp(-mu_g R) * cos(theta) / (4 pi R^2) dV'
  //   cos(theta) = |z' - z_wall| / R
  //
  // В осесимметрии азимутальный интеграл сводится к эллиптическому — замкнутой
  // элементарной формы НЕТ (в отличие от объёмного ядра §7.5), поэтому считаем
  // его один раз при инициализации квадратурой по phi. Подынтегральная функция
  // периодическая и гладкая при R > dr, поэтому средняя точка (= трапеция для
  // периодической функции) сходится спектрально; число узлов адаптируется по
  // ширине пика 1/R^3.
  //
  // Матрица хранится в CSR по строкам-точкам стенки: элементы, дающие менее
  // vfDropRel от максимума строки, отбрасываются (для короткой группы это
  // экономит почти всё, для длинной ell = 238 мкм обрезание невозможно, §7.4).
  // =========================================================================
  _buildViewFactor() {
    const G = this.g, nr = G.nr, nz = G.nz;
    const zL = G.zf[this.jg0];          // нижняя поверхность диэлектрика
    const zR = G.zf[this.jg1 + 1];      // верхняя

    // Зеркальная симметрия газовой z-сетки (профили NUMERICS_2D §1.3 симметричны:
    // газ разбивается симметричным tanh). Тогда матрица верхней стенки — это
    // матрица нижней с отражённым источником, и её можно не хранить: вдвое меньше
    // памяти и вдвое меньше времени сборки. Проверяем численно, не на веру.
    let sym = true;
    for (let jj = 0; jj < this.ngz; jj++) {
      const a = G.zc[this.jg0 + jj] - zL;
      const b = zR - G.zc[this.jg1 - jj];
      if (Math.abs(a - b) > 1e-12 * (zR - zL)) { sym = false; break; }
    }
    this.vfMirror = sym;

    this.vfL = [];
    this.vfR = [];
    for (const grp of this.emitGroups) {
      // Группы с нулевым выходом не строим вовсе — это и экономия, и гарантия
      // строгого нуля вклада (VF7: при E_th = 10.5 группа PE-B обязана дать 0).
      if (!(grp.Y > 0)) { this.vfL.push(null); this.vfR.push(null); continue; }
      this.vfL.push(this._buildVFWall(grp.mu, zL));
      this.vfR.push(sym ? null : this._buildVFWall(grp.mu, zR));
    }
    if (sym) { this._qAm = new Float64Array(this.N); this._qBm = new Float64Array(this.N); }
  }

  /** Отражение источника по z внутри газа: j -> jg0 + jg1 - j. */
  _mirrorSrc(src, dst) {
    const nz = this.g.nz;
    dst.fill(0);
    for (let i = 0; i < this.g.nr; i++) {
      const b = i * nz;
      for (let j = this.jg0; j <= this.jg1; j++) dst[b + (this.jg0 + this.jg1 - j)] = src[b + j];
    }
    return dst;
  }

  _buildVFWall(mu, zw) {
    const G = this.g, nr = G.nr, nz = G.nz;
    const p = this.p;
    const rowPtr = new Int32Array(nr + 1);
    const cols = [];
    const vals = [];
    const rowVal = new Float64Array(nr * this.ngz);
    const rowCol = new Int32Array(nr * this.ngz);
    const invFourPi = 1 / (4 * Math.PI);

    // Предфильтр: upper(ip,jp) = значение ядра при R = R_min — это ВЕРХНЯЯ оценка
    // истинного веса (exp(-mu R) и 1/R^3 монотонно убывают по R). Клетки, у которых
    // даже верхняя оценка ниже порога строки, азимутальной квадратурой не считаем.
    // Для короткой группы (ell = 14.8 мкм) это снимает почти всю работу; для длинной
    // (238 мкм) обрезание невозможно в принципе (PHOTO_PROCESSES §7.4).
    const upperBuf = new Float64Array(nr * this.ngz);

    for (let iw = 0; iw < nr; iw++) {
      const rw = G.rc[iw];
      let cnt = 0, mx = 0, maxUpper = 0;
      for (let ip = 0; ip < nr; ip++) {
        const rp = G.rc[ip];
        const radW = G.Acell[ip] / (2 * Math.PI);
        for (let jj = 0; jj < this.ngz; jj++) {
          const jp = this.jg0 + jj;
          const h = Math.abs(G.zc[jp] - zw);
          const dmin = Math.abs(rw - rp);
          const Rmin = Math.sqrt(dmin * dmin + h * h);
          let u = 0;
          if (Rmin > 0) {
            u = Math.exp(-mu * Rmin) * h * invFourPi / (Rmin * Rmin * Rmin)
                * (2 * Math.PI) * radW * G.dz[jp];
          }
          upperBuf[ip * this.ngz + jj] = u;
          if (u > maxUpper) maxUpper = u;
        }
      }
      const thrPre = 0.1 * p.vfDropRel * maxUpper;
      for (let ip = 0; ip < nr; ip++) {
        const rp = G.rc[ip];
        // r' dr' проинтегрировано точно по ячейке: (rf+^2 - rf-^2)/2 = Acell/(2pi)
        const radW = G.Acell[ip] / (2 * Math.PI);
        for (let jj = 0; jj < this.ngz; jj++) {
          const jp = this.jg0 + jj;
          const h = Math.abs(G.zc[jp] - zw);
          const dzc = G.dz[jp];
          if (upperBuf[ip * this.ngz + jj] <= thrPre) continue;
          const dmin = Math.abs(rw - rp);
          const Rmin = Math.sqrt(dmin * dmin + h * h);
          // Адаптивная азимутальная квадратура: ширина пика по phi ~ Rmin/sqrt(rw*rp).
          let nq = p.vfAzimuthNodes;
          const sc = Math.sqrt(rw * rp);
          if (sc > 0) {
            const width = Math.max(Rmin / sc, 1e-6);
            nq = Math.min(512, Math.max(nq, Math.ceil(8 * Math.PI / width)));
          }
          // Средняя точка на [0, pi], удвоение по симметрии cos(phi).
          let acc = 0;
          const dphi = Math.PI / nq;
          const rr = rw * rw + rp * rp + h * h;
          const bb = 2 * rw * rp;
          for (let k = 0; k < nq; k++) {
            const phi = (k + 0.5) * dphi;
            const R2 = rr - bb * Math.cos(phi);
            const R = Math.sqrt(R2);
            acc += Math.exp(-mu * R) / (R2 * R);
          }
          acc *= 2 * dphi;                       // интеграл по [0, 2pi]
          const wgt = acc * h * invFourPi * radW * dzc;
          if (wgt > 0) {
            rowVal[cnt] = wgt;
            rowCol[cnt] = ip * nz + jp;
            cnt++;
            if (wgt > mx) mx = wgt;
          }
        }
      }
      const thr = mx * p.vfDropRel;
      for (let k = 0; k < cnt; k++) {
        if (rowVal[k] >= thr) { cols.push(rowCol[k]); vals.push(rowVal[k]); }
      }
      rowPtr[iw + 1] = cols.length;
    }
    return { rowPtr, col: Int32Array.from(cols), val: Float32Array.from(vals) };
  }

  // =========================================================================
  // ФОТООТЛИПАНИЕ: прямое ядро 1/(4 pi R^2) на огрублённой сетке (§7.5)
  //
  // Гельмгольц ЗДЕСЬ НЕ ГОДИТСЯ: при kappa -> 0 его функция Грина стремится к
  // 1/(4 pi R), а прозрачный поток имеет 1/(4 pi R^2). Это разные функции, и
  // никакой малой kappa их не примирить.
  //
  // Полная матрица на мелкой сетке — (nr*nz)^2 = 4.8 ГБ, исключено. Psi по
  // построению гладкая (оператор с ядром 1/R^2 сглаживающий), поэтому считаем
  // на огрублённой сетке photoCoarse = [32, 48] и интерполируем, добавляя
  // ближнюю поправку (иначе занижается поток внутри самого филамента).
  // =========================================================================
  _buildTransparentKernel() {
    const G = this.g, p = this.p;
    const ncr = p.photoCoarse[0], ncz = p.photoCoarse[1];
    const Rdom = G.rf[G.nr];
    const z0 = G.zf[this.jg0], z1 = G.zf[this.jg1 + 1];
    const dzc = (z1 - z0) / ncz;

    // Огрублённые ячейки: равномерно по r^2 (равные площади колец) — так каждая
    // огрублённая ячейка несёт одинаковый объём, и приосевая не вырождается.
    const crf = new Float64Array(ncr + 1);
    for (let i = 0; i <= ncr; i++) crf[i] = Rdom * Math.sqrt(i / ncr);
    const crc = new Float64Array(ncr), cdr = new Float64Array(ncr), cA = new Float64Array(ncr);
    for (let i = 0; i < ncr; i++) {
      crc[i] = 0.5 * (crf[i] + crf[i + 1]);
      cdr[i] = crf[i + 1] - crf[i];
      cA[i] = Math.PI * (crf[i + 1] * crf[i + 1] - crf[i] * crf[i]);
    }
    const czc = new Float64Array(ncz);
    for (let j = 0; j < ncz; j++) czc[j] = z0 + (j + 0.5) * dzc;

    const M = ncr * ncz;
    const K = new Float32Array(M * M);
    for (let ic = 0; ic < ncr; ic++) {
      // r' dr' по ячейке = Acell/(2pi); ядро даёт ещё 1/(2 sqrt(P1 P2)).
      const radW = cA[ic] / (2 * Math.PI) * dzc;
      for (let jc = 0; jc < ncz; jc++) {
        const src = ic * ncz + jc;
        for (let io = 0; io < ncr; io++) {
          for (let jo = 0; jo < ncz; jo++) {
            const obs = io * ncz + jo;
            let wgt;
            if (io === ic && jo === jc) {
              wgt = this._selfWeight(ic === 0, cA[ic] * dzc, cdr[ic], dzc);
            } else {
              wgt = radW * transparentKernel(crc[io], czc[jo], crc[ic], czc[jc]);
            }
            K[obs * M + src] = wgt;
          }
        }
      }
    }
    this.coarse = { ncr, ncz, crf, crc, cdr, cA, czc, dzc, z0, z1, K, M,
                    q: new Float64Array(M), psi: new Float64Array(M) };
    // Радиус ближней зоны — диагональ огрублённой ячейки (§7.5 рекомендует до двух;
    // множитель — параметр, стоимость поправки растёт как Rnear^2).
    const cdm = cdr[Math.floor(ncr / 2)];
    this.Rnear = p.photoNearFactor * Math.sqrt(cdm * cdm + dzc * dzc);

    // Индексные окна ближней зоны — считаются ОДИН РАЗ. Без них поправка сканирует
    // все ncr*ncz и все nr*ngz ячеек для каждой точки наблюдения (3.7e8 итераций
    // на обновление при nr=96) — это доминировало над всей стоимостью фотомодуля.
    this._nearWin = this._buildNearWindows();
  }

  _buildNearWindows() {
    const G = this.g, C = this.coarse, Rn = this.Rnear;
    const span = (centers, half, lo, hi, x) => {
      let a = lo, b = hi;
      while (a <= hi && centers[a] < x - Rn - half[a]) a++;
      while (b >= lo && centers[b] > x + Rn + half[b]) b--;
      return [a, b];
    };
    const zeroHalf = (n, v) => { const a = new Float64Array(n); a.fill(v); return a; };
    const fineR = [], fineZ = [], crsR = [], crsZ = [];
    const cHalfR = Float64Array.from(C.cdr), cHalfZ = zeroHalf(C.ncz, C.dzc);
    for (let i = 0; i < G.nr; i++) {
      fineR.push(span(G.rc, G.dr, 0, G.nr - 1, G.rc[i]));
      crsR.push(span(C.crc, cHalfR, 0, C.ncr - 1, G.rc[i]));
    }
    for (let j = this.jg0; j <= this.jg1; j++) {
      fineZ.push(span(G.zc, G.dz, this.jg0, this.jg1, G.zc[j]));
      crsZ.push(span(C.czc, cHalfZ, 0, C.ncz - 1, G.zc[j]));
    }
    return { fineR, fineZ, crsR, crsZ };
  }

  /**
   * Самоячейка прозрачного ядра. Точечное выражение сингулярно, 3D-интеграл
   * конечен — интегрируем по конечному объёму (PHOTO_PROCESSES.md §7.5):
   *   тороидальная ячейка (dr, dz << r), предел стержня:  Psi = C_self*q*sqrt(dr*dz)
   *   приосевая ячейка (r -> 0): предел стержня неприменим, ячейка компактна,
   *   шаровой эквивалент:        Psi = q*R_eq,  R_eq = (3V/(4 pi))^(1/3)
   * Возвращается вес при q (то есть Psi_self / q).
   */
  _selfWeight(isAxis, V, dr, dz) {
    if (isAxis) return Math.cbrt(3 * V / (4 * Math.PI));
    return this.p.Cself * Math.sqrt(dr * dz);
  }

  // =========================================================================
  // ОСНОВНОЙ ШАГ
  // =========================================================================
  /**
   * @param {object} f  {ionizRate: Float64Array(N) [м^-3 с^-1] — S_ion = k_E1*n_e*N,
   *                     powerDens: Float64Array(N) [Вт/м^3]    — P_dep = e*|Gamma_e . E|}
   * @param {object} o  {step: number} — номер шага для nPhotoSubstep
   * @returns {object}  ссылки на выходные массивы (переиспользуются, не копировать)
   */
  update(f = {}, o = {}) {
    const p = this.p;
    const anyOn = p.photoIonization || p.photoEmission || p.photoDetachment;
    // PH1/VF1: при всех трёх флагах false модуль обязан вернуть тождественные нули
    // и не выполнить ни одной операции — никакой «почти нулевой» инициализации.
    if (!anyOn) return this._out();

    const step = o.step ?? (++this._stepCount);
    const ns = Math.max(1, p.nPhotoSubstep | 0);
    if (this._lastUpdateStep >= 0 && (step % ns) !== 0) return this._out(); // поля заморожены
    this._lastUpdateStep = step;

    if (!this._built) this.build();

    const N = this.N;
    const S = f.ionizRate, P = f.powerDens;

    // 3a. Источники фотонов (PHOTO_PROCESSES.md §7.6).
    //     Окно A обслуживает ДВА потребителя — объёмную ионизацию и фотоэмиссию.
    const EA = this.p.hvA * EV, EB = this.p.hvB * EV, EC = this.p.hvC * EV;
    const qA = this._qA, qB = this._qB, qC = this._qC;
    const needA = p.photoIonization || p.photoEmission;
    if (needA) {
      if (S) { const k = p.etaGamma; for (let n = 0; n < N; n++) qA[n] = k * S[n]; }
      else qA.fill(0);
      this._zeroNonGas(qA);
    }
    if (p.photoEmission) {
      if (P) { const k = p.fVUV / EB; for (let n = 0; n < N; n++) qB[n] = k * P[n]; }
      else qB.fill(0);
      this._zeroNonGas(qB);
    }
    if (p.photoDetachment) {
      if (P) { const k = p.fDet / EC; for (let n = 0; n < N; n++) qC[n] = k * P[n]; }
      else qC.fill(0);
      this._zeroNonGas(qC);
    }

    if (p.photoIonization) this._doPhotoIonization(qA);
    else this.photoIonRate.fill(0);

    if (p.photoEmission) this._doPhotoEmission(qA, qB);
    else { this.photoEmitFluxL.fill(0); this.photoEmitFluxR.fill(0); }

    if (p.photoDetachment) this._doPhotoDetachment(qC);
    else { this.psi.fill(0); this.nuPdOm.fill(0); this.nuPdO2m.fill(0); this.nuPdO3m.fill(0); }

    if (p.photoCheck) this._assertFinite();
    return this._out();
  }

  _out() {
    return {
      photoIonRate: this.photoIonRate,
      psi: this.psi,
      nuPd: { Om: this.nuPdOm, O2m: this.nuPdO2m, O3m: this.nuPdO3m },
      photoEmitFluxL: this.photoEmitFluxL,
      photoEmitFluxR: this.photoEmitFluxR,
      gammaI: this.gammaI,
      FphOverFi: this.FphOverFi,
    };
  }

  _zeroNonGas(a) {
    const nz = this.g.nz, nr = this.g.nr;
    for (let i = 0; i < nr; i++) {
      const base = i * nz;
      for (let j = 0; j < this.jg0; j++) a[base + j] = 0;
      for (let j = this.jg1 + 1; j < nz; j++) a[base + j] = 0;
    }
  }

  // ---- 3b. Окно A: три группы Гельмгольца -------------------------------
  _doPhotoIonization(qA) {
    const p = this.p, N = this.N;
    const out = this.photoIonRate;
    out.fill(0);
    if (p.photoIonMode === 'off') return;
    if (p.photoIonMode === 'local') {
      // Отладочный режим: локальное замыкание S_pi = Phi_ion * q_A. Физически
      // неверен (теряет весь нелокальный перенос, ради которого канал и введён),
      // держится только для быстрых прогонов и сравнения.
      const k = p.PhiIon;
      for (let n = 0; n < N; n++) out[n] = k * qA[n];
      return;
    }
    if (!this.solver || typeof this.solver.solveHelmholtz !== 'function') {
      throw new Error('photo2d: photoIonization=true requires solver.solveHelmholtz (see sim2d/PHOTO_API.md)');
    }
    const scr = this._scr, Sj = this._Sj;
    for (let j = 0; j < 3; j++) {
      const Aj = this.Aj[j];
      for (let n = 0; n < N; n++) scr[n] = Aj * qA[n];
      const res = this.solver.solveHelmholtz(scr, this.lamA[j], Sj);
      const Sres = res || Sj;
      for (let n = 0; n < N; n++) out[n] += Sres[n];
    }
    const ph = p.PhiIon;
    for (let n = 0; n < N; n++) out[n] *= ph;
    this._zeroNonGas(out);
  }

  // ---- 3c. Фотоэмиссия: view-factor matvec ------------------------------
  _doPhotoEmission(qA, qB) {
    const nr = this.g.nr;
    const L = this.photoEmitFluxL, R = this.photoEmitFluxR;
    L.fill(0); R.fill(0);
    let qAm = null, qBm = null;
    if (this.vfMirror) { qAm = this._mirrorSrc(qA, this._qAm); qBm = this._mirrorSrc(qB, this._qBm); }
    for (let gi = 0; gi < this.emitGroups.length; gi++) {
      const grp = this.emitGroups[gi];
      if (!(grp.Y > 0)) continue;                 // VF7: строгий ноль ниже порога
      const isA = grp.window === 'A';
      const src = isA ? qA : qB;
      const wsrc = grp.wsrc;                      // доля источника, приходящаяся на группу
      this._csrAccum(this.vfL[gi], src, wsrc * grp.Y, L);
      if (this.vfMirror) this._csrAccum(this.vfL[gi], isA ? qAm : qBm, wsrc * grp.Y, R);
      else this._csrAccum(this.vfR[gi], src, wsrc * grp.Y, R);
    }
  }

  _csrAccum(m, q, scale, out) {
    if (!m) return;
    const rp = m.rowPtr, col = m.col, val = m.val;
    for (let i = 0; i < out.length; i++) {
      let s = 0;
      for (let k = rp[i]; k < rp[i + 1]; k++) s += val[k] * q[col[k]];
      out[i] += scale * s;
    }
  }

  // ---- 3d. Фотоотлипание: огрублённое ядро + ближняя поправка -----------
  _doPhotoDetachment(qC) {
    const p = this.p, G = this.g, N = this.N;
    const psi = this.psi;
    psi.fill(0);
    if (p.photoDetachMode === 'off') { this._applySigmas(); return; }

    if (p.photoDetachMode === 'uniform') {
      // Приближение равномерного засвета: Psi = G_slab * <q>.
      // ⚠ Для ОДИНОЧНОГО филамента 50–100 мкм это плохое приближение: теряется и
      // сильный профиль 1/R^2 вблизи филамента (занижение в 10–30 раз), и широкие
      // хвосты. Допустимо только для усреднённой модели множества микроразрядов.
      let num = 0, den = 0;
      for (let i = 0; i < G.nr; i++) {
        for (let j = this.jg0; j <= this.jg1; j++) {
          const n = i * G.nz + j;
          num += qC[n] * this.V[n]; den += this.V[n];
        }
      }
      const val = den > 0 ? p.GslabRel * num / den : 0;
      for (let i = 0; i < G.nr; i++)
        for (let j = this.jg0; j <= this.jg1; j++) psi[i * G.nz + j] = val;
      this._applySigmas();
      return;
    }

    const C = this.coarse;
    // 1) агрегируем мелкий источник на огрублённую сетку (по объёму)
    C.q.fill(0);
    const cw = new Float64Array(C.M);
    for (let i = 0; i < G.nr; i++) {
      const ic = this._coarseIr(G.rc[i]);
      for (let j = this.jg0; j <= this.jg1; j++) {
        const jc = Math.min(C.ncz - 1, Math.max(0, Math.floor((G.zc[j] - C.z0) / C.dzc)));
        const n = i * G.nz + j, c = ic * C.ncz + jc;
        C.q[c] += qC[n] * this.V[n];
        cw[c] += this.V[n];
      }
    }
    for (let c = 0; c < C.M; c++) C.q[c] = cw[c] > 0 ? C.q[c] / cw[c] : 0;

    // 2) matvec на огрублённой сетке
    const K = C.K, M = C.M, cpsi = C.psi;
    for (let obs = 0; obs < M; obs++) {
      let s = 0;
      const row = obs * M;
      for (let src = 0; src < M; src++) s += K[row + src] * C.q[src];
      cpsi[obs] = s;
    }

    // 3) билинейная интерполяция на мелкую сетку
    for (let i = 0; i < G.nr; i++) {
      const r = G.rc[i];
      for (let j = this.jg0; j <= this.jg1; j++) {
        psi[i * G.nz + j] = this._interpCoarse(r, G.zc[j]);
      }
    }

    // 4) ближняя поправка (§7.5): огрубление занижает поток внутри филамента.
    //    Psi = interp(Psi_coarse) - Psi_coarse_near + Psi_fine_near
    if (p.photoDetachNear) this._nearCorrection(qC, psi);

    this._applySigmas();
  }

  _coarseIr(r) {
    const C = this.coarse;
    // crf[i] = Rdom*sqrt(i/ncr) -> i = ncr*(r/Rdom)^2
    const Rdom = C.crf[C.ncr];
    let i = Math.floor(C.ncr * (r / Rdom) * (r / Rdom));
    if (i < 0) i = 0; if (i >= C.ncr) i = C.ncr - 1;
    return i;
  }

  _interpCoarse(r, z) {
    const C = this.coarse;
    // по r — линейно по узлам crc, по z — линейно по czc, с зажимом на краях
    let i = this._coarseIr(r);
    let i0 = i, i1 = i;
    if (r > C.crc[i] && i + 1 < C.ncr) i1 = i + 1;
    else if (r < C.crc[i] && i - 1 >= 0) i0 = i - 1;
    const tr = (i0 === i1) ? 0 : (r - C.crc[i0]) / (C.crc[i1] - C.crc[i0]);
    let j = Math.min(C.ncz - 1, Math.max(0, Math.floor((z - C.z0) / C.dzc - 0.5)));
    const j1 = Math.min(C.ncz - 1, j + 1);
    const tz = (j1 === j) ? 0 : (z - C.czc[j]) / C.dzc;
    const tzc = Math.min(1, Math.max(0, tz)), trc = Math.min(1, Math.max(0, tr));
    const a = C.psi[i0 * C.ncz + j] * (1 - tzc) + C.psi[i0 * C.ncz + j1] * tzc;
    const b = C.psi[i1 * C.ncz + j] * (1 - tzc) + C.psi[i1 * C.ncz + j1] * tzc;
    return a * (1 - trc) + b * trc;
  }

  _nearCorrection(qC, psi) {
    const G = this.g, C = this.coarse, Rn = this.Rnear, Rn2 = Rn * Rn;
    const nz = G.nz, W = this._nearWin;
    for (let i = 0; i < G.nr; i++) {
      const r = G.rc[i];
      const [ic0, ic1] = W.crsR[i], [if0, if1] = W.fineR[i];
      for (let j = this.jg0; j <= this.jg1; j++) {
        const z = G.zc[j];
        const [jc0, jc1] = W.crsZ[j - this.jg0], [jf0, jf1] = W.fineZ[j - this.jg0];
        let corr = 0;
        // вычитаем вклад огрублённых ячеек ближней зоны
        for (let ic = ic0; ic <= ic1; ic++) {
          const dr0 = Math.abs(C.crc[ic] - r);
          const radW = C.cA[ic] / (2 * Math.PI) * C.dzc;
          for (let jc = jc0; jc <= jc1; jc++) {
            const dz0 = C.czc[jc] - z;
            if (dr0 * dr0 + dz0 * dz0 > Rn2) continue;
            const inside = (r >= C.crf[ic] && r < C.crf[ic + 1] &&
                            z >= C.czc[jc] - 0.5 * C.dzc && z < C.czc[jc] + 0.5 * C.dzc);
            const w = inside
              ? this._selfWeight(ic === 0, C.cA[ic] * C.dzc, C.cdr[ic], C.dzc)
              : radW * transparentKernel(r, z, C.crc[ic], C.czc[jc]);
            corr -= w * C.q[ic * C.ncz + jc];
          }
        }
        // добавляем вклад мелких ячеек ближней зоны
        for (let ip = if0; ip <= if1; ip++) {
          const dr0 = Math.abs(G.rc[ip] - r);
          const radW = G.Acell[ip] / (2 * Math.PI);
          for (let jp = jf0; jp <= jf1; jp++) {
            const dz0 = G.zc[jp] - z;
            if (dr0 * dr0 + dz0 * dz0 > Rn2) continue;
            const w = (ip === i && jp === j)
              ? this._selfWeight(ip === 0, this.V[ip * nz + jp], G.dr[ip], G.dz[jp])
              : radW * G.dz[jp] * transparentKernel(r, z, G.rc[ip], G.zc[jp]);
            corr += w * qC[ip * nz + jp];
          }
        }
        const v = psi[i * nz + j] + corr;
        psi[i * nz + j] = v > 0 ? v : 0;   // PH8: поток не может быть отрицательным
      }
    }
  }

  _applySigmas() {
    const s = this.p.sigmaPd, N = this.N, psi = this.psi;
    const a = this.nuPdOm, b = this.nuPdO2m, c = this.nuPdO3m;
    for (let n = 0; n < N; n++) {
      const P = psi[n];
      a[n] = s[0] * P; b[n] = s[1] * P; c[n] = s[2] * P;
    }
  }

  // -------------------------------------------------------------------------
  // Прямое ядро окна C БЕЗ огрубления — эталон для тестов (PH3/PH4/VF9) и для
  // ближней поправки. O((nr*ngz)^2), использовать только на мелких сетках.
  // -------------------------------------------------------------------------
  psiExact(qC, out) {
    const G = this.g, nz = G.nz;
    const res = out || new Float64Array(this.N);
    res.fill(0);
    for (let i = 0; i < G.nr; i++) {
      const r = G.rc[i];
      for (let j = this.jg0; j <= this.jg1; j++) {
        const z = G.zc[j];
        let s = 0;
        for (let ip = 0; ip < G.nr; ip++) {
          const radW = G.Acell[ip] / (2 * Math.PI);
          for (let jp = this.jg0; jp <= this.jg1; jp++) {
            const q = qC[ip * nz + jp];
            if (q === 0) continue;
            const w = (ip === i && jp === j)
              ? this._selfWeight(ip === 0, this.V[ip * nz + jp], G.dr[ip], G.dz[jp])
              : radW * G.dz[jp] * transparentKernel(r, z, G.rc[ip], G.zc[jp]);
            s += w * q;
          }
        }
        res[i * nz + j] = s;
      }
    }
    return res;
  }

  /**
   * Суммарный темп фотоотлипания и его разложение по сортам.
   * S_e = sum_s n_s * nu_pd,s ; S_s = -n_s * nu_pd,s ;
   * нейтральные продукты: O- + hv -> O + e, O3- + hv -> O3 + e (сохранение вещества),
   * O2- + hv -> O2 + e (фон, не отслеживается). Заряд сохраняется тождественно
   * (каждый акт: один отрицательный ион -> один электрон), проверка A5 ERRATA.
   */
  detachRate(n) {
    const N = this.N, out = this.photoDetachRate;
    const Om = n.Om ?? n['O-'], O2m = n.O2m ?? n['O2-'], O3m = n.O3m ?? n['O3-'];
    for (let k = 0; k < N; k++) {
      let s = 0;
      if (Om)  s += Om[k]  * this.nuPdOm[k];
      if (O2m) s += O2m[k] * this.nuPdO2m[k];
      if (O3m) s += O3m[k] * this.nuPdO3m[k];
      out[k] = s;
    }
    return out;
  }

  _assertFinite() {
    // PH8: ни NaN, ни отрицательных потоков ни при каком разумном входе.
    const chk = (a, name) => {
      for (let n = 0; n < a.length; n++) {
        const v = a[n];
        if (!Number.isFinite(v)) throw new Error(`photo2d: non-finite value in ${name}[${n}] = ${v}`);
        if (v < 0) throw new Error(`photo2d: negative value in ${name}[${n}] = ${v}`);
      }
    };
    chk(this.photoIonRate, 'photoIonRate');
    chk(this.psi, 'psi');
    chk(this.photoEmitFluxL, 'photoEmitFluxL');
    chk(this.photoEmitFluxR, 'photoEmitFluxR');
  }
}

export default PhotoModule;
