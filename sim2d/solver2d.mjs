// solver2d.mjs — главный солвер 2D осесимметричного ДБР в чистом O2.
//
// Опорные документы (в порядке приоритета):
//   docs/ERRATA.md        — A1..A6 (блокеры), B1..B6, E1..E6, F (фото)
//   docs/NUMERICS_2D.md   — §4 (SG-потоки), §5 (шаг), §6 (стенки), §7 (ток), §10.2 (порядок)
//   docs/PHOTO_PROCESSES.md + sim2d/PHOTO_API.md — встраивание фотомодуля
//
// СОГЛАШЕНИЯ (жёстко):
//   * раскладка полей   field[i*nz + j], z — быстрый индекс;
//   * z-грань j лежит между ячейками j-1 и j:  EzF[i*(nz+1)+j], Gz[i*(nz+1)+j];
//   * r-грань i лежит между ячейками i-1 и i:  ErF[i*nz+j],     Gr[i*nz+j], i = 0..nr;
//   * поток положителен вдоль +z / +r;
//   * всё в СИ; комментарии по-русски, идентификаторы английские.
//
// ЧТО СДЕЛАНО ПО БЛОКЕРАМ ERRATA (кратко, подробности — у мест применения):
//   A1  знак sigma — целиком в poisson2d (buildRhs), здесь только накопление dSigma/dt;
//       тест S6 проверяет, что накопленный sigma ГАСИТ поле в зазоре.
//   A2  n_floor = 1 м^-3 ТОЛЬКО как защита от round-off; физическая затравка — явный член S_bg.
//   A3  ток Сато–Морроу по ФИЗИЧЕСКОМУ весовому полю E_L = 1/d_eff (точно одномерно, §7.1),
//       phi_B из полунеявной матрицы не используется нигде.
//   A5  источники — через reaction extents с ограничением ПО ЭКСТЕНТУ (не по плотностям).
//   A6  полный стек диэлектрик|газ|диэлектрик решается как единая r-z задача (poisson2d).
//   E1  |E| перед всеми LFA-таблицами (в physics2d absEN, плюс модуль здесь).
//   E2  dt_field = delta*max(|E|, E_scale)/max|dE/dt| + ограничение по фазе источника.
//   E3  dt_chem по ВСЕМ сортам.
//   E5  полунеявная поправка в полной форме (rhs содержит dt*div(J* - kappa*E^n)).

import {
  SPECIES_IDS, SPECIES_BY_ID, REACTIONS, stoichiometry,
  QE, KB, EPS0, ME, TD, N_REF,
  gasDensity, muE_N, DE_N, meanEnergy, kIoniz, setTuning,
} from './physics2d.mjs';
import { SeparableSolver } from './poisson2d.mjs';
import { PhotoModule } from './photo2d.mjs';

// ───────────────────────────────────────────────────────────── константы/утилиты

const IDX = Object.create(null);
SPECIES_IDS.forEach((id, i) => { IDX[id] = i; });
const NSP = SPECIES_IDS.length;          // 9 транспортируемых сортов
const CHARGED = SPECIES_IDS.filter((id) => SPECIES_BY_ID[id].charged).map((id) => IDX[id]);
const NCH = CHARGED.length;              // 6
const CE = CHARGED.indexOf(IDX.e);       // позиция электронов в списке заряженных

/**
 * B(x) = x/(exp(x)-1) — ядро Шарфеттера–Гуммеля (NUMERICS_2D §4.2).
 * Ряд Тейлора при |x| < 1e-4 (требование ТЗ), асимптотика при |x| > 40:
 * там результат неотличим от предельного в double, и заодно нет overflow exp.
 */
export function bern(x) {
  if (x > 40) return x * Math.exp(-x);
  if (x < -40) return -x;
  const ax = x < 0 ? -x : x;
  if (ax < 1e-4) return 1 - 0.5 * x + (x * x) / 12;
  return x / Math.expm1(x);
}

export const DEFAULTS = {
  // --- геометрия ---
  gapMM: 1.0, dielMM: 0.5, epsR: 9, radiusMM: 0.5,
  nr: 96, nz: 224, nzDiel: null,        // nzDiel = ячеек в КАЖДОМ барьере; null => nz/14
  betaR: 1.6, betaG: 2.63, qB: 1.3147,  // сгущение сетки (0/1 => равномерно)

  // --- источник и газ ---
  U0kV: 10, freqKHz: 10, pressureTorr: 760, tempK: 300,
  gammaIon: 0.02,                       // ЭФФЕКТИВНЫЙ гамма (фотоэмиссия свёрнута, ERRATA F2)

  // --- затравка ---
  seedBackground: 1e13,                 // фон n_e = n_O2p, м^-3
  seedSpotAmp: 1e18,                    // амплитуда гауссова пятна на оси (ERRATA B6)
  seedSpotSigmaUM: 40,                  // сигма пятна, мкм
  bgIonizRate: 1e13,                    // S_bg [м^-3 с^-1], ЯВНЫЙ фоновый источник пар (ERRATA A2)

  // --- фото ---
  photoModule: true,                    // создавать ли PhotoModule вообще
  photoIonization: true, photoEmission: true, photoDetachment: true,
  nPhotoSubstep: 3,
  photoParams: null,                    // прочие параметры PHOTO_DEFAULTS

  // --- шаг по времени ---
  dtMin: 1e-15, dtMax: 2e-8, dtStart: null,
  CFL: 0.4, diffSafety: 0.4, reacSafety: 0.2,
  dEtol: 0.02, Escale: 1e5,             // ERRATA E2
  // ERRATA E3 (см. развёрнутый комментарий у _limits): опорная плотность и фильтр
  // «активных» сортов, без которого сорт, рождающийся ИЗ НУЛЯ, пригвождает dt к dtMin.
  chemSafety: 0.1, nRefChem: 1e10, chemRelRef: 1e-6, chemActiveFrac: 1e-3,
  dtGrowth: 1.25, phaseSteps: 200,
  sigmaSafety: 0.05, sigmaRef: 1e-6,

  // --- численные ---
  nFloor: 1.0,                          // ERRATA A2: 1..1e2 м^-3, только round-off
  // детекторы разноса (см. _stepOk): «конечное и положительное» — НЕ признак
  // корректности, n_e = 8e64 м^-3 и E/N = 3e39 Тд удовлетворяют обоим
  ENmaxTrust: 5e3,                      // Тд; выше — заведомо нефизично (пробой ~150 Тд)
  growthMax: 30,                        // максимальный рост max n_s за ОДИН шаг
  negTol: 1e-6,                         // допустимая отрицательность ОТНОСИТЕЛЬНО максимума ТОГО ЖЕ сорта
  chemFrac: 0.9,                        // доля запаса реагента, доступная за шаг (лимитер extent)
  pcgRelTol: 1e-8, pcgMaxIter: 400,
  strictGauss: true,                    // точный Пуассон в начале шага (гарантия div(eps E) = rho)
  maxRetry: 8,                          // попыток отбраковки шага
  sigmaLeakTau: Infinity,               // §6.3, по умолчанию выключено
  // --- пристеночное ГУ (см. _wallCoeffs, ERRATA «STOP №2», DIVERGENCE_ANALYSIS §8.4/§9) ---
  //   'hagelaar' — ДЕФОЛТ. Потоко-ограничивающее ГУ Хагелаара: дрейфовый член входит
  //                со знаком (2a-1) и ВЫЧИТАЕТСЯ, когда поле отталкивает частицы,
  //                вплоть до полного запирания потока (kw = 0);
  //   'legacy'   — старая форма kw = max(0, v_n) + ¼v_th (поток на стенку НЕ зависит
  //                от поля). Сохранена только для сравнения: на ней получены числа
  //                VALIDATION.md и на ней же держится петля разгона §8.4.
  //   'thermal'  — устаревшее имя 'legacy' (принимается ради старых манифестов).
  // ПРЕДУПРЕЖДЕНИЕ О СТАРОМ ЗАМЕРЕ. В прежней редакции этого комментария стояло, что
  // «hagelaar не помогает» (авария на 0.6 % раньше). Тот замер сделан НЕПОЛНОЙ формой:
  // множитель (1-r)/(1+r) применялся сразу ко ВСЕЙ скобке, дрейфовый член не имел
  // своего множителя 1/(1+r), а линеаризации пристеночного потока в полунеявном
  // Пуассоне не было вовсе — то есть обратной связи по полю ВНУТРИ шага не было.
  // Здесь реализована полная форма 1D (src/solver.js::_wallCoeffs), проверенная там.
  wallBC: 'hagelaar',
  // reflE — коэффициент отражения ЭЛЕКТРОНОВ от диэлектрика, физический диапазон 0..0.2.
  // ВЫБОР МОДЕЛИ, НЕ ПОДГОНКА (обоснование дословно из 1D): r=0 (идеальное поглощение) —
  // идеализация, при которой в пристеночной ячейке остаётся непогасающая «плазма одной
  // ячейки» и атомарный кислород выходит в 140 раз плотнее самого газа ([O]/N = 1.4e2 —
  // заведомый артефакт). При r >= 0.05 артефакт исчезает ([O]/N = 1.0e-2) и
  // восстанавливается нормальная сеточная сходимость озона. Устойчивость от reflE
  // не зависит вообще; ёмкости не меняются в 4-м знаке. Разбор: DIVERGENCE_ANALYSIS §9.
  reflE: 0.05,
  reflI: 0.0,                           // коэффициент отражения ИОНОВ от диэлектрика
  wallReflect: null,                    // УСТАРЕЛО: если число — используется как reflE и reflI
  // Линеаризация пристеночного потока в полунеявном Пуассоне (см. _conductivityAndCurrent).
  // Выключается только для замера её вклада; на физику ГУ не влияет.
  wallLinearize: true,

  // --- отладочные режимы (для тестов транспорта) ---
  transportOnly: false,                 // без химии, стенок, sigma и фото
  frozenEz: null,                       // если число — поле заморожено, Пуассон не решается
  chemistry: true, wallFluxes: true, neutralTransport: true,

  d6MassConvention: 'neutral',          // ERRATA D6
};

// ───────────────────────────────────────────────────────────── класс

export class DBD2D {
  constructor(params = {}) {
    const p = { ...DEFAULTS, ...params };
    if (p.transportOnly) {
      p.chemistry = false; p.wallFluxes = false; p.photoModule = false; p.bgIonizRate = 0;
    }
    this.p = p;
    // конвенция массы для D6 — параметр физического модуля, он глобальный (документировано)
    setTuning({ d6MassConvention: p.d6MassConvention });

    const nzDiel = p.nzDiel ?? Math.max(3, Math.round(p.nz / 14));
    const ngap = p.nz - 2 * nzDiel;
    if (ngap < 4) throw new Error('DBD2D: слишком мало ячеек в газовом зазоре');

    const S = new SeparableSolver({
      nr: p.nr, nz1: nzDiel, ngap, nz2: nzDiel,
      radiusMM: p.radiusMM, dielMM: p.dielMM, gapMM: p.gapMM, epsR: p.epsR,
      betaR: p.betaR, betaG: p.betaG, qB: p.qB,
    });
    this.poisson = S;

    const nr = (this.nr = S.nr), nz = (this.nz = S.nz), n = (this.ncell = nr * nz);
    this.JG0 = S.JG0; this.JG1 = S.JG1;
    this.nzf = nz + 1;

    // метрика цилиндрической дивергенции (NUMERICS_2D §4.3)
    this.crIn = new Float64Array(nr);
    this.crOut = new Float64Array(nr);
    for (let i = 0; i < nr; i++) {
      const d = S.rf[i + 1] * S.rf[i + 1] - S.rf[i] * S.rf[i];
      this.crIn[i] = (2 * S.rf[i]) / d;      // на оси rf[0] = 0 => тождественный ноль (ERRATA B2)
      this.crOut[i] = (2 * S.rf[i + 1]) / d;
    }

    // --- газ ---
    this.N = gasDensity(p.pressureTorr, p.tempK);
    this.T = p.tempK;
    this.omega = 2 * Math.PI * p.freqKHz * 1e3;
    this.period = 1 / (p.freqKHz * 1e3);
    this.U0 = p.U0kV * 1e3;

    // --- сорта ---
    this.spec = SPECIES_IDS.map((id) => {
      const s = SPECIES_BY_ID[id];
      return {
        id, k: IDX[id], z: s.z, charged: !!s.charged, isE: id === 'e',
        mass: s.mass, muN: s.muN || 0, D300: s.D300 ?? 0, gammaWall: s.gammaWall ?? 0,
        mu: s.muN ? s.muN / this.N : 0,
        D: s.muN ? ((s.muN / this.N) * KB * p.tempK) / QE
                 : (s.D300 ? s.D300 * Math.pow(p.tempK / 300, 1.75) * (N_REF / this.N) : 0),
        vth: Math.sqrt((8 * KB * p.tempK) / (Math.PI * s.mass)),
      };
    });
    this.spE = this.spec[IDX.e];

    // --- реакции: предкомпиляция стехиометрии (ERRATA A5) ---
    this.rx = REACTIONS.map((r) => {
      let nBg = 0;
      const cnt = Object.create(null);
      for (const s of r.reagents) { if (s === 'O2') nBg++; else cnt[s] = (cnt[s] || 0) + 1; }
      const reag = Object.keys(cnt).map((id) => ({ k: IDX[id], c: cnt[id] }));
      const nu = stoichiometry(r);
      const prod = [];
      const cons = [];
      for (const id of Object.keys(nu)) {
        if (id === 'O2') continue;                 // фон не транспортируется
        const v = nu[id];
        if (v !== 0) prod.push({ k: IDX[id], v });
        if (v < 0) cons.push({ k: IDX[id], c: -v }); // ограничение — по ЧИСТОМУ расходу
      }
      return { id: r.id, rate: r.rate, nBg, reag, prod, cons };
    });
    // обратный индекс «сорт -> реакции, которые его РАСХОДУЮТ» (лимитер extent'ов)
    this.consBySpec = [];
    for (let s = 0; s < NSP; s++) {
      const list = [];
      this.rx.forEach((r, m) => { for (const c of r.cons) if (c.k === s) list.push({ m, c: c.c }); });
      this.consBySpec.push(list);
    }

    // --- поля ---
    this.n = {};
    this.nArr = [];
    this.nSave = [];
    for (const id of SPECIES_IDS) {
      const a = new Float64Array(n);
      this.n[id] = a; this.nArr.push(a); this.nSave.push(new Float64Array(n));
    }
    this.sigLo = new Float64Array(nr);
    this.sigHi = new Float64Array(nr);
    this.sigLoSave = new Float64Array(nr);
    this.sigHiSave = new Float64Array(nr);

    this.phi = new Float64Array(n);
    this.phiNew = new Float64Array(n);
    this.rho = new Float64Array(n);
    this.EzF = new Float64Array(nr * this.nzf);
    this.ErF = new Float64Array((nr + 1) * nz);
    this.EzF2 = new Float64Array(nr * this.nzf);
    this.ErF2 = new Float64Array((nr + 1) * nz);
    this.Ecell = new Float64Array(n);
    this.EcellPrev = new Float64Array(n);
    this.Ez = new Float64Array(n);
    this.Er = new Float64Array(n);
    this.EN = new Float64Array(n);
    this.ENzF = new Float64Array(nr * this.nzf);
    this.ENrF = new Float64Array((nr + 1) * nz);
    this.muEzF = new Float64Array(nr * this.nzf);
    this.muErF = new Float64Array((nr + 1) * nz);
    this.DezF = new Float64Array(nr * this.nzf);
    this.DerF = new Float64Array((nr + 1) * nz);
    this.kappaCell = new Float64Array(n);
    this.Fz = new Float64Array(nr * this.nzf);
    this.Fr = new Float64Array((nr + 1) * nz);
    this.ionizRate = new Float64Array(n);
    this.powerDens = new Float64Array(n);
    this.photoIonRate = new Float64Array(n);
    this.photoDetachRate = new Float64Array(n);
    this.photoEmitFluxL = new Float64Array(nr);
    this.photoEmitFluxR = new Float64Array(nr);
    this.Ssrc = [];
    for (let s = 0; s < NSP; s++) this.Ssrc.push(new Float64Array(n));

    this.Gz = [];
    this.Gr = [];
    for (let s = 0; s < NCH; s++) {
      this.Gz.push(new Float64Array(nr * this.nzf));
      this.Gr.push(new Float64Array((nr + 1) * nz));
    }
    this._gzTmp = new Float64Array(nr * this.nzf);
    this._grTmp = new Float64Array((nr + 1) * nz);
    this._dR = new Float64Array(this.rx.length);
    this._theta = new Float64Array(NSP);

    // маска газа по z (контракт state.gasMask)
    this.gasMask = new Uint8Array(nz);
    for (let j = 0; j < nz; j++) this.gasMask[j] = S.isGas[j];

    // --- фотомодуль ---
    this.photo = null;
    if (p.photoModule) {
      this.photo = new PhotoModule(
        { nr, nz, rf: S.rf, zf: S.zf, JG0: S.JG0, JG1: S.JG1,
          rc: S.rc, dr: S.dr, Acell: S.Acell, zc: S.zc, dz: S.dz },
        { photoIonization: p.photoIonization, photoEmission: p.photoEmission,
          photoDetachment: p.photoDetachment, nPhotoSubstep: p.nPhotoSubstep,
          gammaEff: p.gammaIon, ...(p.photoParams || {}) },
        S,
      );
    }
    // ERRATA F2 / PHOTO_API §3.2: ионный гамма берётся ТОЛЬКО из фотомодуля,
    // иначе вторичная эмиссия будет посчитана дважды.
    this.gammaI = this.photo ? this.photo.gammaI : p.gammaIon;
    // пристеночное ГУ (DEFAULTS.wallBC, реализация — _wallCoeffs)
    const bc = String(p.wallBC);
    this.wallLegacy = (bc === 'legacy' || bc === 'thermal');
    const clampR = (v) => (v > 0 ? (v < 0.95 ? v : 0.95) : 0);
    const rDep = Number.isFinite(p.wallReflect) ? p.wallReflect : null;
    this.reflE = clampR(rDep !== null ? rDep : p.reflE);
    this.reflI = clampR(rDep !== null ? rDep : p.reflI);
    // множители Хагелаара: th — при ¼v_th, fe — при дрейфовой скорости
    this.wThE = this.wallLegacy ? 0.25 : 0.25 * (1 - this.reflE) / (1 + this.reflE);
    this.wFeE = this.wallLegacy ? 1 : 1 / (1 + this.reflE);
    this.wThI = this.wallLegacy ? 0.25 : 0.25 * (1 - this.reflI) / (1 + this.reflI);
    this.wFeI = this.wallLegacy ? 1 : 1 / (1 + this.reflI);
    // коэффициенты ГУ по радиусу: [(side*NCH + c)*nr + i], side 0 = нижняя стенка
    this.kwWall = new Float64Array(2 * NCH * nr);
    this.linWall = new Float64Array(2 * NCH * nr);
    // диагностика потоков на стенку [м^-2 с^-1], индекс side*nr + i (заполняет _walls)
    this.wallGe = new Float64Array(2 * nr);      // электроны (с вычетом эмиссии)
    this.wallGpos = new Float64Array(2 * nr);    // положительные ионы
    this.wallGneg = new Float64Array(2 * nr);    // отрицательные ионы

    // --- диагностика/состояние ---
    this.t = 0;
    this.dt = p.dtStart ?? Math.min(p.dtMax, this.period / p.phaseSteps);
    this.stepIndex = 0;
    this.Q = 0;
    this.Icond = 0; this.Idisp = 0; this.Itot = 0;
    this.Ugap = 0; this.Uapp = 0;
    this.maxEN = 0; this.o3ppm = 0;
    this.maxdEdt = 0;
    this.nCG = 0; this.pcgResid = 0;
    this.limiter = 'init';
    this.clipCount = 0; this.qClip = 0;
    this.rejects = 0;
    this.dEff = S.dEff;
    this.ELz = 1 / S.dEff;                                  // ERRATA A3: физическое весовое поле
    this.Ccell = (EPS0 * Math.PI * S.R * S.R) / S.dEff;
    this.Vgas = 0;
    for (let i = 0; i < nr; i++) for (let j = S.JG0; j <= S.JG1; j++) this.Vgas += S.Acell[i] * S.dz[j];

    this._seed();
    this._updateRho();
    this._solveField(this.t, this.phi);
    this._fieldsFromPhi(this.phi, this.t, this.EzF, this.ErF, true);
    this.EcellPrev.set(this.Ecell);
    this.Uapp = this._U(this.t);
    this._diagUgap();
    this._diagScalars();          // state должен быть валиден ДО первого step()
  }

  // ───────────────────────────────────────────── затравка

  _seed() {
    const S = this.poisson, p = this.p, nz = this.nz;
    const ne = this.n.e, np = this.n.O2p;
    const sig = p.seedSpotSigmaUM * 1e-6;
    const z0 = S.zf[S.JG0] + 2 * sig;     // пятно у поверхности НИЖНЕГО диэлектрика (ERRATA B6)
    for (let i = 0; i < this.nr; i++) {
      const r = S.rc[i];
      for (let j = S.JG0; j <= S.JG1; j++) {
        const z = S.zc[j];
        let v = p.seedBackground;
        if (p.seedSpotAmp > 0) {
          const e1 = (r * r) / (2 * sig * sig) + ((z - z0) * (z - z0)) / (2 * sig * sig);
          if (e1 < 60) v += p.seedSpotAmp * Math.exp(-e1);
        }
        ne[i * nz + j] = v;
        np[i * nz + j] = v;                // квазинейтральная затравка
      }
    }
  }

  // ───────────────────────────────────────────── поле

  _U(t) { return this.U0 * Math.sin(this.omega * t); }

  _updateRho() {
    const nz = this.nz, rho = this.rho;
    rho.fill(0);
    for (const ci of CHARGED) {
      const sp = this.spec[ci], a = this.nArr[ci], q = QE * sp.z;
      for (let i = 0; i < this.nr; i++) {
        for (let j = this.JG0; j <= this.JG1; j++) { const k = i * nz + j; rho[k] += q * a[k]; }
      }
    }
  }

  /** Точный электростатический Пуассон (рабочая точка «начало шага»). */
  _solveField(t, out) {
    if (this.p.frozenEz !== null) { out.fill(0); return; }
    this.poisson.solvePoisson(this.rho, this.sigLo, this.sigHi, this._U(t), out);
  }

  /**
   * Поля на гранях + приведённое поле + (опц.) ячеечные величины.
   * ERRATA E1: везде |E|, таблицы LFA читаются по модулю поля НА ГРАНИ (§4.1).
   */
  _fieldsFromPhi(phi, t, EzF, ErF, cellToo) {
    const S = this.poisson, nr = this.nr, nz = this.nz, nzf = this.nzf;
    if (this.p.frozenEz !== null) {
      EzF.fill(this.p.frozenEz); ErF.fill(0);
    } else {
      S.computeEz(phi, this._U(t), EzF);
      S.computeEr(phi, ErF);
    }
    if (!cellToo) return;
    const N = this.N;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, tz = i * nzf;
      for (let j = 0; j < nz; j++) {
        const ez = 0.5 * (EzF[tz + j] + EzF[tz + j + 1]);
        const er = 0.5 * (ErF[o + j] + ErF[(i + 1) * nz + j]);
        this.Ez[o + j] = ez; this.Er[o + j] = er;
        const em = Math.hypot(ez, er);
        this.Ecell[o + j] = em;
        this.EN[o + j] = this.gasMask[j] ? em / (N * TD) : 0;
      }
    }
  }

  /** |E| и электронные коэффициенты на гранях (NUMERICS_2D §4.1). */
  _faceCoefs() {
    const nr = this.nr, nz = this.nz, nzf = this.nzf, N = this.N;
    const EzF = this.EzF, ErF = this.ErF;
    const j0 = this.JG0, j1 = this.JG1;
    // z-грани газа: j = j0..j1+1 (включая пристеночные)
    for (let i = 0; i < nr; i++) {
      const tz = i * nzf, o = i * nz, o1 = (i + 1) * nz;
      for (let j = j0; j <= j1 + 1; j++) {
        const jl = Math.max(j - 1, j0), jr = Math.min(j, j1);
        const erAt = 0.25 * (ErF[o + jl] + ErF[o1 + jl] + ErF[o + jr] + ErF[o1 + jr]);
        const em = Math.hypot(EzF[tz + j], erAt);
        const en = em / (N * TD);
        this.ENzF[tz + j] = en;
        this.muEzF[tz + j] = muE_N(en) / N;
        this.DezF[tz + j] = DE_N(en) / N;
      }
    }
    // r-грани газа: i = 1..nr-1
    for (let i = 1; i < nr; i++) {
      const o = i * nz, tzm = (i - 1) * nzf, tz = i * nzf;
      for (let j = j0; j <= j1; j++) {
        const ezAt = 0.25 * (EzF[tzm + j] + EzF[tzm + j + 1] + EzF[tz + j] + EzF[tz + j + 1]);
        const em = Math.hypot(ErF[o + j], ezAt);
        const en = em / (N * TD);
        this.ENrF[o + j] = en;
        this.muErF[o + j] = muE_N(en) / N;
        this.DerF[o + j] = DE_N(en) / N;
      }
    }
  }

  // ───────────────────────────────────────────── SG-потоки

  /**
   * SG-потоки одного сорта по обоим направлениям (NUMERICS_2D §4.3).
   * Цилиндрическая метрика в саму формулу НЕ входит — только в площади граней
   * при сборке дивергенции (иначе теряется консервативность).
   */
  _fluxSpecies(si, dens, Gz, Gr) {
    const S = this.poisson, nr = this.nr, nz = this.nz, nzf = this.nzf;
    const sp = this.spec[si], j0 = this.JG0, j1 = this.JG1;
    const isE = sp.isE, zc = sp.z;
    Gz.fill(0); Gr.fill(0);
    // --- z ---
    for (let i = 0; i < nr; i++) {
      const tz = i * nzf, o = i * nz;
      for (let j = j0 + 1; j <= j1; j++) {
        const h = S.hz[j];
        const mu = isE ? this.muEzF[tz + j] : sp.mu;
        const D = isE ? this.DezF[tz + j] : sp.D;
        const v = zc * mu * this.EzF[tz + j];
        Gz[tz + j] = sgFlux(dens[o + j - 1], dens[o + j], v, D, h);
      }
    }
    // --- r --- (Gr[0] и Gr[nr] тождественно нули: ось и зеркальная стенка)
    for (let i = 1; i < nr; i++) {
      const o = i * nz, om = (i - 1) * nz;
      const h = S.hr[i];
      for (let j = j0; j <= j1; j++) {
        const mu = isE ? this.muErF[o + j] : sp.mu;
        const D = isE ? this.DerF[o + j] : sp.D;
        const v = zc * mu * this.ErF[o + j];
        Gr[o + j] = sgFlux(dens[om + j], dens[o + j], v, D, h);
      }
    }
  }

  /** Консервативное обновление плотности по потокам (§4.3). */
  _advect(dens, Gz, Gr, dt) {
    const S = this.poisson, nr = this.nr, nz = this.nz, nzf = this.nzf;
    const j0 = this.JG0, j1 = this.JG1;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, tz = i * nzf, o1 = (i + 1) * nz;
      const cin = this.crIn[i], cout = this.crOut[i];
      for (let j = j0; j <= j1; j++) {
        dens[o + j] -= dt * ((Gz[tz + j + 1] - Gz[tz + j]) / S.dz[j]
                           + (cout * Gr[o1 + j] - cin * Gr[o + j]));
      }
    }
  }

  // ───────────────────────────────────────────── полунеявный Пуассон

  /** Проводимость плазмы по ячейкам [См/м] и явный ток на гранях. */
  _conductivityAndCurrent() {
    const nr = this.nr, nz = this.nz, nzf = this.nzf, N = this.N;
    const kap = this.kappaCell;
    kap.fill(0);
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      for (let j = this.JG0; j <= this.JG1; j++) {
        let s = 0;
        for (let c = 0; c < NCH; c++) {
          const sp = this.spec[CHARGED[c]];
          const mu = sp.isE ? muE_N(this.EN[o + j]) / N : sp.mu;
          s += mu * this.nArr[CHARGED[c]][o + j];
        }
        kap[o + j] = QE * s;
      }
    }
    // ── ЛИНЕАРИЗАЦИЯ ПРИСТЕНОЧНОГО ПОТОКА в полунеявном Пуассоне ─────────────
    // 1D-эталон: sc[fL] = QE*mu*lin*n на пристеночной ГРАНИ, откуда матрица берёт
    // g_g,eff = (eps0 + dt*kappa_wall)/(dx/2) (ERRATA A4). В 2D проводимость
    // передаётся в poisson2d ТОЛЬКО поячеечно (kappaCell), а транспортный
    // коэффициент грани газ/диэлектрик собирается гармонически из aCell соседей,
    // поэтому единственный доступный здесь эквивалент — задать kappa пристеночной
    // ГАЗОВОЙ ячейки из тех же множителей lin, что и сам поток:
    //     kappa_wall = QE * sum_s lin_s * z_s^2 * mu_s * n_s,
    // то есть d(J·n^)/dE_n на стенке. Смысл: когда ГУ Хагелаара ЗАПИРАЕТ поток
    // (lin = 0), поверхность внутри шага не отвечает на изменение поля вообще, и
    // полунеявная поправка не может «дозаряжать» sigma там, где заряд не идёт;
    // когда поток открыт, отклик равен ровно 1/(1+r) от дрейфового.
    // ЧЕСТНАЯ ОГОВОРКА: это же значение kappa участвует и в первой ВНУТРЕННЕЙ
    // грани (j0+1). Оценка ошибки: dt*kappa при n_e = 1e21, mu_e = 3e-2 и dt = 1e-13
    // равно 5e-16 Ф/м против eps0 = 8.9e-12, то есть вклад проводимости в эту грань
    // на 4 порядка ниже ёмкостного и подмена там не наблюдаема.
    if (!this.wallLegacy && this.p.wallLinearize && this.p.wallFluxes) {
      const nr2 = nr, j0 = this.JG0, j1 = this.JG1;
      for (let i = 0; i < nr2; i++) {
        const o = i * nz;
        for (let side = 0; side < 2; side++) {
          const jg = side === 0 ? j0 : j1;
          const kwB = (side * NCH) * nr2 + i;
          let s = 0;
          for (let c = 0; c < NCH; c++) {
            const sp = this.spec[CHARGED[c]];
            const mu = sp.isE ? muE_N(this.EN[o + jg]) / N : sp.mu;
            s += this.linWall[kwB + c * nr2] * sp.z * sp.z * mu * this.nArr[CHARGED[c]][o + jg];
          }
          kap[o + jg] = QE * s;
        }
      }
    }
    // J* на гранях (векторное, А/м^2) — понадобится для правой части (ERRATA E5)
    const Jz = this.Fz, Jr = this.Fr;
    Jz.fill(0); Jr.fill(0);
    for (let c = 0; c < NCH; c++) {
      const q = QE * this.spec[CHARGED[c]].z, gz = this.Gz[c], gr = this.Gr[c];
      for (let k = 0; k < nr * nzf; k++) Jz[k] += q * gz[k];
      for (let k = 0; k < (nr + 1) * nz; k++) Jr[k] += q * gr[k];
    }
  }

  /**
   * F = J* - kappa_face*E^n на гранях. kappa_face берётся ТОЧНО той же
   * гармонической сборкой, что использует матрица (poisson2d.buildFaceT),
   * иначе линеаризация в матрице и в правой части рассогласована.
   */
  _buildF(dt) {
    const S = this.poisson, nr = this.nr, nz = this.nz, nzf = this.nzf;
    const kap = this.kappaCell, eps = S.eps, dz = S.dz, dr = S.dr;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, tz = i * nzf;
      for (let j = 1; j < nz; j++) {
        const a1 = eps[j - 1] + dt * kap[o + j - 1], a2 = eps[j] + dt * kap[o + j];
        const aF = S.hz[j] / ((0.5 * dz[j - 1]) / a1 + (0.5 * dz[j]) / a2);
        const eF = S.hz[j] / ((0.5 * dz[j - 1]) / eps[j - 1] + (0.5 * dz[j]) / eps[j]);
        this.Fz[tz + j] -= ((aF - eF) / dt) * this.EzF[tz + j];
      }
      this.Fz[tz] = 0; this.Fz[tz + nz] = 0;   // на металле явного тока нет
    }
    for (let i = 1; i < nr; i++) {
      const o = i * nz, om = (i - 1) * nz;
      for (let j = 0; j < nz; j++) {
        const a1 = eps[j] + dt * kap[om + j], a2 = eps[j] + dt * kap[o + j];
        const aF = S.hr[i] / ((0.5 * dr[i - 1]) / a1 + (0.5 * dr[i]) / a2);
        this.Fr[o + j] -= ((aF - eps[j]) / dt) * this.ErF[o + j];
      }
    }
  }

  // ───────────────────────────────────────────── химия (ERRATA A5)

  /**
   * Один проход по реакциям в ячейке k.
   * @param apply false — только накопить чистые скорости S_s в this.Ssrc (для dt_chem),
   *              true  — посчитать extent'ы, ограничить их и применить.
   */
  _chemCell(k, dt, apply) {
    const EN = this.EN[k], N = this.N, T = this.T;
    const dR = this._dR, rx = this.rx, nArr = this.nArr;
    // Ssrc — ВАЛОВЫЙ оборот sum_j |nu_{s,j}|*rate_j, а не чистая скорость: ERRATA E3
    // прямо указывает, что |nu_i - nu_a| пропускает жёсткость при взаимной компенсации.
    if (!apply) for (let s = 0; s < NSP; s++) this.Ssrc[s][k] = 0;
    for (let m = 0; m < rx.length; m++) {
      const r = rx[m];
      let v = r.rate(EN, N, T);
      for (let q = 0; q < r.reag.length; q++) {
        const nn = nArr[r.reag[q].k][k];
        for (let c = 0; c < r.reag[q].c; c++) v *= nn;
      }
      if (r.nBg) v *= r.nBg === 1 ? N : N * N;
      if (!(v > 0)) v = 0;
      if (apply) dR[m] = v * dt;
      else for (let q = 0; q < r.prod.length; q++) this.Ssrc[r.prod[q].k][k] += Math.abs(r.prod[q].v) * v;
    }
    if (!apply) return;

    // --- ограничение ПО ЭКСТЕНТУ (а не по плотностям после факта) ---
    const theta = this._theta;
    theta.fill(1);
    const frac = this.p.chemFrac;
    for (let s = 0; s < NSP; s++) {
      const lst = this.consBySpec[s];
      if (lst.length === 0) continue;
      let cons = 0;
      for (let q = 0; q < lst.length; q++) cons += lst[q].c * dR[lst[q].m];
      if (cons > 0) {
        const avail = frac * nArr[s][k];
        if (cons > avail) theta[s] = avail / cons;
      }
    }
    for (let m = 0; m < rx.length; m++) {
      const r = rx[m];
      let th = 1;
      for (let q = 0; q < r.cons.length; q++) { const t2 = theta[r.cons[q].k]; if (t2 < th) th = t2; }
      if (th < 1) dR[m] *= th;
    }
    // применение ОДНОГО И ТОГО ЖЕ extent'а ко всем участникам => заряд сохраняется точно
    for (let m = 0; m < rx.length; m++) {
      const d = dR[m];
      if (d === 0) continue;
      const r = rx[m];
      for (let q = 0; q < r.prod.length; q++) nArr[r.prod[q].k][k] += r.prod[q].v * d;
    }
  }

  /** Скорости химии для ограничителя dt (ERRATA E3) — по ВСЕМ сортам. */
  _chemRates() {
    if (!this.p.chemistry) return;
    for (let i = 0; i < this.nr; i++) {
      for (let j = this.JG0; j <= this.JG1; j++) this._chemCell(i * this.nz + j, 0, false);
    }
  }

  /** Применение химии + фотоисточников + фонового источника (все — через extent'ы). */
  _chemApply(dt) {
    const p = this.p, nz = this.nz;
    const ne = this.n.e, nO2p = this.n.O2p;
    const nOm = this.n.Om, nO2m = this.n.O2m, nO3m = this.n.O3m;
    const nO = this.n.O, nO3 = this.n.O3;
    const pdOm = this.photo ? this.photo.nuPdOm : null;
    const pdO2m = this.photo ? this.photo.nuPdO2m : null;
    const pdO3m = this.photo ? this.photo.nuPdO3m : null;
    const Spi = this.photoIonRate;
    const Sbg = p.bgIonizRate;
    for (let i = 0; i < this.nr; i++) {
      for (let j = this.JG0; j <= this.JG1; j++) {
        const k = i * nz + j;
        if (p.chemistry) this._chemCell(k, dt, true);

        // --- фотоионизация: ОДИН extent на пару (e, O2+) => заряд сохраняется тождественно
        const dPi = Spi[k] * dt;
        if (dPi !== 0) { ne[k] += dPi; nO2p[k] += dPi; }

        // --- фоновый источник (ERRATA A2: физическая затравка — ЯВНЫЙ член, не floor)
        if (Sbg !== 0) { const d = Sbg * dt; ne[k] += d; nO2p[k] += d; }

        // --- фотоотлипание: extent ограничен запасом соответствующего аниона
        let pd = 0;
        if (pdOm) {
          let d = Math.min(nOm[k] * pdOm[k] * dt, p.chemFrac * nOm[k]);
          if (d > 0) { nOm[k] -= d; ne[k] += d; nO[k] += d; pd += d; }
          d = Math.min(nO2m[k] * pdO2m[k] * dt, p.chemFrac * nO2m[k]);
          if (d > 0) { nO2m[k] -= d; ne[k] += d; pd += d; }
          d = Math.min(nO3m[k] * pdO3m[k] * dt, p.chemFrac * nO3m[k]);
          if (d > 0) { nO3m[k] -= d; ne[k] += d; nO3[k] += d; pd += d; }
        }
        this.photoDetachRate[k] = dt > 0 ? pd / dt : 0;
      }
    }
  }

  // ───────────────────────────────────────────── пристеночное ГУ

  /**
   * ГУ ХАГЕЛААРА для потока заряженных частиц на диэлектрик — 2D-версия.
   *
   * Hagelaar G J M, de Hoog F J, Kroesen G M W, "Boundary conditions in fluid
   * models of gas discharges", Phys. Rev. E 62 (2000) 1452, ур. (18)-(21).
   * Эталон — src/solver.js::_wallCoeffs (1D), проверенный тестами 62/63 и
   * зарядовым балансом 33/33 (DIVERGENCE_ANALYSIS §9).
   *
   *     Gamma_s·n^ = max(0, (1-r)/(1+r)·¼·v_th,s·n_s + (2a-1)/(1+r)·|mu_s·E_n|·n_s )
   *
   *   a = 1, если дрейф направлен К стенке, a = 0 — если ОТ неё; r — коэффициент
   *   отражения (reflE / reflI). Так как (2a-1)·|mu_s·E_n| — это ровно ЗНАКОВАЯ
   *   проекция дрейфовой скорости на внешнюю нормаль, реализация записана в
   *   эквивалентной бесветвевой форме
   *
   *     kw_s = max(0, (1-r)/(1+r)·¼·v_th,s + (v_d·n^)/(1+r)),   Gamma_s·n^ = kw_s·n_s
   *
   *   где v_d = z_s·mu_s·E — вектор вдоль +z, n^ = -z на нижней стенке и +z на верхней.
   *
   * ОСЕСИММЕТРИЧНАЯ ГЕОМЕТРИЯ. Стенки — это границы газ/диэлектрик по z (грани
   * jf = JG0 и jf = JG1+1). Нормаль всюду ±z, но E_z на этих гранях СВОЁ в каждом
   * радиальном узле i, поэтому kw = kw(i, side, сорт): при зажигании филамента
   * приосевые узлы могут быть заперты, а периферийные — открыты.
   *
   * ЧТО ЭТО МЕНЯЕТ. Старое ГУ ('legacy') было kw = max(0, v_d·n^) + ¼v_th: при
   * отталкивающем поле дрейфовые «ворота» просто закрывались, и на стенку шёл
   * ПОЛНЫЙ тепловой поток, НЕ ЗАВИСЯЩИЙ ОТ ПОЛЯ. Дебаевская длина при n_e = 1e24,
   * T_e = 4 эВ равна 15 нм — сетка на порядки грубее, слой не разрешён, и это ГУ
   * замыкает петлю разгона sigma↓ → E_wall↑ → E/N↑ → k_ion↑ → n_e↑ (§8.4).
   * У Хагелаара дрейфовый член при a = 0 ВЫЧИТАЕТСЯ и поток падает до нуля при
   * |mu_s·E_n| >= (1-r)/(1+r)·¼·v_th. Это ЛИНЕЙНОЕ потоко-ограничение, а НЕ
   * больцмановский множитель exp(-e·Phi/kT) (тот пробовали — стало хуже, n_e → 2.45e25).
   *
   * Член -2/(1+r)·D·dn/dn СОЗНАТЕЛЬНО НЕ реализован: схема конечно-объёмная,
   * диффузия между пристеночной ячейкой и её соседом уже посчитана внутренней
   * гранью Шарфеттера–Гуммеля, добавление градиента ещё и на стенке дало бы
   * двойной счёт диффузии в первой ячейке.
   *
   * Побочно: на АНОДЕ положительные ионы отталкиваются и kw = 0 ТОЧНО. Старая форма
   * гнала туда поток ¼v_th·n — вторая ошибка того же ГУ, найденная в 1D.
   *
   * @param {Float64Array} EzF поле на z-гранях, по которому считается ГУ
   * @param {Float64Array} ErF поле на r-гранях (для |E| в LFA-таблицах электронов)
   * @param {boolean} refresh true — электронные mu и v_th пересчитать по ЭТОМУ полю
   *        (нужно, когда ГУ считается по E^{n+1}); false — взять готовые muEzF/ENzF.
   *
   * Заполняет kwWall (>=0, м/с) и linWall — множитель при (z·mu·n) для ЛИНЕАРИЗАЦИИ
   * пристеночного потока по полю: d(Gamma·n^)/dE_z = lin·z·mu·n. При запертом
   * (обрезанном нулём) потоке производная нулевая, при открытом lin = 1/(1+r).
   */
  _wallCoeffs(EzF, ErF, refresh) {
    const nr = this.nr, nz = this.nz, nzf = this.nzf, N = this.N;
    const j0 = this.JG0, j1 = this.JG1;
    const legacy = this.wallLegacy;
    const kw = this.kwWall, lin = this.linWall;
    for (let i = 0; i < nr; i++) {
      const tz = i * nzf, o = i * nz, o1 = (i + 1) * nz;
      for (let side = 0; side < 2; side++) {
        const jf = side === 0 ? j0 : j1 + 1;
        const jg = side === 0 ? j0 : j1;              // пристеночная ГАЗОВАЯ ячейка
        const nrm = side === 0 ? -1 : +1;             // внешняя нормаль (из газа)
        const Ew = EzF[tz + jf];
        let en, muEw;
        if (refresh) {
          // |E| на грани собирается ТОЧНО как в _faceCoefs (§4.1), иначе LFA-таблицы
          // будут прочитаны в другой точке и ГУ рассинхронизируется с транспортом
          const erAt = 0.5 * (ErF[o + jg] + ErF[o1 + jg]);
          en = Math.hypot(Ew, erAt) / (N * TD);
          muEw = muE_N(en) / N;
        } else {
          en = this.ENzF[tz + jf];
          muEw = this.muEzF[tz + jf];
        }
        const vthE = Math.sqrt((8 * QE * (2 / 3) * meanEnergy(en)) / (Math.PI * ME));
        const base = (side * NCH) * nr + i;
        for (let c = 0; c < NCH; c++) {
          const sp = this.spec[CHARGED[c]];
          const isE = sp.isE;
          const mu = isE ? muEw : sp.mu;
          const vth = isE ? vthE : sp.vth;
          const th = isE ? this.wThE : this.wThI;
          const fe = isE ? this.wFeE : this.wFeI;
          const vn = sp.z * mu * Ew * nrm;            // проекция дрейфа на внешнюю нормаль
          const idx = base + c * nr;
          if (legacy) {
            // СТАРАЯ форма, бит-в-бит как до переноса ГУ (порядок операций сохранён)
            kw[idx] = (vn > 0 ? Math.abs(mu * Ew) : 0) + 0.25 * vth;
            lin[idx] = vn > 0 ? 1 : 0;
          } else {
            const k = th * vth + fe * vn;
            if (!(k > 0)) { kw[idx] = 0; lin[idx] = 0; }
            else { kw[idx] = k; lin[idx] = fe; }
          }
        }
      }
    }
  }

  /**
   * Диагностика/тесты: коэффициент ГУ для сорта `id` на стенке `side` (0 — нижняя,
   * 1 — верхняя) у радиального узла i. Значения — от последнего вызова _wallCoeffs.
   * @returns {{kw:number, lin:number}} kw [м/с]: Gamma·n^ = kw*n; lin — множитель
   *          линеаризации d(Gamma·n^)/dE_n = lin*z*mu*n.
   */
  wallCoeff(side, id, i = 0) {
    const c = CHARGED.indexOf(IDX[id]);
    if (c < 0) throw new Error(`wallCoeff: сорт ${id} не заряжен`);
    const k = (side * NCH + c) * this.nr + i;
    return { kw: this.kwWall[k], lin: this.linWall[k] };
  }

  // ───────────────────────────────────────────── стенки (§6.1, §6.2)

  _walls(dt) {
    const S = this.poisson, nr = this.nr, nz = this.nz;
    const j0 = this.JG0, j1 = this.JG1;
    const gamI = this.gammaI;
    const Gs = this._Gwall || (this._Gwall = new Float64Array(NCH));
    const dSigLo = this._dSigLo || (this._dSigLo = new Float64Array(nr));
    const dSigHi = this._dSigHi || (this._dSigHi = new Float64Array(nr));
    // ГУ пересчитывается по НОВОМУ полю E^{n+1}: именно это даёт обратную связь по
    // полю ВНУТРИ шага (в 1D она была ключевой). В legacy — по E^n, как раньше,
    // чтобы старое поведение воспроизводилось бит-в-бит.
    if (this.wallLegacy) this._wallCoeffs(this.EzF, this.ErF, false);
    else this._wallCoeffs(this.EzF2, this.ErF2, true);
    let IcondWall = 0;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, A = S.Acell[i];
      for (let side = 0; side < 2; side++) {
        const jw = side === 0 ? j0 : j1;                 // пристеночная газовая ячейка
        const nrm = side === 0 ? -1 : +1;                // внешняя нормаль (из газа)
        const k = o + jw;
        // ── ПРИСТЕНОЧНЫЙ СТОК — НЕЯВНО (найдено валидацией, дефект B2) ──────────
        // Раньше сток применялся явно: n -= dt*G/dz. Скорость стока
        // nu = (a|mu E| + 0.25 v_th)/dz_wall для электронов достигает 1.1e10 c^-1
        // (постоянная времени 88 пс при dz_wall = 15.6 мкм) и НЕ ограничена CFL по
        // полю — тепловая часть от dt вообще не зависит. В пике тока сумма стока в
        // стенку и оттока через внутреннюю грань превышала содержимое ячейки,
        // n_e уходил в -2.5e18 м^-3, floor-клип возвращал их к +1 и инжектировал
        // 8.3 % перенесённого заряда (V12).
        // Неявная форма n_new = (n + dt*S/dz)/(1 + dt*nu) положительна при ЛЮБОМ dt,
        // а поток, записываемый в sigma и в ток, берётся СОГЛАСОВАННО:
        //   G = nu*dz*n_new - S   =>   n_new = n - dt*G/dz  тождественно,
        // поэтому заряд сохраняется машинно точно (тест S3/S10).
        const dzw = S.dz[jw];
        // Коэффициенты ГУ — из _wallCoeffs (Хагелаар при wallBC='hagelaar').
        // Их kw уже учитывает знак дрейфа, отражение r и запирание потока.
        const kwB = (side * NCH) * nr + i;
        let Gpos = 0, Gneg = 0;
        Gs.fill(0);
        for (let c = 0; c < NCH; c++) {
          const sp = this.spec[CHARGED[c]];
          if (sp.isE) continue;
          const nu = this.kwWall[kwB + c * nr] / dzw;
          const arr = this.nArr[CHARGED[c]];
          const nNew = arr[k] / (1 + dt * nu);
          arr[k] = nNew;
          const g = nu * dzw * nNew;
          Gs[c] = g;
          if (sp.z > 0) Gpos += g; else Gneg += g;
        }
        const Gpe = side === 0 ? this.photoEmitFluxL[i] : this.photoEmitFluxR[i];
        // ЕДИНСТВЕННОЕ правильное место для Gamma_pe — внутри Gamma_e (PHOTO_API §3.1).
        // Вторичная эмиссия берётся с ТОГО ЖЕ Gpos, что и поток ионов на стенку:
        // на аноде ГУ Хагелаара даёт kw_ion = 0 точно, поэтому и эмиссии там нет —
        // иначе электроны эмитировались бы с поверхности, на которую ионы не летят.
        const Se = gamI * Gpos + Gpe;                     // приход электронов со стенки
        const nuE = this.kwWall[kwB + CE * nr] / dzw;
        const neNew = (this.n.e[k] + (dt * Se) / dzw) / (1 + dt * nuE);
        this.n.e[k] = neNew;
        const Ge = nuE * dzw * neNew - Se;
        Gs[CE] = Ge;
        // диагностика (тест W3, условие плавающей стенки Gamma_e + Gamma_neg = Gamma_pos)
        this.wallGe[side * nr + i] = Ge;
        this.wallGpos[side * nr + i] = Gpos;
        this.wallGneg[side * nr + i] = Gneg;
        // накопление sigma (§6.2): уход электрона из диэлектрика = +заряд
        const dsigDt = QE * (Gpos - Ge - Gneg);
        if (side === 0) { this.sigLo[i] += dt * dsigDt; dSigLo[i] = dsigDt; }
        else { this.sigHi[i] += dt * dsigDt; dSigHi[i] = dsigDt; }
        // вклад пристеночной грани в ток Сато–Морроу: Jz = nrm*e*(Gpos - Ge - Gneg)
        IcondWall += nrm * QE * (Gpos - Ge - Gneg) * A * 0.5 * S.dz[jw];
      }
      // нейтралы: прилипание к стенке
      if (this.p.neutralTransport) {
        for (let s = NCH; s < NSP; s++) {
          const sp = this.spec[s], a = this.nArr[s];
          if (!(sp.gammaWall > 0)) continue;
          for (let side = 0; side < 2; side++) {
            const jw = side === 0 ? j0 : j1, kk = o + jw;
            const nu = (0.25 * sp.gammaWall * sp.vth) / S.dz[jw];
            a[kk] /= (1 + dt * nu);                        // тот же неявный сток
          }
        }
      }
    }
    this._IcondWall = IcondWall * this.ELz;
    if (Number.isFinite(this.p.sigmaLeakTau)) {
      const f = Math.exp(-dt / this.p.sigmaLeakTau);
      for (let i = 0; i < nr; i++) { this.sigLo[i] *= f; this.sigHi[i] *= f; }
    }
  }

  // ───────────────────────────────────────────── шаг по времени

  _limits() {
    const S = this.poisson, nr = this.nr, nz = this.nz, N = this.N, p = this.p;
    let invDrift = 0, invDiff = 0, invReac = 0, invChem = 0;
    let nMax = 0;
    for (let s = 0; s < NSP; s++) {
      const a = this.nArr[s];
      for (let i = 0; i < nr; i++) for (let j = this.JG0; j <= this.JG1; j++) {
        const v = a[i * nz + j]; if (v > nMax) nMax = v;
      }
    }
    // ─────────────────────────────────────────────────────────────────────────
    // ERRATA E3, честная оговорка. Буквальная форма dt <= 0.1*min_s (n_s+n_ref)/|S_s|
    // по ВСЕМ сортам НЕРАБОТОСПОСОБНА для сорта, рождающегося ИЗ НУЛЯ: при n_s = 0
    // и n_ref = 1e10 трёхтельное прилипание (S(O2-) = 1.1e26 м^-3с^-1 при n_e = 1e17)
    // даёт dt = 8.3e-17 с, и критерий не самовосстанавливается (за шаг dt плотность
    // набирается как S*dt, то есть предел растёт линейно вместе с dt). Замерено:
    // O4+ 7.0e-17 с, O2- 8.3e-17 с при t = 0 — dt намертво прибит к dtMin.
    // Поэтому критерий предъявляется только к «активным» сортам
    // (n_s > chemActiveFrac * max_s n_s); положительность ОСТАЛЬНЫХ гарантирована
    // не шагом, а ограничением extent'ов (ERRATA A5), которое работает при любом dt.
    // Восстановить буквальную форму: chemActiveFrac = 0.
    // ─────────────────────────────────────────────────────────────────────────
    const nRef = Math.max(p.nRefChem, p.chemRelRef * nMax);
    const nActive = p.chemActiveFrac * nMax;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, dr = S.dr[i];
      for (let j = this.JG0; j <= this.JG1; j++) {
        const k = o + j, dz = S.dz[j];
        const en = this.EN[k];
        const muE = muE_N(en) / N, DE = DE_N(en) / N;
        for (let c = 0; c < NCH; c++) {
          const sp = this.spec[CHARGED[c]];
          const mu = sp.isE ? muE : sp.mu, D = sp.isE ? DE : sp.D;
          const cfl = (Math.abs(mu * this.Ez[k]) / dz) + (Math.abs(mu * this.Er[k]) / dr);
          if (cfl > invDrift) invDrift = cfl;
          const df = (2 * D) / (dz * dz) + (2 * D) / (dr * dr);
          if (df > invDiff) invDiff = df;
        }
        for (let s = NCH; s < NSP; s++) {
          const D = this.spec[s].D;
          const df = (2 * D) / (dz * dz) + (2 * D) / (dr * dr);
          if (df > invDiff) invDiff = df;
        }
        // реакционный: рост лавины за шаг (ТЗ / §5.1)
        if (p.chemistry) {
          const nu = kIoniz(en) * N;
          if (nu > invReac) invReac = nu;
          for (let s = 0; s < NSP; s++) {
            const ns = this.nArr[s][k];
            if (ns < nActive) continue;                     // см. комментарий ниже
            const Sv = this.Ssrc[s][k];
            if (Sv > 0) {
              const q = Sv / (ns + nRef);                   // ERRATA E3
              if (q > invChem) invChem = q;
            }
          }
        }
      }
    }
    // Пристеночный сток отдельного ограничения по dt НЕ требует: он применяется
    // НЕЯВНО (см. _walls, дефект B2) и положителен при любом dt.
    const lim = {};
    lim.drift = invDrift > 0 ? p.CFL / invDrift : Infinity;
    lim.diff = invDiff > 0 ? p.diffSafety / invDiff : Infinity;
    lim.reac = invReac > 0 ? p.reacSafety / invReac : Infinity;
    lim.chem = invChem > 0 ? p.chemSafety / invChem : Infinity;
    // ERRATA E2: защита от нуля на старте синусоиды
    let maxE = 0;
    for (let i = 0; i < nr; i++) for (let j = 0; j < nz; j++) {
      const v = this.Ecell[i * nz + j]; if (v > maxE) maxE = v;
    }
    lim.field = this.maxdEdt > 0
      ? (p.dEtol * Math.max(maxE, p.Escale)) / this.maxdEdt : Infinity;
    lim.phase = this.period / p.phaseSteps;
    let maxSig = 0, maxDSig = 0;
    for (let i = 0; i < nr; i++) {
      maxSig = Math.max(maxSig, Math.abs(this.sigLo[i]), Math.abs(this.sigHi[i]));
      if (this._dSigLo) {
        maxDSig = Math.max(maxDSig, Math.abs(this._dSigLo[i]), Math.abs(this._dSigHi[i]));
      }
    }
    lim.surf = maxDSig > 0 ? (p.sigmaSafety * (maxSig + p.sigmaRef)) / maxDSig : Infinity;
    return lim;
  }

  _chooseDt() {
    const p = this.p;
    const lim = this._limits();
    let best = Infinity, name = 'none';
    for (const key of Object.keys(lim)) if (lim[key] < best) { best = lim[key]; name = key; }
    let dt = Math.min(best, this.dt * p.dtGrowth, p.dtMax);
    if (dt < p.dtMin) dt = p.dtMin;
    this.limiter = name;
    return dt;
  }

  // ───────────────────────────────────────────── диагностика

  _diagUgap() {
    const S = this.poisson, nz = this.nz, nr = this.nr;
    const dz = S.dz, eps = S.eps;
    let acc = 0, area = 0;
    const gdLo = eps[S.JS_LO] / (0.5 * dz[S.JS_LO]), ggLo = eps[S.JG0] / (0.5 * dz[S.JG0]);
    const ggHi = eps[S.JG1] / (0.5 * dz[S.JG1]), gdHi = eps[S.JS_HI] / (0.5 * dz[S.JS_HI]);
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const pLo = (gdLo * this.phi[o + S.JS_LO] + ggLo * this.phi[o + S.JG0] + this.sigLo[i]) / (gdLo + ggLo);
      const pHi = (ggHi * this.phi[o + S.JG1] + gdHi * this.phi[o + S.JS_HI] + this.sigHi[i]) / (ggHi + gdHi);
      acc += (pLo - pHi) * S.Acell[i];
      area += S.Acell[i];
    }
    this.Ugap = acc / area;
    return this.Ugap;
  }

  /** Мгновенное напряжение на зазоре при ЗАДАННЫХ rho/sigma и приложенном U (для теста S6). */
  gapVoltage(U, useSigma = true) {
    const zero = this._zeroSig || (this._zeroSig = new Float64Array(this.nr));
    const tmp = this._phiTmp || (this._phiTmp = new Float64Array(this.ncell));
    this.poisson.solvePoisson(this.rho, useSigma ? this.sigLo : zero,
      useSigma ? this.sigHi : zero, U, tmp);
    const S = this.poisson, nz = this.nz, dz = S.dz, eps = S.eps;
    const gdLo = eps[S.JS_LO] / (0.5 * dz[S.JS_LO]), ggLo = eps[S.JG0] / (0.5 * dz[S.JG0]);
    const ggHi = eps[S.JG1] / (0.5 * dz[S.JG1]), gdHi = eps[S.JS_HI] / (0.5 * dz[S.JS_HI]);
    let acc = 0, area = 0;
    for (let i = 0; i < this.nr; i++) {
      const o = i * nz;
      const sl = useSigma ? this.sigLo[i] : 0, sh = useSigma ? this.sigHi[i] : 0;
      const pLo = (gdLo * tmp[o + S.JS_LO] + ggLo * tmp[o + S.JG0] + sl) / (gdLo + ggLo);
      const pHi = (ggHi * tmp[o + S.JG1] + gdHi * tmp[o + S.JS_HI] + sh) / (ggHi + gdHi);
      acc += (pLo - pHi) * S.Acell[i]; area += S.Acell[i];
    }
    return acc / area;
  }

  /** Полный заряд: объёмный (газ) + поверхностный. Для теста сохранения заряда S3. */
  charge() {
    const S = this.poisson, nz = this.nz;
    let qv = 0, qs = 0;
    for (const ci of CHARGED) {
      const a = this.nArr[ci], q = QE * this.spec[ci].z;
      for (let i = 0; i < this.nr; i++) for (let j = this.JG0; j <= this.JG1; j++) {
        qv += q * a[i * nz + j] * S.Acell[i] * S.dz[j];
      }
    }
    for (let i = 0; i < this.nr; i++) qs += (this.sigLo[i] + this.sigHi[i]) * S.Acell[i];
    return { volume: qv, surface: qs, total: qv + qs };
  }

  /** Ток проводимости по Сато–Морроу (ERRATA A3, §7.2): E_L = 1/d_eff, точно 1D. */
  _sato() {
    const S = this.poisson, nr = this.nr, nz = this.nz, nzf = this.nzf;
    let acc = 0;
    for (let i = 0; i < nr; i++) {
      const tz = i * nzf, A = S.Acell[i];
      for (let j = this.JG0 + 1; j <= this.JG1; j++) {
        let jz = 0;
        for (let c = 0; c < NCH; c++) jz += QE * this.spec[CHARGED[c]].z * this.Gz[c][tz + j];
        acc += jz * A * S.hz[j];
      }
    }
    return acc * this.ELz + (this._IcondWall || 0);
  }

  _floorClip() {
    const S = this.poisson, nz = this.nz, fl = this.p.nFloor;
    for (let s = 0; s < NSP; s++) {
      const a = this.nArr[s], q = QE * this.spec[s].z;
      for (let i = 0; i < this.nr; i++) for (let j = this.JG0; j <= this.JG1; j++) {
        const k = i * nz + j;
        if (a[k] < fl) {
          this.qClip += q * (fl - a[k]) * S.Acell[i] * S.dz[j];
          a[k] = fl; this.clipCount++;
        }
      }
    }
  }

  _save() {
    for (let s = 0; s < NSP; s++) this.nSave[s].set(this.nArr[s]);
    this.sigLoSave.set(this.sigLo); this.sigHiSave.set(this.sigHi);
    this._qClipSave = this.qClip; this._clipSave = this.clipCount;
  }

  _restore() {
    for (let s = 0; s < NSP; s++) this.nArr[s].set(this.nSave[s]);
    this.sigLo.set(this.sigLoSave); this.sigHi.set(this.sigHiSave);
    this.qClip = this._qClipSave; this.clipCount = this._clipSave;
  }

  // ───────────────────────────────────────────── ОДИН ШАГ

  /**
   * Порядок операций — NUMERICS_2D §10.2.
   * @returns {number} фактически сделанный dt
   */
  step() {
    const p = this.p;
    this.Uapp = this._U(this.t);

    // 1. поле в начале шага. strictGauss: точный Пуассон => div(eps E) = rho машинно,
    //    дрейф закона Гаусса за 1e5 шагов исключён по построению (а не «держится на 1e-12»).
    if (p.strictGauss || this.stepIndex === 0) {
      this._updateRho();
      this._solveField(this.t, this.phi);
    }
    this._fieldsFromPhi(this.phi, this.t, this.EzF, this.ErF, true);
    this._diagUgap();
    this._faceCoefs();

    // 2. явные SG-потоки заряженных сортов
    for (let c = 0; c < NCH; c++) this._fluxSpecies(CHARGED[c], this.nArr[CHARGED[c]], this.Gz[c], this.Gr[c]);

    // 3. диагностические поля для фотомодуля (ERRATA E1: ionizRate по |E|)
    this._ionizAndPower();

    // 4. фотомодуль (до источников, после поля — PHOTO_API §2)
    this._photoUpdate();

    // 5. скорости химии для ограничителя dt
    this._chemRates();

    let dt = this._chooseDt();
    this._save();

    for (let attempt = 0; ; attempt++) {
      // потоки пересчитываются на КАЖДОЙ попытке: после отбраковки они уже содержат
      // коррекцию по старому dt, и переиспользовать их нельзя
      if (attempt > 0) {
        for (let c = 0; c < NCH; c++) {
          this._fluxSpecies(CHARGED[c], this.nArr[CHARGED[c]], this.Gz[c], this.Gr[c]);
        }
      }
      // 6. проводимость + J* + правая часть полунеявной задачи (ERRATA E5, полная форма).
      //    Пристеночные коэффициенты по E^n нужны ДО Пуассона: из них строится
      //    линеаризация пристеночного потока (см. _conductivityAndCurrent).
      this._wallCoeffs(this.EzF, this.ErF, false);
      this._conductivityAndCurrent();
      this._buildF(dt);

      // 7. полунеявный Пуассон -> phi^{n+1}
      if (p.frozenEz !== null) {
        this.EzF2.set(this.EzF); this.ErF2.set(this.ErF); this.nCG = 0;
      } else {
        this._updateRho();
        let res = null;
        try {
          res = this.poisson.solveSemiImplicit({
            rho: this.rho, sigLo: this.sigLo, sigHi: this.sigHi, U: this._U(this.t + dt),
            dt, kappaCell: this.kappaCell, Fz: this.Fz, Fr: this.Fr,
            phiGuess: this.phi, relTol: p.pcgRelTol, maxIter: p.pcgMaxIter,
          });
        } catch (e) {
          // отказ факторизации (например «прогонка неустойчива») — это ОТБРАКОВКА шага,
          // а не повод уронить прогон: причина всегда выше по течению, в разносе решения
          this.rejectReason = `poisson: ${e.message}`;
          res = null;
        }
        if (res === null) {
          if (this._retryOrThrow(dt, attempt)) { dt *= 0.4; continue; }
        }
        this.phiNew.set(res.phi);
        this.nCG = res.iters; this.pcgResid = res.resid;
        this._fieldsFromPhi(this.phiNew, this.t + dt, this.EzF2, this.ErF2, false);
      }

      // 8. коррекция потоков по новому полю (§4.4) — до обновления плотностей
      this._correctFluxes();

      // 9. транспорт заряженных
      for (let c = 0; c < NCH; c++) this._advect(this.nArr[CHARGED[c]], this.Gz[c], this.Gr[c], dt);
      // 10. транспорт нейтралов (только диффузия)
      if (p.neutralTransport) {
        for (let s = NCH; s < NSP; s++) {
          this._fluxSpecies(s, this.nArr[s], this._gzTmp, this._grTmp);
          this._advect(this.nArr[s], this._gzTmp, this._grTmp, dt);
        }
      }

      // 11. ток проводимости (по скорректированным потокам)
      this._IcondWall = 0;
      // 12. стенки: вторичная + фотоэмиссия, тепловой сток, накопление sigma
      if (p.wallFluxes) this._walls(dt);

      // 13. источники через reaction extents
      this._chemApply(dt);

      // --- отбраковка шага ---
      // ВАЖНО (найдено валидацией, см. docs/VALIDATION.md §«Что было не так», дефект B1):
      // приёмка ОБЯЗАНА идти ДО floor-клипа. Раньше _floorClip() стоял выше и затирал
      // отрицательные плотности единицей, поэтому детектор `negative` в _stepOk был
      // мёртвым кодом: nMin никогда не был < 0. Замерено на прогоне 2x76: в пике тока
      // электроны уходили в n_e = -2.5e18 м^-3 в пристеночной ячейке, клип возвращал их
      // к +1 м^-3, и это инжектировало 2.7e-11 Кл — 8.3 % от перенесённого за импульс
      // заряда (V12 требует < 1e-4).
      if (this._stepOk()) { this._floorClip(); break; }
      this._retryOrThrow(dt, attempt);
      dt *= 0.4;
    }

    // 15. диагностика
    this.Icond = this._sato();
    this.Idisp = this.Ccell * this.U0 * this.omega * Math.cos(this.omega * (this.t + 0.5 * dt));
    this.Itot = this.Icond + this.Idisp;
    this.Q += this.Itot * dt;

    // dE/dt для ERRATA E2 (по модулю поля в ячейках)
    if (p.frozenEz === null) {
      let mx = 0;
      const nz = this.nz;
      for (let i = 0; i < this.nr; i++) for (let j = 0; j < nz; j++) {
        const k = i * nz + j;
        const d = Math.abs(this.Ecell[k] - this.EcellPrev[k]);
        if (d > mx) mx = d;
      }
      this.maxdEdt = mx / dt;
      this.EcellPrev.set(this.Ecell);
    }

    if (!p.strictGauss) this.phi.set(this.phiNew);
    this.t += dt;
    this.dt = dt;
    this.stepIndex++;
    this._diagScalars();
    return dt;
  }

  /**
   * Критерий приёмки шага. ВАЖНО: проверки «конечное» и «неотрицательное»
   * НЕДОСТАТОЧНО. Замерено на прогоне 12x48: после потери разрешения катодного слоя
   * решение уходило в n_e = 8.3e64 м^-3 и E/N = 3.0e39 Тд, оставаясь при этом
   * конечным и положительным, и шаг принимался — разнос был замечен только когда
   * прогонка Пуассона выдала max|c'| = Infinity, 6000 шагов спустя.
   * Поэтому здесь ещё два детектора: абсолютный потолок E/N и максимальный рост
   * плотности за шаг.
   */
  _stepOk() {
    let nMax = 0;
    const nz = this.nz;
    // Отрицательность проверяется ПОСОРТНО и относительно максимума ЭТОГО ЖЕ сорта.
    // Глобальный максимум по всем сортам (его задают нейтралы O/O3, на 2-3 порядка
    // больше электронов) делал допуск бессмысленно мягким.
    for (let s = 0; s < NSP; s++) {
      const a = this.nArr[s];
      let mx = 0, mn = 0;
      for (let i = 0; i < this.nr; i++) for (let j = this.JG0; j <= this.JG1; j++) {
        const v = a[i * nz + j];
        if (!Number.isFinite(v)) { this.rejectReason = 'nonfinite'; return false; }
        if (v > mx) mx = v;
        if (v < mn) mn = v;
      }
      if (mn < -this.p.negTol * mx) {
        this.rejectReason = `n(${SPECIES_IDS[s]}) = ${mn.toExponential(2)} при max ${mx.toExponential(2)}`;
        return false;
      }
      if (mx > nMax) nMax = mx;
    }
    // E/N берём по ПОЛЮ КОНЦА ШАГА (EzF2/ErF2), а не по this.EN: тот посчитан от phi^n
    // и заметил бы разнос на шаг позже.
    let eMax = 0;
    const nzf = this.nzf;
    for (let i = 0; i < this.nr; i++) {
      const tz = i * nzf, o = i * nz;
      for (let j = this.JG0; j <= this.JG1 + 1; j++) {
        const v = Math.abs(this.EzF2[tz + j]);
        if (!Number.isFinite(v)) { this.rejectReason = 'E nonfinite'; return false; }
        if (v > eMax) eMax = v;
      }
      for (let j = this.JG0; j <= this.JG1; j++) {
        const v = Math.abs(this.ErF2[o + j]);
        if (!Number.isFinite(v)) { this.rejectReason = 'E nonfinite'; return false; }
        if (v > eMax) eMax = v;
      }
    }
    const enMax = eMax / (this.N * TD);
    if (enMax > this.p.ENmaxTrust) { this.rejectReason = `E/N = ${enMax.toExponential(2)} Тд`; return false; }
    if (this._nMaxPrev > 0 && nMax > this.p.growthMax * this._nMaxPrev) {
      this.rejectReason = `рост max n за шаг x${(nMax / this._nMaxPrev).toExponential(2)}`;
      return false;
    }
    this._nMaxPrev = nMax;
    this.rejectReason = null;
    return true;
  }

  /**
   * Откат шага. Если дробить dt больше некуда — падаем ГРОМКО.
   * Молча продолжать с dt = dtMin нельзя: именно так разнос доезжает до n_e = 1e64
   * и «результат» выглядит числами, а не ошибкой.
   */
  _retryOrThrow(dt, attempt) {
    this._restore();
    this.rejects++;
    if (attempt >= this.p.maxRetry || dt * 0.4 < this.p.dtMin) {
      throw new Error(
        `DBD2D: шаг ${this.stepIndex} (t = ${this.t.toExponential(6)} с) не проходит приёмку `
        + `даже при dt = ${dt.toExponential(3)} с. Причина: ${this.rejectReason}. `
        + 'Типовой источник — разгон объёмного заряда в ПРИСТЕНОЧНОЙ ячейке. '
        + 'ВНИМАНИЕ: прежняя формулировка «лечится сгущением z-сетки» ОПРОВЕРГНУТА измерением '
        + '(docs/DIVERGENCE_ANALYSIS §8.3: dz_wall от 15.6 до 0.184 мкм — срыв всегда, момент '
        + 'сдвигается на 6 %). Действующий диагноз (§8.4): тепловой поток электронов ¼v_th·n_e '
        + 'на КАТОДНЫЙ барьер не подавлялся отталкивающим полем, ток на барьер менял знак, '
        + 'и петля σ↓→E↑→ионизация↑ уходила в разнос. Лечится ГУ Хагелаара '
        + '(wallBC: \'hagelaar\', дефолт) — дрейфовый член входит со знаком (2a−1). '
        + 'Если авария всё же случилась при wallBC = \'hagelaar\', это НЕ известный дефект: '
        + 'сгущение сетки и уменьшение dt здесь не помогают, нужна диагностика. '
        + 'Состояние откачено к началу шага и пригодно для чекпойнта.',
      );
    }
    return true;
  }

  _correctFluxes() {
    if (this.p.frozenEz !== null) return;
    const nr = this.nr, nz = this.nz, nzf = this.nzf;
    const j0 = this.JG0, j1 = this.JG1;
    for (let c = 0; c < NCH; c++) {
      const sp = this.spec[CHARGED[c]], gz = this.Gz[c], gr = this.Gr[c];
      const dens = this.nArr[CHARGED[c]];
      for (let i = 0; i < nr; i++) {
        const tz = i * nzf, o = i * nz;
        for (let j = j0 + 1; j <= j1; j++) {
          const mu = sp.isE ? this.muEzF[tz + j] : sp.mu;
          const dE = this.EzF2[tz + j] - this.EzF[tz + j];
          const v = sp.z * mu * this.EzF[tz + j];
          const nup = v >= 0 ? dens[o + j - 1] : dens[o + j];
          gz[tz + j] += sp.z * mu * nup * dE;
        }
      }
      for (let i = 1; i < nr; i++) {
        const o = i * nz, om = (i - 1) * nz;
        for (let j = j0; j <= j1; j++) {
          const mu = sp.isE ? this.muErF[o + j] : sp.mu;
          const dE = this.ErF2[o + j] - this.ErF[o + j];
          const v = sp.z * mu * this.ErF[o + j];
          const nup = v >= 0 ? dens[om + j] : dens[o + j];
          gr[o + j] += sp.z * mu * nup * dE;
        }
      }
    }
  }

  _ionizAndPower() {
    const nr = this.nr, nz = this.nz, nzf = this.nzf, N = this.N;
    const ne = this.n.e;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, tz = i * nzf, o1 = (i + 1) * nz;
      for (let j = this.JG0; j <= this.JG1; j++) {
        const k = o + j;
        this.ionizRate[k] = kIoniz(this.EN[k]) * ne[k] * N;
        const gz = 0.5 * (this.Gz[0][tz + j] + this.Gz[0][tz + j + 1]);
        const gr = 0.5 * (this.Gr[0][o + j] + this.Gr[0][o1 + j]);
        this.powerDens[k] = QE * Math.abs(gz * this.Ez[k] + gr * this.Er[k]);
      }
    }
  }

  _photoUpdate() {
    if (!this.photo) {
      this.photoIonRate.fill(0);
      this.photoEmitFluxL.fill(0); this.photoEmitFluxR.fill(0);
      return;
    }
    const out = this.photo.update(
      { ionizRate: this.ionizRate, powerDens: this.powerDens },
      { step: this.stepIndex },
    );
    this.photoIonRate.set(out.photoIonRate);
    this.photoEmitFluxL.set(out.photoEmitFluxL);
    this.photoEmitFluxR.set(out.photoEmitFluxR);
  }

  _diagScalars() {
    const nz = this.nz, S = this.poisson;
    let mx = 0, o3 = 0, neMax = 0;
    for (let i = 0; i < this.nr; i++) {
      for (let j = this.JG0; j <= this.JG1; j++) {
        const k = i * nz + j;
        if (this.EN[k] > mx) mx = this.EN[k];
        if (this.n.e[k] > neMax) neMax = this.n.e[k];
        o3 += this.n.O3[k] * S.Acell[i] * S.dz[j];
      }
    }
    this.maxEN = mx;
    this.neMax = neMax;
    this.o3ppm = ((o3 / this.Vgas) / this.N) * 1e6;
    let sm = 0, pl = 0, pr = 0;
    for (let i = 0; i < this.nr; i++) {
      sm = Math.max(sm, Math.abs(this.sigLo[i]), Math.abs(this.sigHi[i]));
      pl += this.photoEmitFluxL[i] * S.Acell[i];
      pr += this.photoEmitFluxR[i] * S.Acell[i];
    }
    this.sigmaMax = sm;
    this.photoEmitTotalL = pl;
    this.photoEmitTotalR = pr;
  }

  // ───────────────────────────────────────────── контракт состояния

  get state() {
    const S = this.poisson;
    return {
      t: this.t, dt: this.dt, step: this.stepIndex,
      r: S.rc, z: S.zc, rf: S.rf, zf: S.zf, gasMask: this.gasMask,
      JG0: this.JG0, JG1: this.JG1,
      n: this.n,
      E: this.Ecell, Er: this.Er, Ez: this.Ez, phi: this.phi, rho: this.rho,
      ionizRate: this.ionizRate, EN: this.EN,
      photoIonRate: this.photoIonRate, photoDetachRate: this.photoDetachRate,
      photoEmitFluxL: this.photoEmitFluxL, photoEmitFluxR: this.photoEmitFluxR,
      sigmaL: this.sigLo, sigmaR: this.sigHi,
      Uapp: this.Uapp, Ugap: this.Ugap,
      Icond: this.Icond, Idisp: this.Idisp, Itot: this.Itot, Q: this.Q,
      o3ppm: this.o3ppm, maxEN: this.maxEN,
      // поля, которые ждёт recorder.mjs (SERIES_NAMES + быстрый neMax вместо скана n.e)
      neMax: this.neMax, sigmaMax: this.sigmaMax,
      photoEmitTotalL: this.photoEmitTotalL, photoEmitTotalR: this.photoEmitTotalR,
      nCG: this.nCG, limiter: this.limiter, rejects: this.rejects, clipCount: this.clipCount,
    };
  }
}

/** SG-поток через грань (постоянные v, D на отрезке между центрами). */
function sgFlux(nHere, nThere, v, D, h) {
  if (D < 1e-30) return v > 0 ? v * nHere : v * nThere;
  const X = (v * h) / D;
  if (X > 40 || X < -40) return v > 0 ? v * nHere : v * nThere;   // точный предел SG = апвинд
  const bp = bern(X), bm = bp + X;                                 // B(-x) = B(x) + x
  return (D / h) * (nHere * bm - nThere * bp);
}

export default DBD2D;
