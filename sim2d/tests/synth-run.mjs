// synth-run.mjs — генератор ПРАВДОПОДОБНОГО синтетического прогона в НАСТОЯЩЕМ формате.
//
// Цель: разблокировать разработку плеера ДО готовности реального солвера.
// Это НЕ физика: здесь нет ни Пуассона, ни переноса, ни химии. Это форма сигнала,
// снятая с ожидаемого поведения (NUMERICS_2D §9, PHYSICS.md, ERRATA §C), одетая
// в тот же контейнер, что и настоящий прогон.
//
//   node sim2d/tests/synth-run.mjs [--periods=1] [--out=data] [--run=run-synth] [--profile=default]
//
// Что воспроизводится:
//   * сетка профиля `default` (nr=96, nz=224 = 16/192/16) по генераторам NUMERICS_2D §1.3;
//   * гауссова головка стримера, идущая через 1 мм зазора за ~5 нс (v ~ 2e5 м/с);
//   * импульс тока ~100 нс с САМОГАШЕНИЕМ поверхностным зарядом (ERRATA A1: σ экранирует);
//   * растекание σ(r) по диэлектрику от 80 до 250 мкм за импульс;
//   * замкнутый параллелограмм Лиссажу Q(U) (по антисимметрии полупериодов);
//   * монотонная наработка O3 и накопление O3⁻ между импульсами;
//   * фотоионизация как размытая копия свечения (ℓ = 14.8/90/238 мкм -> два масштаба),
//     фотоотлипание от O3⁻, интегральная фотоэмиссия на обе поверхности.
//
// ЧЕСТНОЕ ОГРАНИЧЕНИЕ (то же, что у настоящей модели): осесимметричная геометрия
// описывает только ЦЕНТРАЛЬНЫЙ канал; азимутальные моды отсутствуют принципиально.

import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';
import { createRecorders } from '../recorder.mjs';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, '..', '..');

// ---- CLI -------------------------------------------------------------------
const argv = Object.fromEntries(process.argv.slice(2).map((a) => {
  const m = /^--([^=]+)=?(.*)$/.exec(a);
  return m ? [m[1], m[2] === '' ? true : m[2]] : [a, true];
}));
const PERIODS = Number(argv.periods ?? 2);
const RUN = String(argv.run ?? 'run-synth');
const OUT = path.resolve(ROOT, String(argv.out ?? 'data'));
const PROFILE = String(argv.profile ?? 'default');

// ---- геометрия и сетка (NUMERICS_2D §1.3) ----------------------------------
const PROFILES = {
  demo:     { nr: 80,  nz1: 14, ngap: 148, nz2: 14 },
  default:  { nr: 96,  nz1: 16, ngap: 192, nz2: 16 },
  accurate: { nr: 160, nz1: 24, ngap: 272, nz2: 24 },
};
const P = PROFILES[PROFILE] || PROFILES.default;
const R = 0.5e-3, D1 = 0.5e-3, LG = 1.0e-3, D2 = 0.5e-3, LZ = D1 + LG + D2;
const NR = P.nr, NZ1 = P.nz1, NGAP = P.ngap, NZ2 = P.nz2, NZ = NZ1 + NGAP + NZ2;
const JG0 = NZ1, JG1 = NZ1 + NGAP - 1;

function makeRGrid(nr, Rmax, betaR = 1.6) {
  const rf = new Float64Array(nr + 1);
  const den = Math.exp(betaR) - 1;
  for (let k = 0; k <= nr; k++) rf[k] = Rmax * Math.expm1(betaR * k / nr) / den;
  return rf;
}
function makeZGrid(nz1, ngap, nz2, d1, Lg, d2, betaG = 2.63, qB = 1.3147) {
  const zf = [0];
  const dz0 = d1 * (qB - 1) / (Math.pow(qB, nz1) - 1);
  for (let k = nz1 - 1; k >= 0; k--) zf.push(zf[zf.length - 1] + dz0 * Math.pow(qB, k));
  const th = Math.tanh(betaG / 2), z0 = zf[zf.length - 1];
  for (let k = 1; k <= ngap; k++) {
    const s = 2 * k / ngap - 1;
    zf.push(z0 + 0.5 * Lg * (1 + Math.tanh(betaG * s / 2) / th));
  }
  for (let k = 0; k < nz2; k++) zf.push(zf[zf.length - 1] + dz0 * Math.pow(qB, k));
  return Float64Array.from(zf);
}

const rf = makeRGrid(NR, R), zf = makeZGrid(NZ1, NGAP, NZ2, D1, LG, D2);
const rc = new Float64Array(NR), zc = new Float64Array(NZ);
for (let i = 0; i < NR; i++) rc[i] = 0.5 * (rf[i] + rf[i + 1]);
for (let j = 0; j < NZ; j++) zc[j] = 0.5 * (zf[j] + zf[j + 1]);
const ZS_LO = zf[JG0], ZS_HI = zf[JG1 + 1];      // поверхности диэлектриков

// ---- параметры прогона ------------------------------------------------------
const U0 = 10e3, FREQ = 10e3, TPER = 1 / FREQ;
const EPS0 = 8.854187817e-12, EPSR = 9;
const AREA = Math.PI * R * R;
const C_G = EPS0 * AREA / LG;                       // 0.695 пФ при R=0.5мм
const C_D = EPS0 * EPSR * AREA / (D1 + D2);         // последовательные барьеры
const C_CELL = C_G * C_D / (C_G + C_D);
const U_B = 3.76e3;                                 // ERRATA C, V2: пробой на газе
const U_HOLD = 0.80 * U_B;                          // напряжение горения
const N_GAS = 2.45e25;                              // м^-3 при 1 атм, 300 K
const V_GAP = AREA * LG;                            // объём зазора, м^3
const O3_YIELD = 3.5e17;                            // молекул O3 на Дж (100 г/кВт*ч)
const GAMMA_SW = (C_G + C_D) / 30e-9;               // проводимость канала: tau = 30 нс

// ---- состояние ---------------------------------------------------------------
const N = NR * NZ;
const f64 = () => new Float64Array(N);
const st = {
  t: 0, dt: 0, r: rc, z: zc, rf, zf,
  gasMask: new Uint8Array(NZ),
  n: { e: f64(), O2p: f64(), O4p: f64(), Om: f64(), O2m: f64(), O3m: f64(), O: f64(), O3: f64(), O2a: f64() },
  E: f64(), Er: f64(), Ez: f64(), phi: f64(), rho: f64(), ionizRate: f64(), EN: f64(),
  photoIonRate: f64(), photoDetachRate: f64(),
  photoEmitFluxL: new Float64Array(NR), photoEmitFluxR: new Float64Array(NR),
  sigmaL: new Float64Array(NR), sigmaR: new Float64Array(NR),
  Uapp: 0, Ugap: 0, Icond: 0, Idisp: 0, Itot: 0, Q: 0,
  maxEN: 0, o3ppm: 0, sigmaMax: 0, photoEmitTotalL: 0, photoEmitTotalR: 0,
  neMax: 0,
};
for (let j = JG0; j <= JG1; j++) st.gasMask[j] = 1;

// динамика разряда
let qSurf = 0;            // перенесённый заряд на диэлектрике, Кл (знак = знак экранирования)
let burning = false;      // горит ли микроразряд
let tIgn = 0;             // момент зажигания
let headZ = ZS_LO;        // положение головки
let headDir = +1;
let o3Total = 0;          // м^-3, средняя по зазору
let neBulk = 1e13;        // «фон» электронов
let nO3mBulk = 1e17;      // O3- между импульсами
let UappPrev = 0;
let tEnd_ = -1;           // момент конца последнего микроразряда
let nPulses = 0;

const T_PULSE = 100e-9;   // длительность микроразряда
const T_DEAD = 300e-9;    // мёртвое время до следующего зажигания
const HEAD_V = 1.0e5;     // м/с — скорость стримера (нижний край 1e5..1e6 м/с)
// Пятно ПОВЕРХНОСТНОГО ЗАРЯДА заметно шире светящегося канала: заряд разносит
// поверхностный стример вдоль диэлектрика. Из условия самогашения
// sigma ~ eps0*E_br*(C_g+C_d)/C_g даёт sigma_max ~ 3e-4 Кл/м^2 (ERRATA V8: E_mem 5.2e6 В/м),
// то есть радиус следа 0.2...0.4 мм при переносимом заряде ~50 пКл.
const SIG_R0 = 2.0e-4, SIG_R1 = 4.0e-4;  // след σ(r): 200 -> 400 мкм за импульс
const HEAD_SR = 6e-5, HEAD_SZ = 4e-5;   // радиус/длина головки

// ---- заполнение 2D-полей (только на кадрах!) --------------------------------
function fillFields() {
  const Ugap = st.Ugap;
  const Ebg = Math.abs(Ugap) / LG;
  const active = burning ? 1 : 0;
  const tau = st.t - tIgn;
  // амплитуда головки спадает после прохода зазора (остаётся канал + катодный слой)
  const crossing = burning && Math.abs(headZ - (headDir > 0 ? ZS_HI : ZS_LO)) > 1e-6;
  const headAmp = crossing ? 1 : Math.exp(-Math.max(0, tau - LG / HEAD_V) / 4e-8);
  const chanAmp = burning ? Math.exp(-Math.max(0, tau - LG / HEAD_V) / 6e-8) : 0;

  const neHead = 1e20 * headAmp * active;
  const neChan = 1e19 * chanAmp * active;
  const sionHead = 3e29 * headAmp * active;
  const rhoAmp = 2.0e-3 * headAmp * active;
  const Ehead = 2.0e7 * headAmp * active;

  let neMax = 0, enMax = 0;
  const zCath = headDir > 0 ? ZS_LO : ZS_HI;
  for (let i = 0; i < NR; i++) {
    const ri = rc[i];
    const gr = Math.exp(-(((ri - spotR) / HEAD_SR) ** 2));
    const grCh = Math.exp(-(((ri - spotR) / (1.4 * HEAD_SR)) ** 2));
    const grPh = Math.exp(-(((ri - spotR) / 1.6e-4) ** 2));            // ядро фотоионизации, ℓ~90 мкм
    const grHalo = Math.exp(-(((ri - spotR) / 4.0e-4) ** 2));          // ореол, ℓ~238 мкм
    for (let j = 0; j < NZ; j++) {
      const k = i * NZ + j;
      if (j < JG0 || j > JG1) {                    // в барьерах — только поле
        st.n.e[k] = 0; st.ionizRate[k] = 0; st.rho[k] = 0;
        st.n.O3m[k] = 0; st.n.O3[k] = 0;
        st.photoIonRate[k] = 0; st.photoDetachRate[k] = 0;
        st.E[k] = Ebg / EPSR; st.EN[k] = 0;
        continue;
      }
      const zj = zc[j];
      const dzh = (zj - headZ) / HEAD_SZ;
      const gz = Math.exp(-dzh * dzh);
      const behind = headDir > 0 ? (zj < headZ) : (zj > headZ);
      const gzCh = behind ? Math.exp(-Math.abs(zj - headZ) / 3e-4) : 0;

      // плотность электронов: головка + послесвечение канала + фон
      const ne = neBulk + neHead * gr * gz + neChan * grCh * gzCh;
      st.n.e[k] = ne;
      if (ne > neMax) neMax = ne;

      // свечение (скорость ионизации) — резко локализовано в головке
      const sion = sionHead * gr * gz * gz + (burning ? 1e24 * grCh * gzCh : 0);
      st.ionizRate[k] = sion;

      // объёмный заряд: ЗНАКОПЕРЕМЕННЫЙ. Положительный заряд в головке катодо-направленного
      // стримера (знак НЕ зависит от направления хода) + отрицательный слой в следе канала,
      // где электроны отстают от ионов. Именно этот знакопеременный профиль обязан пережить
      // квантование без сдвига нуля.
      st.rho[k] = rhoAmp * (gr * gz - 0.45 * grCh * gzCh);

      // поле: фон + усиление в головке + катодный слой у «мгновенного катода»
      const cath = Math.exp(-Math.abs(zj - zCath) / 2e-5) * (burning ? 1 : 0);
      const E = Ebg + Ehead * gr * gz + 0.8 * Ehead * cath * grCh;
      st.E[k] = E;
      const en = E / (N_GAS * 1e-21);
      st.EN[k] = en;
      if (en > enMax) enMax = en;

      // фотоионизация: размытая копия свечения, eta_gamma*Phi_ion ~ 1e-3,
      // два масштаба (узкое ядро + широкий ореол ~10 % потока)
      const gzPh = Math.exp(-(((zj - headZ) / 1.2e-4) ** 2));
      const gzHalo = Math.exp(-(((zj - headZ) / 4.0e-4) ** 2));
      st.photoIonRate[k] = 1e-3 * sionHead * (0.9 * grPh * gzPh + 0.1 * grHalo * gzHalo);

      // O3⁻ — доминирующий отрицательный ион между импульсами (ERRATA F2)
      const o3m = nO3mBulk * (1 + 2 * grCh * gzCh);
      st.n.O3m[k] = o3m;
      // фотоотлипание: ν_pd ~ 10 c^-1 при светящемся импульсе (PHOTO_PROCESSES §5)
      st.photoDetachRate[k] = (burning ? 10 : 1e-3) * o3m * (0.5 + 0.5 * grHalo);

      // озон: накопленный фон + приращение в канале
      st.n.O3[k] = o3Total * (1 + 0.6 * grCh * Math.exp(-(((zj - 0.5 * (ZS_LO + ZS_HI)) / 3e-4) ** 2))) + 1e18;
    }
  }
  st.neMax = neMax;
  st.maxEN = enMax;
  // фотоэмиссия с поверхностей (интегральный поток, VF-критерии PHOTO_PROCESSES)
  let fl = 0, fr = 0;
  for (let i = 0; i < NR; i++) {
    const w = Math.exp(-(((rc[i] - spotR) / 2.2e-4) ** 2));
    const A = Math.PI * (rf[i + 1] ** 2 - rf[i] ** 2);
    st.photoEmitFluxL[i] = 1e-4 * sionHead * 1e-6 * w * (headDir > 0 ? 1 : 0.3);
    st.photoEmitFluxR[i] = 1e-4 * sionHead * 1e-6 * w * (headDir > 0 ? 0.3 : 1);
    fl += st.photoEmitFluxL[i] * A; fr += st.photoEmitFluxR[i] * A;
  }
  st.photoEmitTotalL = fl; st.photoEmitTotalR = fr;
}

// σ(r): НАКОПИТЕЛЬНАЯ модель. Каждый шаг ток проводимости осаждает dq в пятне
// текущего радиуса, затем заряд растекается по поверхности (поверхностная диффузия).
// Так интеграл σ по площади ТОЧНО равен перенесённому заряду, а профиль живёт:
// узкое пятно во время импульса -> расплывание к следующему зажиганию.
// Поверхностная диффузия на Al2O3 МАЛА (заряд не растекается — потому память ДБР и
// локальна). Значит последовательные микроразряды не могут садиться в одно место:
// уже осаждённый заряд экранирует поле локально, и следующий пробой уходит туда,
// где экранирование слабее. В осесимметрии «другое место» = ДРУГОЙ РАДИУС, то есть
// кольцо. Это ровно тот артефакт, о котором предупреждает ERRATA B6 — и он показан
// явно, а не спрятан.
const D_SURF = 2e-8;                      // м^2/с — слабое растекание
const sigDens = new Float64Array(NR);     // Кл/м^2 на нижней поверхности
const sigTmp = new Float64Array(NR);
const cellA = new Float64Array(NR);
for (let i = 0; i < NR; i++) cellA[i] = Math.PI * (rf[i + 1] ** 2 - rf[i] ** 2);
let spotR = 0;                            // радиус текущего канала

// выбор места пробоя: минимум |sigma| (слабейшее экранирование), с лёгким
// предпочтением оси — там поле выше из-за затравочного пятна
function pickSpot() {
  let best = 0, bestVal = Infinity;
  for (let i = 0; i < NR; i++) {
    const v = Math.abs(sigDens[i]) + 3e-5 * (rc[i] / R);
    if (v < bestVal) { bestVal = v; best = i; }
  }
  return rc[best];
}

function depositSigma(dq, dt) {
  if (dq !== 0) {
    const tau = burning ? (st.t - tIgn) : 0;
    const sr = Math.min(SIG_R1, SIG_R0 + (SIG_R1 - SIG_R0) * Math.min(1, tau / 8e-8));
    let norm = 0;
    for (let i = 0; i < NR; i++) norm += Math.exp(-(((rc[i] - spotR) / sr) ** 2)) * cellA[i];
    for (let i = 0; i < NR; i++) sigDens[i] += dq * Math.exp(-(((rc[i] - spotR) / sr) ** 2)) / norm;
  }
  const a = D_SURF * dt;
  if (a <= 0) return;
  sigTmp.set(sigDens);
  for (let i = 0; i < NR; i++) {
    let flux = 0;
    if (i + 1 < NR) flux += 2 * Math.PI * rf[i + 1] * (sigTmp[i + 1] - sigTmp[i]) / (rc[i + 1] - rc[i]);
    if (i > 0) flux -= 2 * Math.PI * rf[i] * (sigTmp[i] - sigTmp[i - 1]) / (rc[i] - rc[i - 1]);
    sigDens[i] = sigTmp[i] + a * flux / cellA[i];
  }
}

function updateSigma() {
  let smax = 0, check = 0;
  for (let i = 0; i < NR; i++) {
    st.sigmaL[i] = sigDens[i];
    st.sigmaR[i] = -sigDens[i];
    check += sigDens[i] * cellA[i];
    const m = Math.abs(sigDens[i]);
    if (m > smax) smax = m;
  }
  st.sigmaMax = smax;
  st.sigmaCheckC = check;                  // должно равняться qSurf
}

// ---- рекордеры ---------------------------------------------------------------
const params = {
  gapMM: LG * 1e3, dielMM: D1 * 1e3, epsR: EPSR, radiusMM: R * 1e3,
  nr: NR, nz: NZ, nz1: NZ1, ngap: NGAP, nz2: NZ2, profile: PROFILE,
  U0kV: U0 / 1e3, freqKHz: FREQ / 1e3, pressureTorr: 760, tempK: 300,
  gammaIon: 0.02, seedBackground: 1e13, seedSpotAmp: 1e17, seedSpotSigmaUM: 40,
  photoIonization: true, photoEmission: true, photoDetachment: true,
  d6MassConvention: 'neutral',
  periods: PERIODS,
};
const meta = {
  synthetic: true,
  generator: 'sim2d/tests/synth-run.mjs',
  warning: 'СИНТЕТИЧЕСКИЙ ПРОГОН: правдоподобная форма сигнала, НЕ результат решения уравнений. ' +
           'Пригоден только для разработки плеера.',
  limitation: 'Осесимметричная модель описывает только центральный канал; ' +
              'азимутальные филаментационные моды отсутствуют принципиально.',
  capacitances: { C_gas_F: C_G, C_diel_F: C_D, C_cell_F: C_CELL },
  breakdown: { U_b_V: U_B, U_hold_V: U_HOLD },
};
const rec = createRecorders({
  dataDir: OUT, runId: RUN, params, meta,
  tEnd: PERIODS * TPER,
  // на пересечение зазора (10 нс) нужно разрешение мельче 1 нс, иначе головку не увидеть
  dtFrameMin: 2e-10, dtFrameMax: 1e-6,
  // синтетика — не научный продукт, полный уровень зажат до 200 МБ
  perLevel: { full: { maxBytes: Number(argv.maxBytes ?? 200e6) } },
});

// ---- интегрирование ----------------------------------------------------------
const T_END = PERIODS * TPER;
const DT_DARK = 2e-8, DT_BURN = 2e-11;
let steps = 0, framesFull = 0;
let qLastPeriod = null, uLastPeriod = null;   // замер замыкания Лиссажу по ПОСЛЕДНЕМУ периоду
const t0wall = Date.now();

st.t = 0;
UappPrev = 0;
updateSigma();
fillFields();
rec.record(st, true);
framesFull++;

while (st.t < T_END) {
  const dt = burning ? DT_BURN : DT_DARK;
  st.dt = dt;
  const t = st.t + dt;

  const Uapp = U0 * Math.sin(2 * Math.PI * FREQ * t);
  // напряжение на газе: делитель ёмкостей МИНУС экранирующее поле памяти (ERRATA A1)
  const Ugap = (C_D / (C_D + C_G)) * Uapp - qSurf / (C_D + C_G);

  // зажигание/гашение. Микроразряд КОНЕЧЕН по времени (~100 нс): канал разрушается
  // накоплением σ и уходом электронов. Пока |Uapp| растёт, за полупериод успевает
  // произойти НЕСКОЛЬКО микроразрядов — как в реальном филаментарном ДБР.
  if (!burning && Math.abs(Ugap) >= U_B && (t - tEnd_ > T_DEAD)) {
    burning = true; tIgn = t; headDir = Math.sign(Ugap) || 1;
    headZ = headDir > 0 ? ZS_LO : ZS_HI;
    spotR = pickSpot();
    nPulses++;
  }
  if (burning && (t - tIgn > T_PULSE || Math.abs(Ugap) <= U_HOLD * 0.999)) {
    burning = false; tEnd_ = t;
  }

  // ток проводимости: канал как нелинейный ключ, гасится ростом qSurf
  const Icond = burning ? Math.sign(Ugap) * GAMMA_SW * Math.max(0, Math.abs(Ugap) - U_HOLD) : 0;
  qSurf += Icond * dt;
  depositSigma(Icond * dt, dt);

  if (burning) {
    headZ += headDir * HEAD_V * dt;
    if (headZ > ZS_HI) headZ = ZS_HI;      // головка достигла анода и встала
    if (headZ < ZS_LO) headZ = ZS_LO;
  }

  const Idisp = C_CELL * (Uapp - UappPrev) / dt;
  const Itot = Icond + Idisp;
  st.Q += Itot * dt;

  // наработка озона от энерговклада: 100 г/(кВт*ч) = 3.5e17 молекул/Дж
  // (ERRATA V10: оценка ПОРЯДКА — эмпирический выход зависит от потока и температуры)
  o3Total += Math.abs(Icond * Ugap) * dt * O3_YIELD / V_GAP;
  // O3- накапливается после импульса, медленно убывает
  nO3mBulk += Math.abs(Icond) * dt * 4e28 - nO3mBulk * dt / 3e-4;
  neBulk = burning ? 1e15 : Math.max(1e13, neBulk * Math.exp(-dt / 5e-6));

  st.t = t; st.Uapp = Uapp; st.Ugap = Ugap;
  st.Icond = Icond; st.Idisp = Idisp; st.Itot = Itot;
  st.o3ppm = o3Total / N_GAS * 1e6;
  if (qLastPeriod === null && t >= T_END - TPER) { qLastPeriod = st.Q; uLastPeriod = Uapp; }
  UappPrev = Uapp;
  steps++;

  // тяжёлые 2D-поля заполняем ТОЛЬКО если кадр действительно будет записан
  if (rec.shouldRecord(st)) { updateSigma(); fillFields(); }
  else { st.neMax = neBulk; }
  const wrote = rec.record(st);
  if (wrote[0]) framesFull++;
}
updateSigma(); fillFields();
rec.record(st, true); framesFull++;

const paths = rec.close({
  checks: {
    // «Лиссажу замкнут» — прямая проверка синтетики; для настоящего прогона тут будут V1..V13
    lissajous_closure_C: st.Q - (qLastPeriod ?? 0),
    microdischarges: nPulses,
    C_cell_F: C_CELL, C_diel_F: C_D, U_b_V: U_B,
  },
});

// ---- отчёт -------------------------------------------------------------------
const wall = (Date.now() - t0wall) / 1e3;
console.log(`synth-run: профиль=${PROFILE} сетка=${NR}x${NZ} (газ ${NGAP}) периодов=${PERIODS}`);
console.log(`шагов: ${steps}, время расчёта: ${wall.toFixed(1)} с`);
for (const r of rec.recorders) {
  const dir = r.dir;
  const sz = (n) => fs.statSync(path.join(dir, n)).size;
  const total = sz('frames.bin') + sz('series.bin') + sz('manifest.json');
  console.log(`  [${r.level.padEnd(7)}] ${dir}`);
  console.log(`      кадров ${String(r.frameCount).padStart(5)}  stride ${r.frameStride} Б  ` +
              `frames.bin ${(sz('frames.bin') / 1e6).toFixed(2)} МБ  ` +
              `series.bin ${(sz('series.bin') / 1e6).toFixed(2)} МБ (${r.seriesCount} записей)  ` +
              `manifest ${(sz('manifest.json') / 1e6).toFixed(2)} МБ  ИТОГО ${(total / 1e6).toFixed(2)} МБ`);
}
const dQ = st.Q - (qLastPeriod ?? 0);
console.log(`микроразрядов: ${nPulses} за ${PERIODS} период(а)`);
console.log(`Лиссажу: dQ за последний период = ${dQ.toExponential(3)} Кл ` +
            `(${(Math.abs(dQ) / Math.max(1e-30, Math.abs(qSurf)) * 100).toFixed(2)} % от |qSurf|; замкнут, если ~0)`);
console.log(`сохранение заряда на поверхности: int(sigma dA) = ${st.sigmaCheckC.toExponential(4)} Кл, ` +
            `qSurf = ${qSurf.toExponential(4)} Кл, невязка ${Math.abs((st.sigmaCheckC - qSurf) / Math.max(1e-30, Math.abs(qSurf))).toExponential(1)}`);
console.log(`o3ppm = ${st.o3ppm.toFixed(1)}, |qSurf| max = ${Math.abs(qSurf).toExponential(3)} Кл`);
console.log('манифесты:', paths.join(', '));
