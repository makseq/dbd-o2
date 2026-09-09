// poisson2d.mjs — разделимый прямой решатель Пуассона/Гельмгольца
// для 2D осесимметричной (r,z) задачи ДБР.
//
// Опорные документы (в порядке приоритета):
//   docs/ERRATA.md      — B1..B5 (симметризация, ось, знак моды, полный стек, трансмиссивность),
//                         A1 (знак sigma), A3 (весовой потенциал), F3 (Гельмгольц с eps=1)
//   docs/NUMERICS_2D.md — §1 (сетка/метрики), §2 (разделимый решатель), §3 (полунеявность, PCG)
//   docs/PHOTO_PROCESSES.md §7 — Гельмгольц, модальное открытое ГУ
//
// СОГЛАШЕНИЯ, зафиксированные раз и навсегда:
//   * раскладка памяти  field[i*nz + j]  (z — быстрый индекс), NUMERICS_2D §1.2;
//   * уравнение         div(eps grad phi) = -rho;
//   * lam[k] — собственные значения оператора (-L_r), то есть lam[k] >= 0 (ERRATA B3),
//     и в модальном уравнении стоит МИНУС:
//         (A_z phi_k)[j] - eps[j]*dz[j]*lam[k]*phi_k[j] = -b_k[j];
//   * поверхностный заряд входит в правую часть С ТЕМ ЖЕ ЗНАКОМ, что и объёмный заряд
//     (ERRATA A1): sigma — это заряд, а не поток. Проверяется тестом P5 (экранирование).

export const EPS0 = 8.8541878128e-12;

// ───────────────────────────────────────────────────────────── сетка

/** Грани по r: сгущение к оси, постоянный коэффициент растяжения exp(betaR/nr). */
export function makeRGrid(nr, R, betaR = 1.6) {
  const rf = new Float64Array(nr + 1);
  if (Math.abs(betaR) < 1e-12) {
    for (let k = 0; k <= nr; k++) rf[k] = (R * k) / nr;
  } else {
    const den = Math.expm1(betaR);
    for (let k = 0; k <= nr; k++) rf[k] = (R * Math.expm1((betaR * k) / nr)) / den;
  }
  rf[0] = 0;
  rf[nr] = R;
  return rf;
}

/** Грани по z: барьер(геом. сгущение к поверхности) | газ(симметричный tanh) | барьер(зеркало). */
export function makeZGrid(nz1, ngap, nz2, d1, Lg, d2, betaG = 2.63, qB = 1.3147) {
  const zf = [0];
  const push = (x) => zf.push(x);
  const last = () => zf[zf.length - 1];

  // --- барьер 1: ячейки от металла к поверхности, ширина убывает как qB^-1
  const geom = Math.abs(qB - 1) > 1e-12;
  const dz0a = geom ? (d1 * (qB - 1)) / (Math.pow(qB, nz1) - 1) : d1 / nz1;
  for (let k = nz1 - 1; k >= 0; k--) push(last() + dz0a * (geom ? Math.pow(qB, k) : 1));
  zf[nz1] = d1; // точное положение границы газ/диэлектрик

  // --- газ: симметричный tanh (густо у обеих поверхностей)
  const z0 = d1;
  if (Math.abs(betaG) < 1e-12) {
    for (let k = 1; k <= ngap; k++) push(z0 + (Lg * k) / ngap);
  } else {
    const th = Math.tanh(betaG / 2);
    for (let k = 1; k <= ngap; k++) {
      const s = (2 * k) / ngap - 1;
      push(z0 + 0.5 * Lg * (1 + Math.tanh((betaG * s) / 2) / th));
    }
  }
  zf[nz1 + ngap] = d1 + Lg; // точное положение второй границы

  // --- барьер 2: зеркало блока 1
  const dz0b = geom ? (d2 * (qB - 1)) / (Math.pow(qB, nz2) - 1) : d2 / nz2;
  for (let k = 0; k < nz2; k++) push(last() + dz0b * (geom ? Math.pow(qB, k) : 1));
  zf[nz1 + ngap + nz2] = d1 + Lg + d2;

  return Float64Array.from(zf);
}

// ───────────────────────────────────────────────── симметричная трёхдиагональная задача

/**
 * Implicit QL с Wilkinson-сдвигом для симметричной трёхдиагональной матрицы.
 * d[0..n-1] — диагональ (на выходе собственные значения),
 * e[0..n-1] — поддиагональ, e[i] связывает i-1 и i (e[0] игнорируется),
 * z[comp*n + mode] — на входе единичная матрица, на выходе собственные векторы (по столбцам).
 */
export function tqli(d, e, n, z) {
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  for (let l = 0; l < n; l++) {
    let iter = 0;
    let m;
    do {
      for (m = l; m < n - 1; m++) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m]) <= Number.EPSILON * dd) break;
      }
      if (m !== l) {
        if (iter++ === 60) throw new Error('tqli: превышено число итераций');
        let g = (d[l + 1] - d[l]) / (2 * e[l]);
        let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? Math.abs(r) : -Math.abs(r)));
        let s = 1;
        let c = 1;
        let p = 0;
        let bailed = false;
        for (let i = m - 1; i >= l; i--) {
          let f = s * e[i];
          const b = c * e[i];
          r = Math.hypot(f, g);
          e[i + 1] = r;
          if (r === 0) {
            d[i + 1] -= p;
            e[m] = 0;
            bailed = true;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          for (let kk = 0; kk < n; kk++) {
            const base = kk * n;
            f = z[base + i + 1];
            z[base + i + 1] = s * z[base + i] + c * f;
            z[base + i] = c * z[base + i] - s * f;
          }
        }
        if (bailed) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0;
      }
    } while (m !== l);
  }
}

// ───────────────────────────────────────────────────────────── решатель

export class SeparableSolver {
  /**
   * @param {object} o
   *   nr, nz1, ngap, nz2 — число ячеек: r, барьер1, газ, барьер2
   *   radiusMM, dielMM, gapMM, epsR
   *   betaR, betaG, qB — параметры сгущения (0 / 1 => равномерно)
   */
  constructor(o = {}) {
    const nr = (this.nr = o.nr ?? 96);
    const nz1 = (this.nz1 = o.nz1 ?? 16);
    const ngap = (this.ngap = o.ngap ?? 192);
    const nz2 = (this.nz2 = o.nz2 ?? 16);
    const nz = (this.nz = nz1 + ngap + nz2);

    const R = (this.R = (o.radiusMM ?? 0.5) * 1e-3);
    const d1 = (this.d1 = (o.dielMM ?? 0.5) * 1e-3);
    const Lg = (this.Lg = (o.gapMM ?? 1.0) * 1e-3);
    const d2 = (this.d2 = o.diel2MM !== undefined ? o.diel2MM * 1e-3 : d1);
    this.Lz = d1 + Lg + d2;
    const epsRd = (this.epsRd = o.epsR ?? 9);

    // --- грани, центры, шаги
    const rf = (this.rf = makeRGrid(nr, R, o.betaR ?? 1.6));
    const zf = (this.zf = makeZGrid(nz1, ngap, nz2, d1, Lg, d2, o.betaG ?? 2.63, o.qB ?? 1.3147));

    const rc = (this.rc = new Float64Array(nr));
    const dr = (this.dr = new Float64Array(nr));
    for (let i = 0; i < nr; i++) {
      rc[i] = 0.5 * (rf[i] + rf[i + 1]);
      dr[i] = rf[i + 1] - rf[i];
    }
    const hr = (this.hr = new Float64Array(nr + 1)); // hr[i] = rc[i]-rc[i-1], hr[0], hr[nr] не используются
    for (let i = 1; i < nr; i++) hr[i] = rc[i] - rc[i - 1];

    const zc = (this.zc = new Float64Array(nz));
    const dz = (this.dz = new Float64Array(nz));
    for (let j = 0; j < nz; j++) {
      zc[j] = 0.5 * (zf[j] + zf[j + 1]);
      dz[j] = zf[j + 1] - zf[j];
    }
    const hz = (this.hz = new Float64Array(nz + 1));
    for (let j = 1; j < nz; j++) hz[j] = zc[j] - zc[j - 1];
    hz[0] = 0.5 * dz[0]; // до металла z=0
    hz[nz] = 0.5 * dz[nz - 1]; // до металла z=Lz

    // --- метрики
    const Acell = (this.Acell = new Float64Array(nr)); // площадь кольца = площадь z-грани
    for (let i = 0; i < nr; i++) Acell[i] = Math.PI * (rf[i + 1] * rf[i + 1] - rf[i] * rf[i]);
    const V = (this.V = new Float64Array(nr * nz));
    for (let i = 0; i < nr; i++) for (let j = 0; j < nz; j++) V[i * nz + j] = Acell[i] * dz[j];

    // --- материальные массивы (функция ТОЛЬКО от j — условие разделимости)
    this.JS_LO = nz1 - 1;
    this.JG0 = nz1;
    this.JG1 = nz1 + ngap - 1;
    this.JS_HI = nz1 + ngap;
    const epsr = (this.epsr = new Float64Array(nz));
    const eps = (this.eps = new Float64Array(nz));
    const isGas = (this.isGas = new Uint8Array(nz));
    for (let j = 0; j < nz; j++) {
      const gas = j >= this.JG0 && j <= this.JG1;
      isGas[j] = gas ? 1 : 0;
      epsr[j] = gas ? 1 : epsRd;
      eps[j] = EPS0 * epsr[j];
    }

    // --- радиальный оператор: K_r (симметричная трёхдиагональная жёсткость), W = diag(Acell)
    // K_r[i][i+1] = 2*pi*rf[i+1]/hr[i+1]; K_r[i][i] = -(сумма внедиагональных)
    // На оси rf[0]=0 => связь тождественно нулевая, деления на r нет (ERRATA B2).
    const trF = (this.trF = new Float64Array(nr + 1)); // внедиагональ K_r по граням
    for (let i = 1; i < nr; i++) trF[i] = (2 * Math.PI * rf[i]) / hr[i];
    trF[0] = 0;
    trF[nr] = 0;

    // симметризация S = W^{-1/2} K_r W^{-1/2}  (ERRATA B1)
    const sd = new Float64Array(nr);
    const se = new Float64Array(nr);
    for (let i = 0; i < nr; i++) sd[i] = -(trF[i] + trF[i + 1]) / Acell[i];
    for (let i = 1; i < nr; i++) se[i] = trF[i] / Math.sqrt(Acell[i - 1] * Acell[i]);

    const Q = (this.Q = new Float64Array(nr * nr));
    for (let i = 0; i < nr; i++) Q[i * nr + i] = 1;
    const mu = new Float64Array(sd);
    tqli(mu, new Float64Array(se), nr, Q);

    // lam[k] = -mu[k] >= 0 (ERRATA B3), сортировка по возрастанию lam
    const idx = Array.from({ length: nr }, (_, k) => k).sort((a, b) => mu[b] - mu[a]);
    const lam = (this.lam = new Float64Array(nr));
    const Qs = new Float64Array(nr * nr);
    for (let k = 0; k < nr; k++) {
      const src = idx[k];
      lam[k] = -mu[src];
      for (let i = 0; i < nr; i++) Qs[i * nr + k] = Q[i * nr + src];
    }
    // фиксация знака (первая компонента положительна) — воспроизводимость
    for (let k = 0; k < nr; k++) {
      if (Qs[0 * nr + k] < 0) for (let i = 0; i < nr; i++) Qs[i * nr + k] = -Qs[i * nr + k];
    }
    this.Q = Qs;
    const lamMax = Math.max(...lam);
    if (Math.abs(lam[0]) > 1e-10 * lamMax) throw new Error('нулевая мода потеряна');
    lam[0] = 0; // нулевая мода точна: K_r * 1 = 0
    this.lamMax = lamMax;

    // sqrt(Acell) — метрика прямого/обратного преобразования
    const sqA = (this.sqA = new Float64Array(nr));
    const isqA = (this.isqA = new Float64Array(nr));
    for (let i = 0; i < nr; i++) {
      sqA[i] = Math.sqrt(Acell[i]);
      isqA[i] = 1 / sqA[i];
    }
    // собственные векторы в исходной метрике: v = W^{-1/2} Q, v^T W v = I
    const vr = (this.vr = new Float64Array(nr * nr));
    for (let i = 0; i < nr; i++) for (let k = 0; k < nr; k++) vr[i * nr + k] = Qs[i * nr + k] * isqA[i];

    // --- электростатическая сборка (a = eps(z))
    this.aCellEs = new Float64Array(nr * nz);
    for (let i = 0; i < nr; i++) for (let j = 0; j < nz; j++) this.aCellEs[i * nz + j] = eps[j];
    this.TZ = new Float64Array(nr * (nz + 1));
    this.TR = new Float64Array((nr + 1) * nz);
    this.buildFaceT(this.aCellEs, this.TZ, this.TR);

    // трансмиссивности на единицу площади для z-направления (зависят только от j) — база (★★)
    this.tz = new Float64Array(nz + 1);
    for (let j = 1; j < nz; j++) this.tz[j] = 1 / ((0.5 * dz[j - 1]) / eps[j - 1] + (0.5 * dz[j]) / eps[j]);
    this.tz[0] = eps[0] / (0.5 * dz[0]);
    this.tz[nz] = eps[nz - 1] / (0.5 * dz[nz - 1]);
    this.azEs = new Float64Array(eps); // a(z) для разделимого решателя

    // веса дельта-источника sigma (§2.6): пересчитываются при смене dt/kappa
    this.wLo = { d: 0, g: 0 };
    this.wHi = { d: 0, g: 0 };
    this.computeSigmaWeights(this.azEs);

    // --- факторизация трёхдиагональных систем (★★) для всех мод
    this.fac = null;
    this.facAz = null;
    this.factorize(this.azEs);

    // --- рабочие буферы (ноль аллокаций в горячем цикле)
    this.phi = new Float64Array(nr * nz);
    this._b = new Float64Array(nr * nz);
    this._hat = new Float64Array(nr * nz);
    this._hat2 = new Float64Array(nr * nz);
    this._res = new Float64Array(nr * nz);
    this._p = new Float64Array(nr * nz);
    this._Ap = new Float64Array(nr * nz);
    this._zv = new Float64Array(nr * nz);
    this._w1 = new Float64Array(nr * nz);
    this._w2 = new Float64Array(nr * nz);
    this._w3 = new Float64Array(nr * nz);
    this._w4 = new Float64Array(nr * nz);
    this._w5 = new Float64Array(nr * nz);
    this._Ez = new Float64Array(nr * (nz + 1));
    this._Er = new Float64Array((nr + 1) * nz);
    this._helm = new Map();
    this._lineFac = null;

    // --- весовой потенциал Сато–Морроу (ERRATA A3): считается ОДИН РАЗ, БЕЗ плазмы
    this.dEff = Lg + (d1 + d2) / epsRd;
    this.phiL = new Float64Array(nr * nz);
    {
      const zero = new Float64Array(nr * nz);
      const zs = new Float64Array(nr);
      this.solvePoisson(zero, zs, zs, 1.0, this.phiL);
    }
    // аналитическое весовое поле (одномерно ТОЧНО, NUMERICS_2D §7.1)
    this.EzL_gas = 1 / this.dEff;
    this.EzL_diel = 1 / (epsRd * this.dEff);
  }

  // ────────────────────────────────────────── сборка граневых трансмиссивностей

  /**
   * TZ[i*(nz+1)+j] — трансмиссивность z-грани j (между j-1 и j), включая площадь Acell.
   * TR[i*nz+j]     — трансмиссивность r-грани i (между i-1 и i), включая площадь 2*pi*rf[i]*dz[j].
   * Формула — последовательное сложение сопротивлений полуячеек (ERRATA B5):
   *   T = A / ( 0.5*h1/a1 + 0.5*h2/a2 ).
   * При a = eps она тождественно совпадает с гармоническим средним eps на грани.
   */
  buildFaceT(aCell, TZ, TR) {
    const { nr, nz, dz, dr, rf, Acell } = this;
    for (let i = 0; i < nr; i++) {
      const ci = i * nz;
      const ti = i * (nz + 1);
      TZ[ti + 0] = (Acell[i] * aCell[ci + 0]) / (0.5 * dz[0]);
      for (let j = 1; j < nz; j++) {
        TZ[ti + j] = Acell[i] / ((0.5 * dz[j - 1]) / aCell[ci + j - 1] + (0.5 * dz[j]) / aCell[ci + j]);
      }
      TZ[ti + nz] = (Acell[i] * aCell[ci + nz - 1]) / (0.5 * dz[nz - 1]);
    }
    TR.fill(0);
    for (let i = 1; i < nr; i++) {
      const ti = i * nz;
      for (let j = 0; j < nz; j++) {
        const A = 2 * Math.PI * rf[i] * dz[j];
        TR[ti + j] = A / ((0.5 * dr[i - 1]) / aCell[(i - 1) * nz + j] + (0.5 * dr[i]) / aCell[i * nz + j]);
      }
    }
  }

  /** Веса дельта-источника sigma (§2.6) для текущего a(z). w_d + w_g == 1 тождественно. */
  computeSigmaWeights(az) {
    const { dz, JS_LO, JG0, JG1, JS_HI } = this;
    {
      const gd = az[JS_LO] / (0.5 * dz[JS_LO]);
      const gg = az[JG0] / (0.5 * dz[JG0]);
      this.wLo.d = gd / (gd + gg);
      this.wLo.g = gg / (gd + gg);
    }
    {
      const gg = az[JG1] / (0.5 * dz[JG1]);
      const gd = az[JS_HI] / (0.5 * dz[JS_HI]);
      this.wHi.g = gg / (gd + gg);
      this.wHi.d = gd / (gd + gg);
    }
  }

  /** LU-факторизация трёхдиагональных систем (★★) для каждой моды при заданном a(z). */
  factorize(az) {
    const { nr, nz, dz, lam } = this;
    const tz = new Float64Array(nz + 1);
    for (let j = 1; j < nz; j++) tz[j] = 1 / ((0.5 * dz[j - 1]) / az[j - 1] + (0.5 * dz[j]) / az[j]);
    tz[0] = az[0] / (0.5 * dz[0]);
    tz[nz] = az[nz - 1] / (0.5 * dz[nz - 1]);
    this.tzSep = tz;
    this.facAz = az;

    const cp = new Float64Array(nr * nz);
    const inv = new Float64Array(nr * nz);
    let maxCp = 0;
    for (let k = 0; k < nr; k++) {
      const off = k * nz;
      const lk = lam[k];
      let prevCp = 0;
      for (let j = 0; j < nz; j++) {
        const a = -tz[j]; // поддиагональ
        const c = -tz[j + 1]; // наддиагональ
        const b = tz[j] + tz[j + 1] + az[j] * dz[j] * lk; // диагональ (>0, M-матрица)
        const den = b - a * prevCp;
        const iv = 1 / den;
        inv[off + j] = iv;
        prevCp = c * iv;
        cp[off + j] = prevCp;
        const ac = Math.abs(prevCp);
        if (ac > maxCp) maxCp = ac;
      }
    }
    this.fac = { cp, inv, tz };
    this.maxCprime = maxCp; // ассерт |c'| < 1 => пивотинг не нужен (NUMERICS_2D §2.4)
    if (!(maxCp < 1)) throw new Error(`прогонка неустойчива: max|c'| = ${maxCp}`);
  }

  // ────────────────────────────────────────── прямое/обратное радиальное преобразование

  /** hat[k*nz+j] = sum_i Q[i][k] * sqrt(Acell[i]) * src[i*nz+j] */
  _forward(src, hat) {
    const { nr, nz, Q, sqA } = this;
    hat.fill(0);
    for (let i = 0; i < nr; i++) {
      const gi = i * nz;
      const qi = i * nr;
      const s = sqA[i];
      for (let k = 0; k < nr; k++) {
        const q = Q[qi + k] * s;
        if (q === 0) continue;
        const ko = k * nz;
        for (let j = 0; j < nz; j++) hat[ko + j] += q * src[gi + j];
      }
    }
  }

  /** dst[i*nz+j] = (sum_k Q[i][k] * hat[k*nz+j]) / sqrt(Acell[i]) */
  _backward(hat, dst) {
    const { nr, nz, Q, isqA } = this;
    dst.fill(0);
    for (let i = 0; i < nr; i++) {
      const gi = i * nz;
      const qi = i * nr;
      const s = isqA[i];
      for (let k = 0; k < nr; k++) {
        const q = Q[qi + k] * s;
        if (q === 0) continue;
        const ko = k * nz;
        for (let j = 0; j < nz; j++) dst[gi + j] += q * hat[ko + j];
      }
    }
  }

  // ────────────────────────────────────────── правая часть

  /**
   * b[i*nz+j] = rho*dz + вклад sigma + вклад Дирихле,
   * так что уравнение имеет вид  A_z phi - eps*dz*lam*phi = -b.
   */
  buildRhs(rho, sigLo, sigHi, U, b, tz = this.tzSep) {
    const { nr, nz, dz, JS_LO, JG0, JG1, JS_HI, wLo, wHi } = this;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      for (let j = 0; j < nz; j++) b[o + j] = rho[o + j] * dz[j];
    }
    // sigma — ЗАРЯД, знак тот же, что у объёмного rho (ERRATA A1)
    if (sigLo) {
      for (let i = 0; i < nr; i++) {
        const s = sigLo[i];
        b[i * nz + JS_LO] += wLo.d * s;
        b[i * nz + JG0] += wLo.g * s;
      }
    }
    if (sigHi) {
      for (let i = 0; i < nr; i++) {
        const s = sigHi[i];
        b[i * nz + JG1] += wHi.g * s;
        b[i * nz + JS_HI] += wHi.d * s;
      }
    }
    if (U !== 0) {
      const t0 = tz[0] * U;
      for (let i = 0; i < nr; i++) b[i * nz + 0] += t0;
    }
  }

  // ────────────────────────────────────────── основной решатель

  /**
   * Точный (прямой) электростатический решатель.
   * @param {Float64Array} rho  плотность заряда [Кл/м^3], nr*nz (в барьерах 0)
   * @param {Float64Array} sigLo поверхностный заряд на нижней границе газ/диэлектрик [Кл/м^2], nr
   * @param {Float64Array} sigHi поверхностный заряд на верхней границе, nr
   * @param {number} U напряжение на активном электроде (z=0); z=Lz заземлён
   * @returns {Float64Array} phi, nr*nz
   */
  solvePoisson(rho, sigLo, sigHi, U, out = this.phi) {
    if (this.facAz !== this.azEs) {
      this.computeSigmaWeights(this.azEs);
      this.factorize(this.azEs);
    }
    const b = this._b;
    this.buildRhs(rho, sigLo, sigHi, U, b, this.fac.tz);
    this._forward(b, this._hat);
    this._solveModes(this._hat, this._hat2);
    this._backward(this._hat2, out);
    return out;
  }

  /** Прогонка Томаса по всем модам: решает (тридиаг) x = rhs, где rhs = bhat. */
  _solveModes(bhat, xhat) {
    const { nr, nz } = this;
    const { cp, inv, tz } = this.fac;
    for (let k = 0; k < nr; k++) {
      const off = k * nz;
      // прямой ход
      let prev = bhat[off] * inv[off];
      xhat[off] = prev;
      for (let j = 1; j < nz; j++) {
        const a = -tz[j];
        prev = (bhat[off + j] - a * prev) * inv[off + j];
        xhat[off + j] = prev;
      }
      // обратный ход
      for (let j = nz - 2; j >= 0; j--) xhat[off + j] -= cp[off + j] * xhat[off + j + 1];
    }
  }

  // ────────────────────────────────────────── полная матрица (для PCG/тестов)

  /** out = A*phi (несделённая симметричная форма, однородные Дирихле-призраки). */
  applyMatrix(phi, out, TZ = this.TZ, TR = this.TR) {
    const { nr, nz } = this;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const tzo = i * (nz + 1);
      const trL = i * nz;
      const trR = (i + 1) * nz;
      for (let j = 0; j < nz; j++) {
        const p = phi[o + j];
        const up = j + 1 < nz ? phi[o + j + 1] : 0;
        const dn = j > 0 ? phi[o + j - 1] : 0;
        let s = TZ[tzo + j + 1] * (up - p) - TZ[tzo + j] * (p - dn);
        if (i + 1 < nr) s += TR[trR + j] * (phi[o + nz + j] - p);
        if (i > 0) s -= TR[trL + j] * (p - phi[o - nz + j]);
        out[o + j] = s;
      }
    }
    return out;
  }

  /**
   * Правая часть в несделённой форме: RHS = -(rho*V + w*sigma*Acell + T_dir*U).
   * Веса w_d/w_g считаются по ФАКТИЧЕСКИМ коэффициентам ячеек aCell (ERRATA B5 / §3.6):
   * в полунеявной схеме они радиально неоднородны, но w_d+w_g == 1 в каждой колонке.
   */
  buildRhsFull(rho, sigLo, sigHi, U, out, TZ = this.TZ, aCell = this.aCellEs) {
    const { nr, nz, dz, Acell, JS_LO, JG0, JG1, JS_HI } = this;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const A = Acell[i];
      for (let j = 0; j < nz; j++) out[o + j] = -rho[o + j] * A * dz[j];
      if (sigLo) {
        const gd = aCell[o + JS_LO] / (0.5 * dz[JS_LO]);
        const gg = aCell[o + JG0] / (0.5 * dz[JG0]);
        const s = sigLo[i] * A;
        out[o + JS_LO] -= (gd / (gd + gg)) * s;
        out[o + JG0] -= (gg / (gd + gg)) * s;
      }
      if (sigHi) {
        const gg = aCell[o + JG1] / (0.5 * dz[JG1]);
        const gd = aCell[o + JS_HI] / (0.5 * dz[JS_HI]);
        const s = sigHi[i] * A;
        out[o + JG1] -= (gg / (gd + gg)) * s;
        out[o + JS_HI] -= (gd / (gd + gg)) * s;
      }
      out[o] -= TZ[i * (nz + 1)] * U;
    }
    return out;
  }

  // ────────────────────────────────────────── поля

  /** Ez на z-гранях: EzOut[i*(nz+1)+j], Ez = -dphi/dz. */
  computeEz(phi, U, out = this._Ez) {
    const { nr, nz, hz } = this;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const t = i * (nz + 1);
      out[t] = -(phi[o] - U) / hz[0];
      for (let j = 1; j < nz; j++) out[t + j] = -(phi[o + j] - phi[o + j - 1]) / hz[j];
      out[t + nz] = -(0 - phi[o + nz - 1]) / hz[nz];
    }
    return out;
  }

  /** Er на r-гранях: ErOut[i*nz+j], i=0..nr; Er[0]=Er[nr]=0 (Нейман). */
  computeEr(phi, out = this._Er) {
    const { nr, nz, hr } = this;
    out.fill(0);
    for (let i = 1; i < nr; i++) {
      const o = i * nz;
      for (let j = 0; j < nz; j++) out[o + j] = -(phi[o + j] - phi[o - nz + j]) / hr[i];
    }
    return out;
  }

  // ────────────────────────────────────────── ГЕЛЬМГОЛЬЦ (ERRATA F3)

  /**
   * Решает  lap(S) - kappa^2 * S = -q  с КОЭФФИЦИЕНТОМ 1 (НЕ eps(z)).
   * Сборка ОТДЕЛЬНАЯ от Пуассона (ERRATA F3): переиспользуется только радиальная
   * диагонализация и прогонка.
   *
   * @param {Float64Array} q   источник [ед./м^3], nr*nz (вне области решения игнорируется)
   * @param {number} kappa     обратная длина поглощения [1/м]
   * @param {object} opts      { domain: 'gas'|'full', bc: 'open'|'dirichlet', out }
   *   domain='gas' (дефолт): только газовый зазор, диэлектрик для ВУФ непрозрачен;
   *   bc='open' (дефолт): модальное открытое условие dS/dn + sqrt(lam_k + kappa^2) S = 0
   *                       (PHOTO_PROCESSES §7.3).
   * @returns {Float64Array} S, nr*nz (нули вне области решения)
   */
  solveHelmholtz(q, kappa, opts = {}) {
    const domain = opts.domain ?? 'gas';
    const bc = opts.bc ?? 'open';
    const out = opts.out ?? new Float64Array(this.nr * this.nz);
    const key = `${domain}|${bc}|${kappa}`;
    let h = this._helm.get(key);
    if (!h) {
      h = this._buildHelmholtz(kappa, domain, bc);
      this._helm.set(key, h);
    }
    const { nr, nz, dz, sqA, isqA, Q } = this;
    const { j0, nzh, cp, inv, sub } = h;

    // прямое преобразование источника: qhat[k][l] = sum_i Q[i][k]*sqrt(A_i)*q[i][j0+l]*dz[j0+l]
    const hat = h.hat;
    hat.fill(0);
    for (let i = 0; i < nr; i++) {
      const gi = i * nz + j0;
      const qi = i * nr;
      const s = sqA[i];
      for (let k = 0; k < nr; k++) {
        const w = Q[qi + k] * s;
        if (w === 0) continue;
        const ko = k * nzh;
        for (let l = 0; l < nzh; l++) hat[ko + l] += w * q[gi + l] * dz[j0 + l];
      }
    }
    // прогонки
    const xh = h.xhat;
    for (let k = 0; k < nr; k++) {
      const off = k * nzh;
      let prev = hat[off] * inv[off];
      xh[off] = prev;
      for (let l = 1; l < nzh; l++) {
        prev = (hat[off + l] - sub[l] * prev) * inv[off + l];
        xh[off + l] = prev;
      }
      for (let l = nzh - 2; l >= 0; l--) xh[off + l] -= cp[off + l] * xh[off + l + 1];
    }
    // обратное преобразование
    out.fill(0);
    for (let i = 0; i < nr; i++) {
      const gi = i * nz + j0;
      const qi = i * nr;
      const s = isqA[i];
      for (let k = 0; k < nr; k++) {
        const w = Q[qi + k] * s;
        if (w === 0) continue;
        const ko = k * nzh;
        for (let l = 0; l < nzh; l++) out[gi + l] += w * xh[ko + l];
      }
    }
    return out;
  }

  _buildHelmholtz(kappa, domain, bc) {
    const { nr, nz, dz, hz, lam } = this;
    const j0 = domain === 'gas' ? this.JG0 : 0;
    const j1 = domain === 'gas' ? this.JG1 : nz - 1;
    const nzh = j1 - j0 + 1;
    const k2 = kappa * kappa;

    // трансмиссивности на единицу площади, коэффициент 1
    const t = new Float64Array(nzh + 1);
    for (let l = 1; l < nzh; l++) t[l] = 1 / hz[j0 + l];
    if (bc === 'dirichlet') {
      t[0] = 1 / (0.5 * dz[j0]);
      t[nzh] = 1 / (0.5 * dz[j1]);
    } else {
      t[0] = 0;
      t[nzh] = 0;
    }
    const sub = new Float64Array(nzh);
    for (let l = 0; l < nzh; l++) sub[l] = -t[l];

    const cp = new Float64Array(nr * nzh);
    const inv = new Float64Array(nr * nzh);
    for (let k = 0; k < nr; k++) {
      const c0 = lam[k] + k2;
      if (bc === 'open' && !(c0 > 0)) {
        throw new Error('solveHelmholtz: kappa=0 с открытым ГУ даёт вырожденную нулевую моду');
      }
      const cop = bc === 'open' ? Math.sqrt(c0) : 0;
      const off = k * nzh;
      let prevCp = 0;
      for (let l = 0; l < nzh; l++) {
        const j = j0 + l;
        let b = t[l] + t[l + 1] + dz[j] * c0;
        if (bc === 'open' && (l === 0 || l === nzh - 1)) b += cop;
        const a = -t[l];
        const c = -t[l + 1];
        const den = b - a * prevCp;
        const iv = 1 / den;
        inv[off + l] = iv;
        prevCp = c * iv;
        cp[off + l] = prevCp;
      }
    }
    return {
      j0,
      j1,
      nzh,
      cp,
      inv,
      sub,
      hat: new Float64Array(nr * nzh),
      xhat: new Float64Array(nr * nzh),
    };
  }

  // ────────────────────────────────────────── ПОЛУНЕЯВНАЯ ПОПРАВКА + PCG

  /**
   * Полунеявный Пуассон: div[(eps + dt*kappa) grad phi] = -rho + dt*div(J* + kappa*grad phi^n).
   * SPD-форма, разделимый решатель S1 + z-линейный блок-Якоби как предобуславливатель PCG.
   *
   * @param {object} a
   *   rho, sigLo, sigHi, U      — как в solvePoisson
   *   dt                        — шаг [с]
   *   kappaCell                 — проводимость плазмы по ячейкам [См/м], nr*nz
   *   Fz, Fr                    — (опц.) поверхностные плотности (J* - kappa*E^n) на z/r-гранях
   *                               [А/м^2], размеры nr*(nz+1) и (nr+1)*nz
   *   phiGuess                  — (опц.) начальное приближение (тёплый старт)
   *   tol, maxIter
   * @returns {{phi:Float64Array, iters:number, resid:number}}
   */
  solveSemiImplicit(a) {
    const { nr, nz, Acell, dz, rf, eps } = this;
    const n = nr * nz;
    const dt = a.dt ?? 0;
    const kap = a.kappaCell;

    // a(r,z) = eps(z) + dt*kappa(r,z)
    const aCell = this._aCell || (this._aCell = new Float64Array(n));
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      for (let j = 0; j < nz; j++) aCell[o + j] = eps[j] + (kap ? dt * kap[o + j] : 0);
    }
    const TZ = this._TZs || (this._TZs = new Float64Array(nr * (nz + 1)));
    const TR = this._TRs || (this._TRs = new Float64Array((nr + 1) * nz));
    this.buildFaceT(aCell, TZ, TR);

    // разделимый предобуславливатель: ar(r) ≡ 1, az(z) = eps + dt*<kappa>_A(z)  (§3.4)
    const az = this._az || (this._az = new Float64Array(nz));
    const Atot = Math.PI * this.R * this.R;
    for (let j = 0; j < nz; j++) {
      let s = 0;
      if (kap) for (let i = 0; i < nr; i++) s += kap[i * nz + j] * Acell[i];
      az[j] = eps[j] + (dt * s) / Atot;
    }
    this.computeSigmaWeights(az);
    this.factorize(az);

    // правая часть
    const rhs = this._res;
    this.buildRhsFull(a.rho, a.sigLo, a.sigHi, a.U, rhs, TZ, aCell);
    if (a.Fz || a.Fr) {
      const Fz = a.Fz;
      const Fr = a.Fr;
      for (let i = 0; i < nr; i++) {
        const o = i * nz;
        const t = i * (nz + 1);
        for (let j = 0; j < nz; j++) {
          let d = 0;
          if (Fz) d += (Fz[t + j + 1] - Fz[t + j]) * Acell[i];
          if (Fr) d += Fr[(i + 1) * nz + j] * 2 * Math.PI * rf[i + 1] * dz[j] - Fr[i * nz + j] * 2 * Math.PI * rf[i] * dz[j];
          rhs[o + j] += dt * d;
        }
      }
    }

    // z-линейный блок-Якоби с ИСТИННЫМИ коэффициентами
    this._buildLineJacobi(TZ, TR);

    // CG на (-A) x = (-rhs) : SPD
    const x = this.phi;
    if (a.phiGuess) x.set(a.phiGuess);
    else x.fill(0);
    const bvec = this._w1;
    for (let i = 0; i < n; i++) bvec[i] = -rhs[i];

    const Ap = this._Ap;
    const r = this._w2;
    const p = this._p;
    const zv = this._zv;
    this.applyMatrix(x, Ap, TZ, TR);
    for (let i = 0; i < n; i++) r[i] = bvec[i] + Ap[i]; // r = b - (-A)x
    let rn = Math.sqrt(dot(r, r));
    // критерий: абсолютный tol, либо относительный relTol от ||b||
    const tol = a.tol ?? (a.relTol ?? 1e-10) * Math.sqrt(dot(bvec, bvec));
    const maxIter = a.maxIter ?? 60;
    let iters = 0;
    if (rn > tol) {
      this._applyPrec(r, zv, TZ, TR);
      p.set(zv);
      let rz = dot(r, zv);
      for (; iters < maxIter; iters++) {
        this.applyMatrix(p, Ap, TZ, TR);
        for (let i = 0; i < n; i++) Ap[i] = -Ap[i];
        const alpha = rz / dot(p, Ap);
        for (let i = 0; i < n; i++) {
          x[i] += alpha * p[i];
          r[i] -= alpha * Ap[i];
        }
        rn = Math.sqrt(dot(r, r));
        if (rn <= tol) {
          iters++;
          break;
        }
        this._applyPrec(r, zv, TZ, TR);
        const rzNew = dot(r, zv);
        const beta = rzNew / rz;
        rz = rzNew;
        for (let i = 0; i < n; i++) p[i] = zv[i] + beta * p[i];
      }
    }
    return { phi: x, iters, resid: rn };
  }

  /** M^{-1}: симметричная мультипликативная комбинация B (z-линия) и S1 (разделимый). */
  _applyPrec(b, x, TZ, TR) {
    const n = this.nr * this.nz;
    const t = this._w3;
    const r2 = this._w4;
    this._applyLineJacobi(b, x);
    this.applyMatrix(x, t, TZ, TR);
    for (let i = 0; i < n; i++) r2[i] = b[i] + t[i];
    this._applyS1(r2, t);
    for (let i = 0; i < n; i++) x[i] += t[i];
    this.applyMatrix(x, t, TZ, TR);
    for (let i = 0; i < n; i++) r2[i] = b[i] + t[i];
    this._applyLineJacobi(r2, t);
    for (let i = 0; i < n; i++) x[i] += t[i];
  }

  /** S1: разделимый прямой решатель для (-A_sep) x = b (b — несделённая правая часть). */
  _applyS1(b, x) {
    const { nr, nz, Acell } = this;
    const tmp = this._w5;
    // перевод в «поделённую» форму: b_div = -b/Acell (уравнение A_z phi - ... = -b_div)
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const iA = 1 / Acell[i];
      for (let j = 0; j < nz; j++) tmp[o + j] = b[o + j] * iA;
    }
    this._forward(tmp, this._hat2);
    this._solveModes(this._hat2, this._hat2);
    this._backward(this._hat2, x);
  }

  _buildLineJacobi(TZ, TR) {
    const { nr, nz } = this;
    if (!this._lineFac) {
      this._lineFac = { cp: new Float64Array(nr * nz), inv: new Float64Array(nr * nz), sub: new Float64Array(nr * nz) };
    }
    const { cp, inv, sub } = this._lineFac;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      const t = i * (nz + 1);
      let prevCp = 0;
      for (let j = 0; j < nz; j++) {
        const diag = TZ[t + j] + TZ[t + j + 1] + TR[i * nz + j] + TR[(i + 1) * nz + j];
        const am = j > 0 ? -TZ[t + j] : 0;
        const ap = j < nz - 1 ? -TZ[t + j + 1] : 0;
        sub[o + j] = am;
        const den = diag - am * prevCp;
        const iv = 1 / den;
        inv[o + j] = iv;
        prevCp = ap * iv;
        cp[o + j] = prevCp;
      }
    }
  }

  _applyLineJacobi(b, x) {
    const { nr, nz } = this;
    const { cp, inv, sub } = this._lineFac;
    for (let i = 0; i < nr; i++) {
      const o = i * nz;
      let prev = b[o] * inv[o];
      x[o] = prev;
      for (let j = 1; j < nz; j++) {
        prev = (b[o + j] - sub[o + j] * prev) * inv[o + j];
        x[o + j] = prev;
      }
      for (let j = nz - 2; j >= 0; j--) x[o + j] -= cp[o + j] * x[o + j + 1];
    }
  }

  /** Сброс к чисто электростатической сборке (после полунеявных вызовов). */
  resetElectrostatic() {
    this.computeSigmaWeights(this.azEs);
    this.factorize(this.azEs);
  }
}

function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
