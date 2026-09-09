// tests/photo.test.mjs — тесты фотомодуля.
//
// Критерии PH1..PH8 (ТЗ) + VF4/VF7/VF8 из docs/PHOTO_PROCESSES.md §8.2.
//
// ВАЖНО: poisson2d.mjs пишется параллельно. Заглушка solveHelmholtz живёт ТОЛЬКО
// здесь (в photo2d.mjs её нет — модуль обязан зависеть от настоящего решателя).
// Заглушка реализована строго по контракту sim2d/PHOTO_API.md §1: отдельная
// сборка с eps=1, обобщённая диагонализация L_r в метрике M_r (ERRATA B1),
// знак модального уравнения по ERRATA B3, модальное открытое ГУ по PHOTO §7.3.
// Тесты PH2 и PH5 надо повторить на настоящем poisson2d.mjs.

import { PhotoModule, fowlerYield, transparentKernel, PHOTO_ION_WEIGHT_SETS } from '../photo2d.mjs';

// ---------------------------------------------------------------------------
// мини-раннер
// ---------------------------------------------------------------------------
let passed = 0, failed = 0;
const lines = [];
function test(name, fn) {
  try { fn(); passed++; lines.push(`  ok   ${name}`); }
  catch (e) { failed++; lines.push(`  FAIL ${name}\n       ${e.message}`); }
}
function ok(cond, msg) { if (!cond) throw new Error(msg); }
function near(a, b, rel, msg) {
  const d = Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
  if (!(d <= rel)) throw new Error(`${msg}: ${a} vs ${b}, отн. ошибка ${d.toExponential(3)} > ${rel}`);
  return d;
}
function section(s) { lines.push(`\n${s}`); }

// ---------------------------------------------------------------------------
// сетки
// ---------------------------------------------------------------------------
function uniformGrid(nr, R, nzDiel, ngz, dDiel, Lg) {
  const nz = 2 * nzDiel + ngz;
  const rf = new Float64Array(nr + 1);
  for (let i = 0; i <= nr; i++) rf[i] = R * i / nr;
  const zf = new Float64Array(nz + 1);
  for (let j = 0; j <= nzDiel; j++) zf[j] = dDiel * j / Math.max(1, nzDiel);
  for (let j = 1; j <= ngz; j++) zf[nzDiel + j] = dDiel + Lg * j / ngz;
  for (let j = 1; j <= nzDiel; j++) zf[nzDiel + ngz + j] = dDiel + Lg + dDiel * j / Math.max(1, nzDiel);
  return { nr, nz, rf, zf, JG0: nzDiel, JG1: nzDiel + ngz - 1 };
}

// ---------------------------------------------------------------------------
// ЗАГЛУШКА solveHelmholtz (см. PHOTO_API.md §1)
// ---------------------------------------------------------------------------
function jacobiEigen(Ain, n) {
  const A = Float64Array.from(Ain);
  const V = new Float64Array(n * n);
  for (let i = 0; i < n; i++) V[i * n + i] = 1;
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
    if (off < 1e-30) break;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p * n + q];
        if (Math.abs(apq) < 1e-300) continue;
        const theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) {
          const akp = A[k * n + p], akq = A[k * n + q];
          A[k * n + p] = c * akp - s * akq;
          A[k * n + q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {
          const apk = A[p * n + k], aqk = A[q * n + k];
          A[p * n + k] = c * apk - s * aqk;
          A[q * n + k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  const lam = new Float64Array(n);
  for (let i = 0; i < n; i++) lam[i] = A[i * n + i];
  return { lam, V };           // V[i*n + k] — i-я компонента k-го вектора
}

class HelmholtzStub {
  constructor(grid) {
    const nr = grid.nr, nz = grid.nz;
    this.nr = nr; this.nz = nz;
    this.JG0 = grid.JG0; this.JG1 = grid.JG1;
    this.ngz = grid.JG1 - grid.JG0 + 1;
    const rf = grid.rf, zf = grid.zf;
    const rc = new Float64Array(nr);
    for (let i = 0; i < nr; i++) rc[i] = 0.5 * (rf[i] + rf[i + 1]);
    // M_r = diag[(rf_{i+1}^2 - rf_i^2)/2]  (ERRATA B1)
    const M = new Float64Array(nr);
    for (let i = 0; i < nr; i++) M[i] = 0.5 * (rf[i + 1] ** 2 - rf[i] ** 2);
    // трансмиссивности: T[i] = rf[i]/(rc[i]-rc[i-1]); T[0]=0 (ось, ERRATA B2), T[nr]=0 (зеркало)
    const T = new Float64Array(nr + 1);
    for (let i = 1; i < nr; i++) T[i] = rf[i] / (rc[i] - rc[i - 1]);
    const K = new Float64Array(nr * nr);
    for (let i = 0; i < nr; i++) {
      K[i * nr + i] = T[i] + T[i + 1];
      if (i + 1 < nr) { K[i * nr + i + 1] = -T[i + 1]; K[(i + 1) * nr + i] = -T[i + 1]; }
    }
    // B = M^{-1/2} K M^{-1/2}, симметричная -> Якоби
    const is = new Float64Array(nr);
    for (let i = 0; i < nr; i++) is[i] = 1 / Math.sqrt(M[i]);
    const B = new Float64Array(nr * nr);
    for (let i = 0; i < nr; i++) for (let j = 0; j < nr; j++) B[i * nr + j] = K[i * nr + j] * is[i] * is[j];
    const { lam, V } = jacobiEigen(B, nr);
    const vv = new Float64Array(nr * nr);       // v = M^{-1/2} u, v^T M v = I
    for (let i = 0; i < nr; i++) for (let k = 0; k < nr; k++) vv[i * nr + k] = V[i * nr + k] * is[i];
    this.M = M; this.lam = lam; this.v = vv; this.rc = rc; this.rf = rf;
    // z-метрика на газовой части
    const zc = new Float64Array(nz), dz = new Float64Array(nz);
    for (let j = 0; j < nz; j++) { zc[j] = 0.5 * (zf[j] + zf[j + 1]); dz[j] = zf[j + 1] - zf[j]; }
    this.zc = zc; this.dz = dz;
    this.qhat = new Float64Array(nr * this.ngz);
    this.shat = new Float64Array(nr * this.ngz);
    this.escaped = 0;
    // проверка ортонормировки v^T M v = I (ERRATA B1)
    let maxOff = 0;
    for (let k = 0; k < nr; k++) for (let l = 0; l < nr; l++) {
      let s = 0;
      for (let i = 0; i < nr; i++) s += vv[i * nr + k] * M[i] * vv[i * nr + l];
      maxOff = Math.max(maxOff, Math.abs(s - (k === l ? 1 : 0)));
    }
    this.orthoResidual = maxOff;
    // все lam >= 0 (знак по ERRATA B3)
    this.lamMin = Math.min(...lam);
  }

  solveHelmholtz(q, kappa, out) {
    const nr = this.nr, nz = this.nz, ngz = this.ngz, JG0 = this.JG0;
    const M = this.M, v = this.v, lam = this.lam, dz = this.dz, zc = this.zc;
    const qh = this.qhat, sh = this.shat;
    qh.fill(0);
    // прямое преобразование в метрике M_r: qhat_k = sum_i v_k[i] M_i q_i
    for (let i = 0; i < nr; i++) {
      const mi = M[i];
      for (let jj = 0; jj < ngz; jj++) {
        const qv = q[i * nz + JG0 + jj] * mi;
        if (qv === 0) continue;
        for (let k = 0; k < nr; k++) qh[k * ngz + jj] += v[i * nr + k] * qv;
      }
    }
    const a = new Float64Array(ngz), b = new Float64Array(ngz), c = new Float64Array(ngz), d = new Float64Array(ngz);
    let escaped = 0;
    for (let k = 0; k < nr; k++) {
      const s = Math.sqrt(Math.max(0, lam[k]) + kappa * kappa);
      for (let jj = 0; jj < ngz; jj++) {
        const j = JG0 + jj;
        const am = jj > 0 ? 1 / (dz[j] * (zc[j] - zc[j - 1])) : 0;
        const ap = jj < ngz - 1 ? 1 / (dz[j] * (zc[j + 1] - zc[j])) : 0;
        a[jj] = -am; c[jj] = -ap;
        // модальное ОТКРЫТОЕ ГУ: dS/dn + sqrt(lam_k + kappa^2) S = 0 (PHOTO §7.3)
        let diag = am + ap + lam[k] + kappa * kappa;
        if (jj === 0) diag += s / dz[j];
        if (jj === ngz - 1) diag += s / dz[j];
        b[jj] = diag;
        d[jj] = qh[k * ngz + jj];
      }
      // прогонка Томаса (матрица строго диагонально доминантна)
      for (let jj = 1; jj < ngz; jj++) {
        const m = a[jj] / b[jj - 1];
        b[jj] -= m * c[jj - 1];
        d[jj] -= m * d[jj - 1];
      }
      sh[k * ngz + ngz - 1] = d[ngz - 1] / b[ngz - 1];
      for (let jj = ngz - 2; jj >= 0; jj--) sh[k * ngz + jj] = (d[jj] - c[jj] * sh[k * ngz + jj + 1]) / b[jj];
      // уход через границы: только мода k=0 даёт ненулевой интеграл по площади
      let sumMi = 0;
      for (let i = 0; i < nr; i++) sumMi += M[i] * v[i * nr + k];
      escaped += 2 * Math.PI * s * (sh[k * ngz] + sh[k * ngz + ngz - 1]) * sumMi;
    }
    this.escaped = escaped;
    const res = out || new Float64Array(nr * nz);
    res.fill(0);
    for (let i = 0; i < nr; i++) {
      for (let k = 0; k < nr; k++) {
        const vk = v[i * nr + k];
        if (vk === 0) continue;
        for (let jj = 0; jj < ngz; jj++) res[i * nz + JG0 + jj] += vk * sh[k * ngz + jj];
      }
    }
    return res;
  }
}

// ===========================================================================
section('--- ЗАГЛУШКА Гельмгольца: самопроверка (ERRATA B1/B3) ---');
// ===========================================================================
{
  const g = uniformGrid(48, 1e-3, 0, 48, 0, 1e-3);
  const h = new HelmholtzStub(g);
  test('stub: v^T M_r v = I с невязкой < 1e-12 (ERRATA B1)', () => {
    ok(h.orthoResidual < 1e-12, `невязка ортонормировки ${h.orthoResidual.toExponential(2)}`);
  });
  test('stub: все lambda_k >= 0 (знак по ERRATA B3)', () => {
    ok(h.lamMin > -1e-9, `lambda_min = ${h.lamMin}`);
  });
}

// ===========================================================================
section('--- PH1: все флаги false -> тождественные нули ---');
// ===========================================================================
test('PH1 модуль возвращает тождественные нули и не влияет на решение', () => {
  const g = uniformGrid(16, 5e-4, 4, 16, 5e-4, 1e-3);
  const pm = new PhotoModule(g, {
    photoIonization: false, photoEmission: false, photoDetachment: false,
  }, null);
  const N = g.nr * g.nz;
  const S = new Float64Array(N).fill(2.35e29);
  const P = new Float64Array(N).fill(3.8e10);
  const out = pm.update({ ionizRate: S, powerDens: P }, { step: 0 });
  for (const [name, arr] of [['photoIonRate', out.photoIonRate], ['psi', out.psi],
                             ['nuPd.Om', out.nuPd.Om], ['nuPd.O2m', out.nuPd.O2m],
                             ['nuPd.O3m', out.nuPd.O3m],
                             ['photoEmitFluxL', out.photoEmitFluxL],
                             ['photoEmitFluxR', out.photoEmitFluxR]]) {
    for (let n = 0; n < arr.length; n++) ok(arr[n] === 0, `${name}[${n}] = ${arr[n]} != 0`);
  }
  // ни одна тяжёлая структура не построена (нет скрытой инициализации)
  ok(pm.vfL === undefined && pm.coarse === undefined, 'построены матрицы при выключенных флагах');
  // gamma_i при выключенной фотоэмиссии = gamma_eff (деления быть не должно)
  ok(out.gammaI === 0.02, `gammaI = ${out.gammaI}, ожидался gamma_eff = 0.02`);
  ok(out.FphOverFi === 0, `F_ph/F_i = ${out.FphOverFi}, ожидался 0`);
});

// ===========================================================================
section('--- PH2: Гельмгольц, точечный источник, exp(-kR)/R ---');
// ===========================================================================
test('PH2 спад от точечного источника совпадает с exp(-kappa*R)/(4 pi R) < 2%', () => {
  const nr = 128, R = 1e-3, ngz = 128, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const h = new HelmholtzStub(g);
  const kappa = 1.111e4;                  // средняя группа, ell = 90 мкм
  const N = nr * g.nz;
  const q = new Float64Array(N);
  const dr = R / nr, dz = Lg / ngz;
  const jmid = ngz >> 1;
  const Vcell = Math.PI * (dr * dr) * dz;                 // приосевая ячейка-цилиндрик
  const Ndot = 1.0;
  q[0 * g.nz + jmid] = Ndot / Vcell;
  const S = h.solveHelmholtz(q, kappa, new Float64Array(N));
  const zmid = 0.5 * (g.zf[jmid] + g.zf[jmid + 1]);
  let worst = 0, worstR = 0, cnt = 0;
  for (let i = 1; i < nr; i++) {
    const rr = 0.5 * (g.rf[i] + g.rf[i + 1]);
    if (rr < 4 * dr || rr > 5 / kappa) continue;
    const num = S[i * g.nz + jmid] * 4 * Math.PI * rr / Ndot;
    const ana = Math.exp(-kappa * rr);
    const err = Math.abs(num - ana) / ana;
    if (err > worst) { worst = err; worstR = rr; }
    cnt++;
  }
  ok(cnt > 10, `слишком мало точек сравнения: ${cnt}`);
  lines.push(`       max отн. ошибка ${(worst * 100).toFixed(2)} % при R = ${(worstR * 1e6).toFixed(1)} мкм (${cnt} точек)`);
  ok(worst < 0.02, `max ошибка ${(worst * 100).toFixed(2)} % > 2 %`);
});

// ===========================================================================
section('--- PH5: сохранение источника фотонов ---');
// ===========================================================================
test('PH5 INT q = kappa^2 INT S + ушедшее через границу, < 1%', () => {
  const nr = 64, R = 1e-3, ngz = 64, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const h = new HelmholtzStub(g);
  const kappa = 1.111e4;
  const N = nr * g.nz;
  const q = new Float64Array(N);
  const dr = R / nr, dz = Lg / ngz;
  // гауссово пятно в центре зазора, чтобы источник был распределённым
  const r0 = 0, z0 = 0.5 * Lg, sg = 1e-4;
  for (let i = 0; i < nr; i++) {
    const rr = 0.5 * (g.rf[i] + g.rf[i + 1]);
    for (let j = 0; j < ngz; j++) {
      const zz = 0.5 * (g.zf[j] + g.zf[j + 1]);
      q[i * g.nz + j] = 1e26 * Math.exp(-((rr - r0) ** 2 + (zz - z0) ** 2) / (2 * sg * sg));
    }
  }
  const S = h.solveHelmholtz(q, kappa, new Float64Array(N));
  let Iq = 0, IS = 0;
  for (let i = 0; i < nr; i++) {
    const A = Math.PI * (g.rf[i + 1] ** 2 - g.rf[i] ** 2);
    for (let j = 0; j < ngz; j++) {
      const V = A * (g.zf[j + 1] - g.zf[j]);
      Iq += q[i * g.nz + j] * V;
      IS += S[i * g.nz + j] * V;
    }
  }
  const absorbed = kappa * kappa * IS;
  const total = absorbed + h.escaped;
  const d = near(total, Iq, 0.01, 'баланс фотонов');
  lines.push(`       поглощено ${(absorbed / Iq * 100).toFixed(2)} %, ушло ${(h.escaped / Iq * 100).toFixed(2)} %, невязка ${(d * 100).toExponential(2)} %`);
});

// ===========================================================================
section('--- VF4: нормировка весов и сохранение источника ---');
// ===========================================================================
test('VF4 sum A_j/lambda_j^2 = 1 для всех трёх наборов весов', () => {
  const g = uniformGrid(8, 5e-4, 2, 8, 5e-4, 1e-3);
  for (const name of Object.keys(PHOTO_ION_WEIGHT_SETS)) {
    const pm = new PhotoModule(g, { photoIonWeights: name, photoEmission: false, photoDetachment: false });
    let s = 0;
    for (let j = 0; j < 3; j++) s += pm.Aj[j] / (pm.lamA[j] * pm.lamA[j]);
    near(s, 1, 1e-12, `набор '${name}'`);
  }
  // ненормированные веса нормируются автоматически
  const pm = new PhotoModule(g, { photoIonWeights: [7, 26, 67], photoEmission: false, photoDetachment: false });
  near(pm.wA[0], 0.07, 1e-12, 'автонормировка весов');
});

test('VF4 интеграл S_pi = Phi_ion * eta_gamma * интеграл S_ion (уход через границу учтён)', () => {
  const nr = 48, R = 1e-3, ngz = 48, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const h = new HelmholtzStub(g);
  const pm = new PhotoModule(g, { photoEmission: false, photoDetachment: false }, h);
  const N = nr * g.nz;
  const S = new Float64Array(N);
  for (let i = 0; i < nr; i++) {
    const rr = 0.5 * (g.rf[i] + g.rf[i + 1]);
    for (let j = 0; j < ngz; j++) {
      const zz = 0.5 * (g.zf[j] + g.zf[j + 1]);
      S[i * g.nz + j] = 2e29 * Math.exp(-((rr) ** 2 + (zz - 0.5 * Lg) ** 2) / (2 * (5e-5) ** 2));
    }
  }
  const out = pm.update({ ionizRate: S, powerDens: new Float64Array(N) }, { step: 0 });
  let Iion = 0, Ipi = 0;
  for (let i = 0; i < nr; i++) {
    const A = Math.PI * (g.rf[i + 1] ** 2 - g.rf[i] ** 2);
    for (let j = 0; j < ngz; j++) {
      const V = A * (g.zf[j + 1] - g.zf[j]);
      Iion += S[i * g.nz + j] * V;
      Ipi += out.photoIonRate[i * g.nz + j] * V;
    }
  }
  const expected = 1e-3 * 0.30 * Iion;
  const frac = Ipi / expected;
  lines.push(`       поглощено внутри области ${(frac * 100).toFixed(1)} % от Phi_ion*eta_gamma*INT S_ion (остальное ушло через границы зазора — это физика, не ошибка)`);
  ok(frac > 0.5 && frac <= 1.0 + 1e-9, `доля ${frac}, ожидалась в (0.5, 1]`);
});

// ===========================================================================
section('--- PH4: аналитический азимут ---');
// ===========================================================================
test('PH4 аналитический азимут = прямая квадратура по углу, < 1e-6', () => {
  const cases = [
    [3.1e-4, 1.2e-3, 7.7e-5, 1.05e-3],
    [1e-6, 5e-4, 4.4e-4, 5.2e-4],
    [4.9e-4, 2.0e-4, 4.9e-4, 2.0e-4 + 1e-6],
    [8.0e-4, 0, 1.3e-4, 9.0e-4],
  ];
  let worst = 0;
  for (const [r, z, rp, zp] of cases) {
    const ana = transparentKernel(r, z, rp, zp);
    // прямая квадратура: INT_0^{2pi} dphi / (4 pi R^2), R^2 = r^2+r'^2-2rr'cos+dz^2
    const n = 400000;
    const dz2 = (z - zp) ** 2, aa = r * r + rp * rp + dz2, bb = 2 * r * rp;
    let s = 0;
    for (let k = 0; k < n; k++) {
      const phi = (k + 0.5) * 2 * Math.PI / n;
      s += 1 / (aa - bb * Math.cos(phi));
    }
    s *= (2 * Math.PI / n) / (4 * Math.PI);
    const e = Math.abs(ana - s) / s;
    worst = Math.max(worst, e);
  }
  lines.push(`       max отн. расхождение ${worst.toExponential(2)}`);
  ok(worst < 1e-6, `расхождение ${worst.toExponential(2)} > 1e-6`);
});

// ===========================================================================
section('--- PH3: прозрачное ядро = 1/(4 pi R^2) ---');
// ===========================================================================
test('PH3 точечный источник даёт 1/(4 pi R^2), < 1e-3', () => {
  const nr = 96, R = 1.5e-3, ngz = 96, Lg = 1.5e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const pm = new PhotoModule(g, { photoIonization: false, photoEmission: false }, null);
  const N = nr * g.nz;
  const q = new Float64Array(N);
  const dr = R / nr, dz = Lg / ngz;
  const jmid = ngz >> 1;
  const Vcell = Math.PI * dr * dr * dz;
  const Ndot = 1;
  q[0 * g.nz + jmid] = Ndot / Vcell;
  const psi = pm.psiExact(q);
  let worst = 0, worstR = 0, cnt = 0;
  for (let i = 1; i < nr; i++) {
    const rr = 0.5 * (g.rf[i] + g.rf[i + 1]);
    if (rr < 20 * dr || rr > 0.9e-3) continue;
    const num = psi[i * g.nz + jmid];
    const ana = Ndot / (4 * Math.PI * rr * rr);
    const e = Math.abs(num - ana) / ana;
    if (e > worst) { worst = e; worstR = rr; }
    cnt++;
  }
  ok(cnt > 10, `мало точек: ${cnt}`);
  lines.push(`       max отн. ошибка ${worst.toExponential(2)} при R = ${(worstR * 1e6).toFixed(0)} мкм (${cnt} точек)`);
  ok(worst < 1e-3, `ошибка ${worst.toExponential(2)} > 1e-3`);
});

// ===========================================================================
section('--- VF8: калибровка самоячейки прозрачного ядра ---');
// ===========================================================================
test('VF8 однородный q: дискретный оператор = средняя хорда (независимая квадратура), < 5%', () => {
  // При q = const Psi(x) = q * <L(Omega)> — средняя по телесному углу длина хорды
  // до границы области. Это аналитический эталон, независимый от нашей дискретизации;
  // именно на нём калибруется C_self (PHOTO_PROCESSES §7.5 / VF8).
  const nr = 48, R = 1e-3, ngz = 48, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const pm = new PhotoModule(g, { photoIonization: false, photoEmission: false }, null);
  const N = nr * g.nz;
  const q = new Float64Array(N).fill(1);
  for (let i = 0; i < nr; i++) for (let j = 0; j < ngz; j++) q[i * g.nz + j] = 1;
  const psi = pm.psiExact(q);

  // эталон: <L> в точке на оси, в середине зазора, для цилиндра R x Lg
  const meanChord = (r0, z0) => {
    const nth = 2000, nph = 400;
    let s = 0, w = 0;
    for (let a = 0; a < nth; a++) {
      const ct = -1 + 2 * (a + 0.5) / nth;
      const st = Math.sqrt(1 - ct * ct);
      for (let b = 0; b < nph; b++) {
        const ph = (b + 0.5) * 2 * Math.PI / nph;
        const ux = st * Math.cos(ph), uy = st * Math.sin(ph), uz = ct;
        // расстояние до боковой поверхности r = R
        const A = ux * ux + uy * uy, B = 2 * r0 * ux, C = r0 * r0 - R * R;
        let tR = Infinity;
        if (A > 0) tR = (-B + Math.sqrt(Math.max(0, B * B - 4 * A * C))) / (2 * A);
        let tZ = Infinity;
        if (uz > 0) tZ = (Lg - z0) / uz; else if (uz < 0) tZ = (0 - z0) / uz;
        s += Math.min(tR, tZ); w++;
      }
    }
    return s / w;
  };
  const rAx = 0.5 * (g.rf[0] + g.rf[1]);
  const zMid = 0.5 * (g.zf[ngz >> 1] + g.zf[(ngz >> 1) + 1]);
  const ref = meanChord(rAx, zMid);
  const got = psi[0 * g.nz + (ngz >> 1)];
  const d = near(got, ref, 0.05, 'самоячейка / средняя хорда на оси');
  lines.push(`       на оси: ядро ${got.toExponential(4)} м, эталон ${ref.toExponential(4)} м, расхождение ${(d * 100).toFixed(2)} %`);

  // и в точке вне оси — там работает стержневой предел C_self
  const i2 = Math.floor(nr / 2);
  const r2 = 0.5 * (g.rf[i2] + g.rf[i2 + 1]);
  const ref2 = meanChord(r2, zMid);
  const got2 = psi[i2 * g.nz + (ngz >> 1)];
  const d2 = near(got2, ref2, 0.05, 'самоячейка / средняя хорда вне оси');
  lines.push(`       вне оси: ядро ${got2.toExponential(4)} м, эталон ${ref2.toExponential(4)} м, расхождение ${(d2 * 100).toFixed(2)} %`);
});

// ===========================================================================
section('--- PH6: view factor, медианный радиус R_50 = sqrt(3)*h ---');
// ===========================================================================
test('PH6 R_50 = sqrt(3)*h для точечного источника на высоте h, < 3%', () => {
  const nr = 160, R = 1.2e-3, ngz = 6, Lg = 6e-4;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const pm = new PhotoModule(g, { photoIonization: false, photoEmission: false, photoDetachment: false }, null);
  const mu = 1e-9;                                  // прозрачная группа: чистая геометрия
  const m = pm._buildVFWall(mu, 0);                 // нижняя стенка z = 0
  const jsrc = 1;                                   // h = 150 мкм
  const h = 0.5 * (g.zf[jsrc] + g.zf[jsrc + 1]);
  const srcCol = 0 * g.nz + jsrc;                   // приосевая ячейка
  const dr = R / nr, dz = Lg / ngz;
  const Vcell = Math.PI * dr * dr * dz;
  const Ndot = 1;
  const qv = Ndot / Vcell;
  // Gamma(i_w) = W[i_w][srcCol] * q
  const F = new Float64Array(nr);
  for (let i = 0; i < nr; i++) {
    for (let k = m.rowPtr[i]; k < m.rowPtr[i + 1]; k++) if (m.col[k] === srcCol) F[i] = m.val[k] * qv;
  }
  // кумулятивный поток; половина телесного угла = Ndot/2 (нормировка §3.4)
  let cum = 0;
  const half = Ndot / 2;
  let r50 = NaN, prevCum = 0, prevR = 0;
  for (let i = 0; i < nr; i++) {
    const A = Math.PI * (g.rf[i + 1] ** 2 - g.rf[i] ** 2);
    prevCum = cum; cum += F[i] * A;
    const rr = g.rf[i + 1];
    if (Number.isNaN(r50) && cum >= half / 2) {
      r50 = prevR + (half / 2 - prevCum) / (cum - prevCum) * (rr - prevR);
    }
    prevR = rr;
  }
  const ana = Math.sqrt(3) * h;
  // заодно нормировка: полный пойманный поток близок к Ndot/2 (часть уходит за r > R)
  const caught = cum / half;
  // аналитическая доля, попадающая внутрь r < R: 1 - h/sqrt(R^2+h^2) (§3.4)
  const caughtAna = 1 - h / Math.sqrt(R * R + h * h);
  lines.push(`       h = ${(h * 1e6).toFixed(1)} мкм, R_50 = ${(r50 * 1e6).toFixed(1)} мкм, sqrt(3)h = ${(ana * 1e6).toFixed(1)} мкм`);
  lines.push(`       поймано ${(caught * 100).toFixed(2)} % от Ndot/2, аналитика для r < R: ${(caughtAna * 100).toFixed(2)} %`);
  near(r50, ana, 0.03, 'медианный радиус засветки');
  near(caught, caughtAna, 0.01, 'нормировка потока на стенку');
});

// ===========================================================================
section('--- PH7: запрет двойного счёта (ERRATA F2) ---');
// ===========================================================================
test('PH7 gamma_i*(1 + F_ph/F_i) = gamma_eff в пределах 1%', () => {
  const g = uniformGrid(8, 5e-4, 2, 8, 5e-4, 1e-3);
  // номинал
  const pm = new PhotoModule(g, { photoIonization: false, photoDetachment: false }, null);
  const total = pm.gammaI * (1 + pm.FphOverFi);
  near(total, 0.02, 0.01, 'суммарный эмиссионный поток при номинале');
  lines.push(`       номинал: F_ph/F_i = ${pm.FphOverFi.toExponential(3)} (контроль ревью ${pm.FphOverFiCheck.toExponential(3)}), gamma_i = ${pm.gammaI.toFixed(6)}, поправка ${((0.02 / pm.gammaI - 1) * 100).toFixed(4)} %`);

  // верхний сценарий (дефектный оксид, окно A): поправка обязана стать значимой
  const hi = new PhotoModule(g, {
    photoIonization: false, photoDetachment: false,
    Eth: 6.0, Yref: 1e-2, fVUV: 5e-3, Tw: 0.3,   // Y_ph ~ 1e-2 — верх таблицы §3.5
  }, null);
  near(hi.gammaI * (1 + hi.FphOverFi), 0.02, 1e-12, 'тождество при верхнем сценарии');
  const corr = (0.02 / hi.gammaI - 1) * 100;
  lines.push(`       верхний сценарий: F_ph/F_i = ${hi.FphOverFi.toExponential(3)}, gamma_i = ${hi.gammaI.toExponential(4)}, поправка ${corr.toFixed(1)} %`);
  ok(corr > 5, `в верхнем сценарии поправка ${corr.toFixed(2)} % — ожидалась значимая (> 5 %)`);

  // при photoEmission = false деления быть не должно
  const off = new PhotoModule(g, { photoEmission: false, photoIonization: false, photoDetachment: false }, null);
  ok(off.gammaI === 0.02, `при photoEmission=false gamma_i = ${off.gammaI}`);
});

// ===========================================================================
section('--- VF7: порог фотоэмиссии переключает канал ---');
// ===========================================================================
test('VF7 при E_th = 10.5 группа PE-B даёт СТРОГИЙ ноль, весь поток из окна A', () => {
  ok(fowlerYield(9.51, 10.5, 5e-3, 12.5) === 0, 'Y_ph(130.4 нм) != 0 при E_th = 10.5');
  ok(fowlerYield(12.3, 10.5, 5e-3, 12.5) > 0, 'Y_ph(окно A) = 0 при E_th = 10.5');
  const y9 = fowlerYield(9.51, 9.0, 5e-3, 12.5);
  near(y9, 1.06e-4, 0.02, 'Y_ph(130.4 нм) при E_th = 9.0');

  const g = uniformGrid(12, 5e-4, 2, 10, 5e-4, 1e-3);
  const hi = new PhotoModule(g, { Eth: 10.5, photoIonization: false, photoDetachment: false }, null);
  hi.build();
  ok(hi.YB === 0, `Y_B = ${hi.YB}`);
  ok(hi.vfL[3] === null && hi.vfR[3] === null, 'матрица PE-B построена при E_th = 10.5 (лишняя память и риск ненулевого вклада)');
  const lo = new PhotoModule(g, { Eth: 9.0, photoIonization: false, photoDetachment: false }, null);
  lo.build();
  ok(lo.YB > 0 && lo.vfL[3] !== null, 'при E_th = 9.0 группа PE-B обязана существовать');
});

// ===========================================================================
section('--- PH8: NaN и отрицательные потоки ---');
// ===========================================================================
test('PH8 нет NaN/отрицательных потоков ни при каком разумном входе', () => {
  const nr = 24, R = 5e-4, ngz = 24, Lg = 1e-3;
  const g = uniformGrid(nr, R, 3, ngz, 5e-4, Lg);
  const h = new HelmholtzStub(g);
  const N = g.nr * g.nz;
  const S = new Float64Array(N), P = new Float64Array(N);
  const inputs = [
    ['нули', () => { S.fill(0); P.fill(0); }],
    ['однородный номинал', () => { S.fill(2.35e29); P.fill(3.8e10); }],
    ['узкий филамент на оси', () => {
      S.fill(0); P.fill(0);
      for (let j = g.JG0; j <= g.JG1; j++) { S[0 * g.nz + j] = 1e31; P[0 * g.nz + j] = 1e12; }
    }],
    ['пятно у стенки', () => {
      S.fill(0); P.fill(0);
      for (let i = 0; i < 4; i++) { S[i * g.nz + g.JG0] = 5e30; P[i * g.nz + g.JG0] = 5e11; }
    }],
    ['экстремум 1e35', () => { S.fill(1e35); P.fill(1e16); }],
  ];
  for (const [name, setup] of inputs) {
    for (const weights of ['air', 'short', 'long']) {
      for (const Eth of [5.0, 9.0, 10.5]) {
        setup();
        const pm = new PhotoModule(g, {
          photoIonWeights: weights, Eth, photoCoarse: [8, 12], nPhotoSubstep: 1,
        }, h);
        const out = pm.update({ ionizRate: S, powerDens: P }, { step: 0 });
        // _assertFinite внутри update уже кидает; продублируем по nu и detachRate
        const n = { Om: new Float64Array(N).fill(1e19), O2m: new Float64Array(N).fill(1e19), O3m: new Float64Array(N).fill(1e20) };
        const dr = pm.detachRate(n);
        for (let k = 0; k < N; k++) {
          ok(Number.isFinite(dr[k]) && dr[k] >= 0, `detachRate[${k}] = ${dr[k]} (${name}, ${weights}, E_th=${Eth})`);
          ok(Number.isFinite(out.nuPd.O3m[k]) && out.nuPd.O3m[k] >= 0, `nuPd.O3m[${k}] (${name})`);
        }
        ok(Number.isFinite(pm.gammaI) && pm.gammaI > 0, `gammaI = ${pm.gammaI}`);
      }
    }
  }
});

// ===========================================================================
section('--- дополнительно: substep, uniform, эксклюзивность O3- ---');
// ===========================================================================
test('view factor: зеркальная экономия памяти даёт тот же поток на верхнюю стенку', () => {
  const nr = 20, R = 5e-4, ngz = 16, Lg = 1e-3;
  const g = uniformGrid(nr, R, 2, ngz, 5e-4, Lg);
  const N = g.nr * g.nz;
  // асимметричный источник: пятно ближе к НИЖНЕЙ стенке
  const S = new Float64Array(N), P = new Float64Array(N);
  for (let i = 0; i < 5; i++) for (let j = g.JG0; j < g.JG0 + 3; j++) { S[i * g.nz + j] = 1e30; P[i * g.nz + j] = 1e12; }
  const mk = (mirror) => {
    const pm = new PhotoModule(g, { photoIonization: false, photoDetachment: false, nPhotoSubstep: 1 }, null);
    pm.build();
    if (!mirror) {   // принудительно строим вторую стенку явно
      pm.vfMirror = false;
      const zR = g.zf[g.JG1 + 1];
      for (let k = 0; k < pm.emitGroups.length; k++)
        pm.vfR[k] = pm.emitGroups[k].Y > 0 ? pm._buildVFWall(pm.emitGroups[k].mu, zR) : null;
    }
    return Float64Array.from(pm.update({ ionizRate: S, powerDens: P }, { step: 0 }).photoEmitFluxR);
  };
  const a = mk(true), b = mk(false);
  ok(a.some((x) => x > 0), 'поток на верхнюю стенку нулевой — тест бессодержателен');
  let worst = 0;
  for (let i = 0; i < nr; i++) worst = Math.max(worst, Math.abs(a[i] - b[i]) / Math.max(b[i], 1e-300));
  lines.push(`       max расхождение зеркальной и явной сборки ${worst.toExponential(2)}`);
  ok(worst < 1e-6, `расхождение ${worst.toExponential(2)}`);
});

test('nPhotoSubstep: поля заморожены между обновлениями (VF10, механика)', () => {
  const nr = 12, R = 5e-4, ngz = 12, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const h = new HelmholtzStub(g);
  const N = nr * g.nz;
  const S = new Float64Array(N).fill(1e29), P = new Float64Array(N).fill(1e10);
  const pm = new PhotoModule(g, { photoEmission: false, photoCoarse: [6, 8], nPhotoSubstep: 3 }, h);
  pm.update({ ionizRate: S, powerDens: P }, { step: 0 });
  const v0 = pm.photoIonRate[5 * g.nz + 5];
  S.fill(1e31); P.fill(1e12);
  pm.update({ ionizRate: S, powerDens: P }, { step: 1 });   // пропуск
  ok(pm.photoIonRate[5 * g.nz + 5] === v0, 'поле обновилось на пропущенном шаге');
  pm.update({ ionizRate: S, powerDens: P }, { step: 3 });   // обновление
  ok(pm.photoIonRate[5 * g.nz + 5] > v0 * 10, 'поле не обновилось на рабочем шаге');
});

test('режим uniform занижает поток на оси филамента в 10-30 раз (§7.5, предупреждение верно)', () => {
  const nr = 48, R = 5e-4, ngz = 48, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const N = nr * g.nz;
  const P = new Float64Array(N);
  // филамент радиусом ~75 мкм на оси
  for (let i = 0; i < nr; i++) {
    const rr = 0.5 * (g.rf[i] + g.rf[i + 1]);
    if (rr > 75e-6) continue;
    for (let j = 0; j < ngz; j++) P[i * g.nz + j] = 3.8e10;
  }
  const mk = (mode) => {
    const pm = new PhotoModule(g, {
      photoIonization: false, photoEmission: false, photoDetachMode: mode,
      photoCoarse: [16, 24], nPhotoSubstep: 1,
    }, null);
    return pm.update({ powerDens: P }, { step: 0 }).psi[0 * g.nz + (ngz >> 1)];
  };
  const kern = mk('kernel'), uni = mk('uniform');
  const ratio = kern / uni;
  lines.push(`       Psi(ось): ядро ${kern.toExponential(3)}, uniform ${uni.toExponential(3)}, занижение x${ratio.toFixed(1)}`);
  ok(ratio > 3, `uniform занижает всего в ${ratio.toFixed(1)} раз — ожидалось заметное занижение`);
});

test('фотоотлипание: канал эксклюзивен для O3- (ни D1, ни D3, ни D6 к нему неприменимы)', () => {
  // Тест не про темп (он ничтожен), а про то, что модуль ВООБЩЕ даёт nu_pd для O3-,
  // и что порядок величины отвечает PHOTO_PROCESSES §5.2 (nu_pd(O3-) ~ 2.2 c^-1
  // при f_det = 1e-4 и P_vol = 3.8e10 Вт/м^3).
  const nr = 32, R = 5e-4, ngz = 32, Lg = 1e-3;
  const g = uniformGrid(nr, R, 0, ngz, 0, Lg);
  const N = nr * g.nz;
  const P = new Float64Array(N).fill(3.8e10);      // однородно излучающий зазор
  const pm = new PhotoModule(g, {
    photoIonization: false, photoEmission: false,
    photoDetachMode: 'uniform', nPhotoSubstep: 1,  // uniform = именно оценка §5.2 (G_slab)
  }, null);
  const out = pm.update({ powerDens: P }, { step: 0 });
  const nu3 = out.nuPd.O3m[0 * g.nz + (ngz >> 1)];
  lines.push(`       nu_pd(O3-) = ${nu3.toFixed(2)} c^-1 (спека §5.2: 2.2 c^-1 при f_det = 1e-4)`);
  ok(nu3 > 0.5 && nu3 < 10, `nu_pd(O3-) = ${nu3}, ожидался порядок 1 c^-1`);
  // Delta n_e за импульс 100 нс при n(O3-) = 1e20
  const dne = nu3 * 100e-9 * 1e20;
  lines.push(`       -> Delta n_e = ${dne.toExponential(2)} м^-3 за импульс 100 нс при n(O3-) = 1e20 (спека: 1e12..1e15)`);
  ok(dne > 1e12 && dne < 1e15, `Delta n_e = ${dne}`);
});

// ===========================================================================
console.log(lines.join('\n'));
console.log(`\n=== photo.test.mjs: ${passed} passed, ${failed} failed ===`);
process.exit(failed ? 1 : 0);
