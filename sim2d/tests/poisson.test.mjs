// tests/poisson.test.mjs — P1..P10 (обязательные) + P11..P12 (полунеявность/PCG)
// для разделимого решателя Пуассона/Гельмгольца.
// Каждый тест печатает PASS/FAIL; при любом провале процесс завершается с кодом 1.
//
// Запуск:  node sim2d/tests/poisson.test.mjs

import { SeparableSolver, EPS0 } from '../poisson2d.mjs';

let failed = 0;
const results = [];

function check(name, ok, detail) {
  results.push({ name, ok, detail });
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}   ${detail}`);
  if (!ok) failed++;
}

const maxAbs = (a) => {
  let m = 0;
  for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i]));
  return m;
};

// ═══════════════════════════════════════════════════════════════ P1

function P1() {
  // Вакуумный конденсатор (eps_r = 1 везде), без заряда -> phi линеен по z.
  const s = new SeparableSolver({
    nr: 24, nz1: 8, ngap: 24, nz2: 8, epsR: 1,
    betaR: 1.6, betaG: 2.63, qB: 1.3147,
  });
  const U = 1e4;
  const rho = new Float64Array(s.nr * s.nz);
  const sig = new Float64Array(s.nr);
  const phi = s.solvePoisson(rho, sig, sig, U);

  let errMax = 0;
  for (let i = 0; i < s.nr; i++)
    for (let j = 0; j < s.nz; j++) {
      const ex = U * (1 - s.zc[j] / s.Lz);
      errMax = Math.max(errMax, Math.abs(phi[i * s.nz + j] - ex) / U);
    }

  // невязка собранной матрицы
  const rhs = new Float64Array(s.nr * s.nz);
  const Aphi = new Float64Array(s.nr * s.nz);
  s.buildRhsFull(rho, sig, sig, U, rhs);
  s.applyMatrix(phi, Aphi);
  let rn = 0;
  for (let i = 0; i < rhs.length; i++) rn = Math.max(rn, Math.abs(Aphi[i] - rhs[i]));
  const rrel = rn / maxAbs(rhs);

  check('P1 вакуумный конденсатор: линейный phi(z)',
    errMax < 1e-10 && rrel < 1e-10,
    `max|phi-phi_ex|/U = ${errMax.toExponential(3)}, невязка отн = ${rrel.toExponential(3)}`);
}

// ═══════════════════════════════════════════════════════════════ P2

function P2() {
  // Ёмкостный делитель: диэлектрик | газ | диэлектрик, без заряда.
  const s = new SeparableSolver({ nr: 32, nz1: 12, ngap: 40, nz2: 12, epsR: 9 });
  const U = 1e4;
  const rho = new Float64Array(s.nr * s.nz);
  const sig = new Float64Array(s.nr);
  const phi = s.solvePoisson(rho, sig, sig, U);
  const Ez = s.computeEz(phi, U);

  const dEff = s.Lg + (s.d1 + s.d2) / s.epsRd;      // 1.11111e-3 м
  const EgasEx = U / dEff;                            // 9.0e6 В/м при U=10 кВ
  const EdielEx = U / (s.epsRd * dEff);

  let eGas = 0, eDiel = 0, eR = 0;
  for (let i = 0; i < s.nr; i++) {
    for (let j = 1; j < s.nz; j++) {
      const inGas = s.isGas[j] && s.isGas[j - 1];
      const inDiel = !s.isGas[j] && !s.isGas[j - 1];
      const e = Ez[i * (s.nz + 1) + j];
      if (inGas) eGas = Math.max(eGas, Math.abs(e - EgasEx) / EgasEx);
      if (inDiel) eDiel = Math.max(eDiel, Math.abs(e - EdielEx) / EdielEx);
    }
  }
  const Er = s.computeEr(phi);
  eR = maxAbs(Er) / EgasEx;

  // ёмкости: заряд на активном электроде
  const A = Math.PI * s.R * s.R;
  let Qel = 0;
  for (let i = 0; i < s.nr; i++) Qel += s.eps[0] * Ez[i * (s.nz + 1)] * s.Acell[i];
  const Ccell = Qel / U;
  const CcellEx = (EPS0 * A) / dEff;
  const Cd = (EPS0 * s.epsRd * A) / (s.d1 + s.d2);   // последовательно два барьера
  const Cg = (EPS0 * A) / s.Lg;
  const errC = Math.abs(Ccell - CcellEx) / CcellEx;

  // сверка с ERRATA C/V4 (на 1 см²)
  const k = 1e-4 / A;
  const docOk =
    Math.abs(Cd * k * 1e12 - 7.97) < 0.01 &&
    Math.abs(Cg * k * 1e12 - 0.885) < 0.002 &&
    Math.abs(CcellEx * k * 1e12 - 0.797) < 0.002;

  check('P2 ёмкостный делитель + ёмкости',
    eGas < 1e-6 && eDiel < 1e-6 && eR < 1e-10 && errC < 1e-6 && docOk,
    `dE_gas=${eGas.toExponential(2)} dE_diel=${eDiel.toExponential(2)} Er/E=${eR.toExponential(2)} ` +
    `dC=${errC.toExponential(2)}; на 1см²: C_d=${(Cd * k * 1e12).toFixed(3)}пФ ` +
    `C_g=${(Cg * k * 1e12).toFixed(3)}пФ C_cell=${(CcellEx * k * 1e12).toFixed(3)}пФ`);
}

// ═══════════════════════════════════════════════════════════════ P3 (SOR)

/**
 * НЕЗАВИСИМЫЙ медленный SOR-решатель. Вся матрица собирается здесь заново, с нуля,
 * без единого обращения к внутренностям SeparableSolver (из него берутся только
 * геометрия сетки и eps(z), то есть исходные данные задачи).
 * aCell = eps + dt*kappa (по умолчанию — чистая электростатика).
 */
function sorSolve(s, rho, sigLo, sigHi, U, aCell = null, tolRel = 1e-12, maxSweeps = 400_000) {
  const { nr, nz, rf, dr, hr, dz, eps, Acell, JS_LO, JG0, JG1, JS_HI } = s;
  const a = aCell || (() => {
    const t = new Float64Array(nr * nz);
    for (let i = 0; i < nr; i++) for (let j = 0; j < nz; j++) t[i * nz + j] = eps[j];
    return t;
  })();

  // --- трансмиссивности граней (последовательное сложение сопротивлений полуячеек)
  const TZ = new Float64Array(nr * (nz + 1));
  const TR = new Float64Array((nr + 1) * nz);
  for (let i = 0; i < nr; i++) {
    const o = i * nz, t = i * (nz + 1);
    TZ[t] = (Acell[i] * a[o]) / (0.5 * dz[0]);
    for (let j = 1; j < nz; j++)
      TZ[t + j] = Acell[i] / ((0.5 * dz[j - 1]) / a[o + j - 1] + (0.5 * dz[j]) / a[o + j]);
    TZ[t + nz] = (Acell[i] * a[o + nz - 1]) / (0.5 * dz[nz - 1]);
  }
  for (let i = 1; i < nr; i++)
    for (let j = 0; j < nz; j++)
      TR[i * nz + j] = (2 * Math.PI * rf[i] * dz[j]) /
        ((0.5 * dr[i - 1]) / a[(i - 1) * nz + j] + (0.5 * dr[i]) / a[i * nz + j]);

  // --- правая часть (несделённая): A*phi = RHS
  const RHS = new Float64Array(nr * nz);
  for (let i = 0; i < nr; i++) {
    const o = i * nz, A = Acell[i];
    for (let j = 0; j < nz; j++) RHS[o + j] = -rho[o + j] * A * dz[j];
    const gdLo = a[o + JS_LO] / (0.5 * dz[JS_LO]), ggLo = a[o + JG0] / (0.5 * dz[JG0]);
    RHS[o + JS_LO] -= (gdLo / (gdLo + ggLo)) * sigLo[i] * A;
    RHS[o + JG0] -= (ggLo / (gdLo + ggLo)) * sigLo[i] * A;
    const ggHi = a[o + JG1] / (0.5 * dz[JG1]), gdHi = a[o + JS_HI] / (0.5 * dz[JS_HI]);
    RHS[o + JG1] -= (ggHi / (ggHi + gdHi)) * sigHi[i] * A;
    RHS[o + JS_HI] -= (gdHi / (ggHi + gdHi)) * sigHi[i] * A;
    RHS[o] -= TZ[i * (nz + 1)] * U;
  }

  const phi = new Float64Array(nr * nz);
  const omega = 1.92;
  const rhsNorm = Math.sqrt(RHS.reduce((x, v) => x + v * v, 0));
  const resid = () => {
    let s2 = 0;
    for (let i = 0; i < nr; i++) {
      const o = i * nz, t = i * (nz + 1);
      for (let j = 0; j < nz; j++) {
        const p = phi[o + j];
        let v = -(TZ[t + j] + TZ[t + j + 1] + TR[i * nz + j] + TR[(i + 1) * nz + j]) * p;
        if (j > 0) v += TZ[t + j] * phi[o + j - 1];
        if (j < nz - 1) v += TZ[t + j + 1] * phi[o + j + 1];
        if (i > 0) v += TR[i * nz + j] * phi[o - nz + j];
        if (i < nr - 1) v += TR[(i + 1) * nz + j] * phi[o + nz + j];
        const d = v - RHS[o + j];
        s2 += d * d;
      }
    }
    return Math.sqrt(s2) / rhsNorm;
  };

  let sweeps = 0, relres = Infinity;
  for (; sweeps < maxSweeps; sweeps++) {
    for (let i = 0; i < nr; i++) {
      const o = i * nz, t = i * (nz + 1);
      for (let j = 0; j < nz; j++) {
        const D = TZ[t + j] + TZ[t + j + 1] + TR[i * nz + j] + TR[(i + 1) * nz + j];
        let num = 0;
        if (j > 0) num += TZ[t + j] * phi[o + j - 1];
        if (j < nz - 1) num += TZ[t + j + 1] * phi[o + j + 1];
        if (i > 0) num += TR[i * nz + j] * phi[o - nz + j];
        if (i < nr - 1) num += TR[(i + 1) * nz + j] * phi[o + nz + j];
        phi[o + j] += omega * ((num - RHS[o + j]) / D - phi[o + j]);
      }
    }
    if (sweeps % 50 === 49) {
      relres = resid();
      if (relres < tolRel) break;
    }
  }
  return { phi, sweeps, relres };
}

function P3() {
  const s = new SeparableSolver({ nr: 20, nz1: 6, ngap: 20, nz2: 6, epsR: 9 });
  const n = s.nr * s.nz;
  // детерминированный «случайный» источник
  let seed = 12345;
  const rnd = () => {
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    return seed / 0x7fffffff - 0.5;
  };
  const rho = new Float64Array(n);
  for (let i = 0; i < s.nr; i++)
    for (let j = 0; j < s.nz; j++) if (s.isGas[j]) rho[i * s.nz + j] = rnd() * 1e-3;
  const sigLo = new Float64Array(s.nr);
  const sigHi = new Float64Array(s.nr);
  for (let i = 0; i < s.nr; i++) { sigLo[i] = rnd() * 1e-5; sigHi[i] = rnd() * 1e-5; }
  const U = 3.3e3;

  const phiFast = Float64Array.from(s.solvePoisson(rho, sigLo, sigHi, U));
  const { phi: phiSor, sweeps, relres } = sorSolve(s, rho, sigLo, sigHi, U, null, 1e-13);

  let dmax = 0;
  for (let i = 0; i < n; i++) dmax = Math.max(dmax, Math.abs(phiFast[i] - phiSor[i]));
  const rel = dmax / maxAbs(phiSor);

  check('P3 сверка с независимым SOR',
    rel < 1e-8 && relres < 1e-12,
    `отн. расхождение = ${rel.toExponential(3)} (SOR: ${sweeps} свипов, невязка ${relres.toExponential(2)})`);
}

// ═══════════════════════════════════════════════════════════════ P4

function P4() {
  // Теорема Гаусса на «таблетке», охватывающей границу газ/диэлектрик:
  //   ∮ D·dS  ==  sigma * Acell   (объёмного заряда внутри нет).
  // Индукция на гранях восстанавливается НЕЗАВИСИМО (гармоническое eps, здесь, а не из
  // решателя). Тест ловит: неверный знак sigma, w_d + w_g != 1, рассогласование сборки.
  // Заряд задан радиально НЕОДНОРОДНЫМ, поэтому радиальные потоки не нули и тоже учтены.
  const s = new SeparableSolver({ nr: 24, nz1: 8, ngap: 32, nz2: 8, epsR: 9 });
  const n = s.nr * s.nz;
  const rho = new Float64Array(n);
  const sigLo = new Float64Array(s.nr);
  const sigHi = new Float64Array(s.nr);
  for (let i = 0; i < s.nr; i++) {
    sigLo[i] = -2e-5 * Math.exp(-((s.rc[i] / 8e-5) ** 2));
    sigHi[i] = 1.3e-5 * Math.exp(-((s.rc[i] / 1.5e-4) ** 2));
  }
  const U = 5e3;
  const phi = s.solvePoisson(rho, sigLo, sigHi, U);
  const { nr, nz, JS_LO, JG0, JG1, JS_HI, eps, hz, dz, dr, hr, rf, Acell } = s;

  // независимые трансмиссивности граней
  const tz = new Float64Array(nz + 1);
  for (let j = 1; j < nz; j++) tz[j] = 1 / ((0.5 * dz[j - 1]) / eps[j - 1] + (0.5 * dz[j]) / eps[j]);
  tz[0] = eps[0] / (0.5 * dz[0]);
  tz[nz] = eps[nz - 1] / (0.5 * dz[nz - 1]);
  const trF = new Float64Array(nr + 1);
  for (let i = 1; i < nr; i++) trF[i] = (2 * Math.PI * rf[i]) / hr[i];

  const P = (i, j) => (i < 0 || i >= nr ? 0 : j < 0 ? U : j >= nz ? 0 : phi[i * nz + j]);
  // поток D наружу через «таблетку» из двух ячеек (i, ja) и (i, jb=ja+1)
  const boxFlux = (i, ja, jb) => {
    let f = 0;
    f += -tz[ja] * Acell[i] * (P(i, ja - 1) - P(i, ja));       // низ, n = -z
    f += tz[jb + 1] * Acell[i] * (P(i, jb) - P(i, jb + 1));    // верх, n = +z
    for (const j of [ja, jb]) {
      f += -trF[i] * eps[j] * dz[j] * (P(i - 1, j) - P(i, j)); // внутренняя r-грань, n = -r
      f += trF[i + 1] * eps[j] * dz[j] * (P(i, j) - P(i + 1, j)); // внешняя r-грань, n = +r
    }
    return f;
  };

  let errLo = 0, errHi = 0, scale = 0;
  for (let i = 0; i < nr; i++) {
    errLo = Math.max(errLo, Math.abs(boxFlux(i, JS_LO, JG0) - sigLo[i] * Acell[i]));
    errHi = Math.max(errHi, Math.abs(boxFlux(i, JG1, JS_HI) - sigHi[i] * Acell[i]));
    scale = Math.max(scale, Math.abs(sigLo[i]) * Acell[i], Math.abs(sigHi[i]) * Acell[i]);
  }
  const rel = Math.max(errLo, errHi) / scale;
  check('P4 скачок D_n == sigma (закон Гаусса на таблетке)', rel < 1e-8,
    `max|∮D·dS - sigma*A| / max|sigma*A| = ${rel.toExponential(3)}`);
}

// ═══════════════════════════════════════════════════════════════ P5 (ЗНАК sigma)

function P5() {
  // ERRATA A1 + C/V8. Разряд при U>0 (анод z=0): электроны уходят к НИЖНЕЙ поверхности
  // (sigma_lo < 0), положительные ионы — к ВЕРХНЕЙ (sigma_hi > 0).
  // При U=0 остаточный заряд обязан создавать поле, НАПРАВЛЕННОЕ ПРОТИВ приложенного,
  // с |E_mem| = 0.1 * sigma/eps0.
  const s = new SeparableSolver({ nr: 24, nz1: 10, ngap: 36, nz2: 10, epsR: 9 });
  const n = s.nr * s.nz;
  const rho = new Float64Array(n);
  const sg = 1e-5; // Кл/м^2
  const sigLo = new Float64Array(s.nr).fill(-sg);
  const sigHi = new Float64Array(s.nr).fill(+sg);
  const phi = s.solvePoisson(rho, sigLo, sigHi, 0);
  const Ez = s.computeEz(phi, 0);

  // ожидание: E_z в зазоре = -0.1*sigma/eps0 (экранирует положительное приложенное поле)
  const Eex = -0.1 * sg / EPS0; // -5.2e6*(sg/1e-5)... при sg=1e-5: -1.1295e5? проверяем численно
  let eMax = 0, eGapMin = Infinity, eGapMax = -Infinity;
  for (let i = 0; i < s.nr; i++)
    for (let j = s.JG0 + 1; j <= s.JG1; j++) {
      const e = Ez[i * (s.nz + 1) + j];
      eGapMin = Math.min(eGapMin, e);
      eGapMax = Math.max(eGapMax, e);
      eMax = Math.max(eMax, Math.abs(e - Eex) / Math.abs(Eex));
    }
  // знак: приложенное U>0 даёт E_z>0, значит память обязана быть отрицательной
  const screening = eGapMax < 0;

  // и контрольный физический смысл: с этим sigma напряжение на зазоре ПАДАЕТ
  const U = 5e3;
  const phi0 = Float64Array.from(s.solvePoisson(rho, new Float64Array(s.nr), new Float64Array(s.nr), U));
  const phiS = s.solvePoisson(rho, sigLo, sigHi, U);
  const ug0 = phi0[0 * s.nz + s.JG0] - phi0[0 * s.nz + s.JG1];
  const ugS = phiS[0 * s.nz + s.JG0] - phiS[0 * s.nz + s.JG1];
  const damps = Math.abs(ugS) < Math.abs(ug0);

  check('P5 ЗНАК sigma: поле памяти экранирует, |E_mem| = 0.1*sigma/eps0',
    eMax < 0.01 && screening && damps,
    `E_gap = ${eGapMin.toExponential(4)}..${eGapMax.toExponential(4)} В/м, ожидание ${Eex.toExponential(4)}, ` +
    `отн.ошибка ${(eMax * 100).toFixed(4)}%, U_gap: ${ug0.toFixed(1)} -> ${ugS.toFixed(1)} В`);
}

// ═══════════════════════════════════════════════════════════════ P6

function P6() {
  const s = new SeparableSolver({ nr: 64, nz1: 8, ngap: 16, nz2: 8 });
  const nr = s.nr;

  // --- K_r строится ЗАНОВО из определения K_r = diag(Acell) * L_r
  const Kr = new Float64Array(nr * nr);
  const u = new Float64Array(nr);
  for (let c = 0; c < nr; c++) {
    u.fill(0); u[c] = 1;
    for (let i = 0; i < nr; i++) {
      const fp = i + 1 < nr ? (s.rf[i + 1] * (u[i + 1] - u[i])) / s.hr[i + 1] : 0;
      const fm = i > 0 ? (s.rf[i] * (u[i] - u[i - 1])) / s.hr[i] : 0; // rf[0]=0 => на оси нуля не делим
      Kr[i * nr + c] = 2 * Math.PI * (fp - fm);
    }
  }
  let asym = 0, rowsum = 0, knorm = 0;
  for (let i = 0; i < nr; i++) {
    let rs = 0;
    for (let j = 0; j < nr; j++) {
      asym = Math.max(asym, Math.abs(Kr[i * nr + j] - Kr[j * nr + i]));
      knorm = Math.max(knorm, Math.abs(Kr[i * nr + j]));
      rs += Kr[i * nr + j];
    }
    rowsum = Math.max(rowsum, Math.abs(rs));
  }

  // --- ортонормировка в метрике W = diag(Acell) (= 2*pi * M_r из ERRATA B1)
  let orth = 0;
  const W = s.Acell;
  for (let a = 0; a < nr; a++)
    for (let b = 0; b < nr; b++) {
      let v = 0;
      for (let i = 0; i < nr; i++) v += s.vr[i * nr + a] * W[i] * s.vr[i * nr + b];
      orth = Math.max(orth, Math.abs(v - (a === b ? 1 : 0)));
    }

  // --- собственная задача: K_r v = -lam W v
  let eig = 0;
  for (let k = 0; k < nr; k++) {
    for (let i = 0; i < nr; i++) {
      let kv = 0;
      for (let j = 0; j < nr; j++) kv += Kr[i * nr + j] * s.vr[j * nr + k];
      eig = Math.max(eig, Math.abs(kv + s.lam[k] * W[i] * s.vr[i * nr + k]));
    }
  }
  const eigRel = eig / (knorm * maxAbs(s.vr));

  check('P6 ортонормировка / симметрия K_r / K_r*1 = 0',
    orth < 1e-12 && asym === 0 && rowsum < 1e-9 * knorm && eigRel < 1e-11,
    `max|v^T W v - I| = ${orth.toExponential(2)}, асимметрия = ${asym}, ` +
    `max|K_r*1|/||K_r|| = ${(rowsum / knorm).toExponential(2)}, невязка задачи = ${eigRel.toExponential(2)}, ` +
    `max|c'| = ${s.maxCprime.toFixed(6)}`);
}

// ═══════════════════════════════════════════════════════════════ P7

function P7() {
  // Neumann-Neumann в цилиндре: k_n*R — нули J_1. lam_n = (j_{1,n}/R)^2.
  const J1z = [3.8317059702, 7.0155866698, 10.1734681351, 13.3236919363, 16.4706300509];
  const s = new SeparableSolver({ nr: 192, nz1: 4, ngap: 8, nz2: 4, betaR: 1.6 });
  let worst = 0;
  const got = [];
  for (let m = 0; m < 5; m++) {
    const ex = (J1z[m] / s.R) ** 2;
    const num = s.lam[m + 1];
    got.push(Math.sqrt(num) * s.R);
    worst = Math.max(worst, Math.abs(num - ex) / ex);
  }
  check('P7 собственные значения = нули J_1',
    worst < 0.01,
    `k*R числ. = [${got.map((x) => x.toFixed(4)).join(', ')}], ` +
    `аналит. = [${J1z.map((x) => x.toFixed(4)).join(', ')}], max отн.ошибка = ${(worst * 100).toFixed(3)}%`);
}

// ═══════════════════════════════════════════════════════════════ P8

function P8() {
  // Точечный источник в бесконечной поглощающей среде: S = Q*exp(-kappa*R)/(4*pi*R).
  // Равномерная сетка, зазор 1 мм, kappa = 1e4 1/м (ell = 100 мкм).
  const s = new SeparableSolver({
    nr: 128, nz1: 4, ngap: 256, nz2: 4, betaR: 0, betaG: 0, qB: 1,
  });
  const kappa = 1e4;
  const n = s.nr * s.nz;
  const q = new Float64Array(n);
  const jc = s.JG0 + (s.ngap >> 1);
  const Qtot = 1.0;
  q[0 * s.nz + jc] = Qtot / (s.Acell[0] * s.dz[jc]);
  const S = s.solveHelmholtz(q, kappa);

  const zsrc = s.zc[jc];
  const anal = (R) => (Qtot * Math.exp(-kappa * R)) / (4 * Math.PI * R);
  let worstZ = 0, worstR = 0;
  const sample = [];
  for (let j = jc + 1; j < s.JG1; j++) {
    const R = s.zc[j] - zsrc;
    if (R < 80e-6 || R > 250e-6) continue;
    const num = S[0 * s.nz + j];
    const e = Math.abs(num - anal(R)) / anal(R);
    worstZ = Math.max(worstZ, e);
    if (sample.length < 3) sample.push(`z:R=${(R * 1e6).toFixed(0)}мкм ${(e * 100).toFixed(2)}%`);
  }
  for (let i = 1; i < s.nr; i++) {
    const R = s.rc[i];
    if (R < 80e-6 || R > 250e-6) continue;
    const num = S[i * s.nz + jc];
    const e = Math.abs(num - anal(R)) / anal(R);
    worstR = Math.max(worstR, e);
  }
  check('P8 Гельмгольц: точечный источник ~ exp(-kR)/(4piR)',
    worstZ < 0.01 && worstR < 0.01,
    `max отн.ошибка: по z ${(worstZ * 100).toFixed(3)}%, по r ${(worstR * 100).toFixed(3)}% ` +
    `(R = 80..250 мкм); ${sample.join(', ')}`);
}

// ═══════════════════════════════════════════════════════════════ P9

function P9() {
  // При kappa=0 и Дирихле на полной области Гельмгольц обязан воспроизвести
  // решение Пуассона с eps = eps0 (коэффициент 1 против eps0 => q = rho/eps0).
  const s = new SeparableSolver({ nr: 24, nz1: 8, ngap: 24, nz2: 8, epsR: 1 });
  const n = s.nr * s.nz;
  let seed = 999;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };
  const rho = new Float64Array(n);
  for (let i = 0; i < n; i++) rho[i] = rnd() * 1e-3;
  const zs = new Float64Array(s.nr);
  const phi = Float64Array.from(s.solvePoisson(rho, zs, zs, 0));
  const q = new Float64Array(n);
  for (let i = 0; i < n; i++) q[i] = rho[i] / EPS0;
  const S = s.solveHelmholtz(q, 0, { domain: 'full', bc: 'dirichlet' });
  let d = 0;
  for (let i = 0; i < n; i++) d = Math.max(d, Math.abs(S[i] - phi[i]));
  const rel = d / maxAbs(phi);
  check('P9 Гельмгольц(kappa=0) == Пуассон(eps=1)', rel < 1e-10,
    `отн. расхождение = ${rel.toExponential(3)}`);
}

// ═══════════════════════════════════════════════════════════════ P10

function P10() {
  // Весовой потенциал (ERRATA A3): считается БЕЗ проводимости плазмы, один раз.
  // В нашей геометрии он ТОЧНО одномерен (NUMERICS_2D §7.1).
  const s = new SeparableSolver({ nr: 48, nz1: 12, ngap: 48, nz2: 12, epsR: 9 });
  const { nr, nz } = s;
  let radial = 0, scale = 0;
  for (let j = 0; j < nz; j++) {
    const ref = s.phiL[0 * nz + j];
    scale = Math.max(scale, Math.abs(ref));
    for (let i = 1; i < nr; i++) radial = Math.max(radial, Math.abs(s.phiL[i * nz + j] - ref));
  }
  const radRel = radial / scale;

  const dEff = s.Lg + (s.d1 + s.d2) / s.epsRd;
  const EgasEx = 1 / dEff, EdielEx = 1 / (s.epsRd * dEff);
  const Ez = s.computeEz(s.phiL, 1.0);
  let eG = 0, eD = 0;
  for (let i = 0; i < nr; i++)
    for (let j = 1; j < nz; j++) {
      const e = Ez[i * (nz + 1) + j];
      if (s.isGas[j] && s.isGas[j - 1]) eG = Math.max(eG, Math.abs(e - EgasEx) / EgasEx);
      if (!s.isGas[j] && !s.isGas[j - 1]) eD = Math.max(eD, Math.abs(e - EdielEx) / EdielEx);
    }
  // аналитический профиль phi_L(z)
  let errPhi = 0;
  for (let j = 0; j < nz; j++) {
    let acc = 0;
    for (let jj = 0; jj < j; jj++) acc += (s.isGas[jj] ? EgasEx : EdielEx) * s.dz[jj];
    acc += (s.isGas[j] ? EgasEx : EdielEx) * 0.5 * s.dz[j];
    errPhi = Math.max(errPhi, Math.abs(s.phiL[0 * nz + j] - (1 - acc)));
  }
  const ErL = maxAbs(s.computeEr(s.phiL)) / EgasEx;

  check('P10 весовой потенциал phi_L: строго одномерен + аналитика делителя',
    radRel < 1e-12 && ErL < 1e-10 && eG < 1e-10 && eD < 1e-10 && errPhi < 1e-12,
    `радиальная зависимость = ${radRel.toExponential(2)}, max|E_r|/E_gas = ${ErL.toExponential(2)}, ` +
    `E_gas = ${EgasEx.toFixed(4)} (ош. ${eG.toExponential(2)}), E_diel = ${EdielEx.toFixed(4)} ` +
    `(ош. ${eD.toExponential(2)}), max|dphi_L| = ${errPhi.toExponential(2)}, d_eff = ${dEff.toExponential(6)} м`);
}

// ═══════════════════════════════════════════ P11/P12 — полунеявная поправка + PCG
// (сверх обязательного списка: путь solveSemiImplicit иначе остался бы непроверенным)

function makeCase(s) {
  let seed = 4242;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };
  const rho = new Float64Array(s.nr * s.nz);
  for (let i = 0; i < s.nr; i++)
    for (let j = 0; j < s.nz; j++) if (s.isGas[j]) rho[i * s.nz + j] = rnd() * 1e-3;
  const sigLo = new Float64Array(s.nr);
  const sigHi = new Float64Array(s.nr);
  for (let i = 0; i < s.nr; i++) { sigLo[i] = -1e-5 * Math.exp(-((s.rc[i] / 1e-4) ** 2)); sigHi[i] = 0.7e-5 * rnd(); }
  return { rho, sigLo, sigHi, U: 4e3 };
}

function P11() {
  // При dt*kappa = 0 полунеявный решатель обязан ТОЧНО воспроизвести прямой.
  const s = new SeparableSolver({ nr: 20, nz1: 6, ngap: 20, nz2: 6, epsR: 9 });
  const c = makeCase(s);
  const direct = Float64Array.from(s.solvePoisson(c.rho, c.sigLo, c.sigHi, c.U));
  const { phi, iters } = s.solveSemiImplicit({ ...c, dt: 0, kappaCell: null, relTol: 1e-13, maxIter: 20 });
  let d = 0;
  for (let i = 0; i < phi.length; i++) d = Math.max(d, Math.abs(phi[i] - direct[i]));
  s.resetElectrostatic();
  check('P11 полунеявный решатель при dt*kappa=0 == прямой', d / maxAbs(direct) < 1e-12 && iters <= 1,
    `отн. расхождение = ${(d / maxAbs(direct)).toExponential(3)}, итераций CG = ${iters}`);
}

function P12() {
  // Филамент: радиальный контраст dt*kappa/eps0 ~ 300 на оси и ~0 на периферии.
  // PCG обязан сойтись к решению НЕЗАВИСИМОГО SOR с той же матрицей.
  const s = new SeparableSolver({ nr: 20, nz1: 6, ngap: 20, nz2: 6, epsR: 9 });
  const c = makeCase(s);
  const dt = 1e-11;
  const kap = new Float64Array(s.nr * s.nz);
  for (let i = 0; i < s.nr; i++)
    for (let j = 0; j < s.nz; j++)
      if (s.isGas[j]) kap[i * s.nz + j] = (300 * EPS0 / dt) * Math.exp(-((s.rc[i] / 6e-5) ** 2));
  const aCell = new Float64Array(s.nr * s.nz);
  for (let i = 0; i < s.nr; i++)
    for (let j = 0; j < s.nz; j++) aCell[i * s.nz + j] = s.eps[j] + dt * kap[i * s.nz + j];

  const { phi, iters, resid } = s.solveSemiImplicit({ ...c, dt, kappaCell: kap, relTol: 1e-13, maxIter: 200 });
  const pcg = Float64Array.from(phi);
  const { phi: ref, relres } = sorSolve(s, c.rho, c.sigLo, c.sigHi, c.U, aCell, 1e-13);
  s.resetElectrostatic();
  let d = 0;
  for (let i = 0; i < pcg.length; i++) d = Math.max(d, Math.abs(pcg[i] - ref[i]));
  const rel = d / maxAbs(ref);
  check('P12 PCG с филаментом (контраст ~3e2) vs независимый SOR',
    rel < 1e-8 && relres < 1e-12,
    `отн. расхождение = ${rel.toExponential(3)}, итераций CG = ${iters}, невязка PCG = ${resid.toExponential(2)}`);
}

// ═══════════════════════════════════════════════════════════════

const t0 = Date.now();
P1(); P2(); P3(); P4(); P5(); P6(); P7(); P8(); P9(); P10(); P11(); P12();
console.log(`\n${results.length - failed}/${results.length} тестов пройдено за ${((Date.now() - t0) / 1000).toFixed(1)} с`);
if (failed) {
  console.log(`ПРОВАЛЕНО: ${results.filter((r) => !r.ok).map((r) => r.name).join('; ')}`);
  process.exit(1);
}
