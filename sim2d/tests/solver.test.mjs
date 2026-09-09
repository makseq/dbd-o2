// tests/solver.test.mjs — тесты главного солвера DBD2D (S1..S8).
//
// S3 и S6 — БЛОКЕРЫ (ERRATA A5 и A1). Их допуски не ослабляются ни при каких условиях.
// Запуск: node sim2d/tests/solver.test.mjs

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { DBD2D, bern } from '../solver2d.mjs';
import { SPECIES_IDS, QE, ME, meanEnergy, bernoulli } from '../physics2d.mjs';
import { Recorder } from '../recorder.mjs';

let passed = 0, failed = 0;
const t0 = Date.now();

function test(name, fn) {
  try {
    const msg = fn();
    passed++;
    console.log(`PASS  ${name}${msg ? '   ' + msg : ''}`);
  } catch (e) {
    failed++;
    console.log(`FAIL  ${name}\n      ${e.message}`);
    if (process.env.TRACE) console.log(e.stack);
  }
}
function ok(cond, msg) { if (!cond) throw new Error(msg); }
function fmt(x, d = 3) { return Number(x).toExponential(d); }

// ─────────────────────────────────────────────────────────── вспомогательное

/** Полная «масса» сорта в газе: sum n*Acell*dz. */
function mass(s, id) {
  const P = s.poisson, a = s.n[id];
  let m = 0;
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) m += a[i * s.nz + j] * P.Acell[i] * P.dz[j];
  }
  return m;
}
/** Центр масс по z. */
function centerZ(s, id) {
  const P = s.poisson, a = s.n[id];
  let m = 0, mz = 0;
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) {
      const w = a[i * s.nz + j] * P.Acell[i] * P.dz[j];
      m += w; mz += w * P.zc[j];
    }
  }
  return mz / m;
}
function clearAll(s) {
  for (const id of SPECIES_IDS) s.n[id].fill(0);
}

// ─────────────────────────────────────────────────────────── S0: функция Бернулли

test('S0  B(x)=x/expm1(x): ряд Тейлора, пределы, тождество B(-x)=B(x)+x', () => {
  ok(Math.abs(bern(0) - 1) < 1e-15, `B(0) = ${bern(0)}`);
  let worst = 0, worstX = 0;
  for (const x of [1e-12, 1e-8, 1e-5, 1e-4, 1e-3, 0.1, 1, 5, 20, 39, 41, 100, 700, 1e6]) {
    for (const sgn of [1, -1]) {
      const v = bern(sgn * x);
      ok(Number.isFinite(v), `B(${sgn * x}) = ${v}`);
      // B(x) -> 0 при x -> +inf (при x > 745 это машинный ноль — это ПРАВИЛЬНО),
      // B(x) -> -x при x -> -inf. Ни NaN, ни Infinity, ни отрицательных значений.
      ok(v >= 0, `B(${sgn * x}) = ${v} < 0`);
      if (sgn * x < 700) ok(v > 0, `B(${sgn * x}) = ${v}: преждевременный ноль`);
      // тождество
      const id = Math.abs(bern(-sgn * x) - (bern(sgn * x) + sgn * x));
      const rel = id / Math.max(1, Math.abs(bern(sgn * x)));
      if (rel > worst) { worst = rel; worstX = sgn * x; }
    }
  }
  ok(worst < 1e-12, `тождество B(-x)=B(x)+x нарушено при x=${worstX}: ${fmt(worst)}`);
  // сверка с эталоном physics2d на «безопасном» диапазоне
  let dmax = 0;
  for (let k = -200; k <= 200; k++) {
    const x = k * 0.15;
    const a = bern(x), b = bernoulli(x);
    dmax = Math.max(dmax, Math.abs(a - b) / Math.max(1e-300, Math.abs(a)));
  }
  ok(dmax < 1e-12, `расхождение с physics2d.bernoulli = ${fmt(dmax)}`);
  return `max|B(-x)-B(x)-x|отн = ${fmt(worst)}, расхождение с physics2d = ${fmt(dmax)}`;
});

// ─────────────────────────────────────────────────────────── S1: чистый дрейф по z

test('S1  чистый дрейф гауссова сгустка по z: масса < 1e-10, центр = v*t < 2%', () => {
  const Ez = 1e6;
  const s = new DBD2D({
    nr: 8, nz: 40, nzDiel: 4, transportOnly: true, frozenEz: Ez,
    betaR: 0, betaG: 0, qB: 1, neutralTransport: false,
    seedBackground: 0, seedSpotAmp: 0, dtMax: 1e-7,
  });
  clearAll(s);
  const P = s.poisson;
  const zc0 = 0.5 * (P.zf[s.JG0] + P.zf[s.JG1 + 1]) - 1.2e-4;
  const sg = 1.0e-4;
  const a = s.n.O2p;
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) {
      const dz = P.zc[j] - zc0;
      a[i * s.nz + j] = 1e18 * Math.exp(-(dz * dz) / (2 * sg * sg));
    }
  }
  const m0 = mass(s, 'O2p'), c0 = centerZ(s, 'O2p');
  const mu = s.spec[1].mu, v = mu * Ez;
  let T = 0;
  for (let k = 0; k < 400 && T < 4e-7; k++) T += s.step();
  const m1 = mass(s, 'O2p'), c1 = centerZ(s, 'O2p');
  const dm = Math.abs(m1 - m0) / m0;
  const dcNum = c1 - c0, dcEx = v * T;
  const errC = Math.abs(dcNum - dcEx) / Math.abs(dcEx);
  ok(dm < 1e-10, `масса не сохранена: ${fmt(dm)}`);
  ok(errC < 0.02, `центр: числ. ${fmt(dcNum)} против v*t = ${fmt(dcEx)}, ошибка ${(errC * 100).toFixed(2)}%`);
  return `dm/m = ${fmt(dm)}, смещение ${fmt(dcNum)} против v*t = ${fmt(dcEx)} (${(errC * 100).toFixed(2)}%), v = ${v.toFixed(1)} м/с, шагов ${s.stepIndex}`;
});

// ─────────────────────────────────────────────────────────── S2: чистая диффузия по r

test('S2  чистая диффузия в r: масса сохранена, поток через ось и r=R = ТОЧНЫЙ ноль', () => {
  // Диффундируют ЭЛЕКТРОНЫ: D_e ~ 4e-2 м²/с против 6e-6 у O2+, иначе за разумное
  // число шагов радиальное расплывание попросту не успевает развиться (замерено:
  // sqrt(2*D_O2p*t) = 5.5 мкм за 200 шагов — тест был бы пустым).
  const s = new DBD2D({
    nr: 24, nz: 24, nzDiel: 3, transportOnly: true, frozenEz: 0,
    betaR: 0, betaG: 0, qB: 1, neutralTransport: false,
    seedBackground: 0, seedSpotAmp: 0, dtMax: 1e-7,
  });
  clearAll(s);
  const P = s.poisson, a = s.n.e;
  const sg = 6e-5;
  for (let i = 0; i < s.nr; i++) {
    const w = 1e18 * Math.exp(-(P.rc[i] * P.rc[i]) / (2 * sg * sg));
    for (let j = s.JG0; j <= s.JG1; j++) a[i * s.nz + j] = w;
  }
  // <r^2> для свободной 2D-диффузии растёт строго как 4*D*t
  const r2 = () => {
    let m = 0, mr = 0;
    for (let i = 0; i < s.nr; i++) for (let j = s.JG0; j <= s.JG1; j++) {
      const w = a[i * s.nz + j] * P.Acell[i] * P.dz[j];
      m += w; mr += w * P.rc[i] * P.rc[i];
    }
    return mr / m;
  };
  const m0 = mass(s, 'e'), R20 = r2();
  let maxAxis = 0, maxWall = 0, T = 0;
  for (let k = 0; k < 200; k++) {
    T += s.step();
    const Gr = s.Gr[0];
    for (let j = 0; j < s.nz; j++) {
      maxAxis = Math.max(maxAxis, Math.abs(Gr[0 * s.nz + j]));
      maxWall = Math.max(maxWall, Math.abs(Gr[s.nr * s.nz + j]));
    }
  }
  const m1 = mass(s, 'e'), R21 = r2();
  const D = s.DerF[1 * s.nz + s.JG0];      // читать ПОСЛЕ шагов: до них грани не заполнены
  const dm = Math.abs(m1 - m0) / m0;
  const growth = (R21 - R20) / T, expect = 4 * D;
  const errG = Math.abs(growth - expect) / expect;
  ok(maxAxis === 0, `поток через ось не ноль: ${maxAxis}`);
  ok(maxWall === 0, `поток через r=R не ноль: ${maxWall}`);
  ok(dm < 1e-10, `масса не сохранена: ${fmt(dm)}`);
  ok(errG < 0.05, `d<r^2>/dt = ${fmt(growth)} против 4D = ${fmt(expect)} (${(errG * 100).toFixed(1)} %)`);
  return `dm/m = ${fmt(dm)}, Gr[ось] = ${maxAxis}, Gr[R] = ${maxWall}, `
       + `d<r²>/dt = ${fmt(growth)} против 4D = ${fmt(expect)} (${(errG * 100).toFixed(2)} %), t = ${fmt(T)} с`;
});

// ─────────────────────────────────────────────────────────── S3: сохранение заряда (БЛОКЕР)

test('S3  СОХРАНЕНИЕ ЗАРЯДА за шаг (ERRATA A5): dQ_vol + dQ_surf = 0, отн. < 1e-6', () => {
  const s = new DBD2D({
    nr: 12, nz: 48, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18, seedSpotSigmaUM: 40,
    photoModule: true, nPhotoSubstep: 1, dtMax: 1e-9,
  });
  s.t = s.period / 4;                      // сразу в область большого напряжения
  let worst = 0, worstStep = -1;
  let q = s.charge(), qc = s.qClip;
  for (let k = 0; k < 60; k++) {
    const dt = s.step();
    const q2 = s.charge(), qc2 = s.qClip;
    const dV = q2.volume - q.volume;
    const dS = q2.surface - q.surface;
    const dC = qc2 - qc;                   // заряд, добавленный клипом floor (ERRATA A2)
    const resid = dV + dS - dC;
    const scale = Math.max(Math.abs(dV), Math.abs(dS), 1e-300);
    const rel = Math.abs(resid) / scale;
    if (rel > worst) { worst = rel; worstStep = k; }
    q = q2; qc = qc2;
    void dt;
  }
  ok(worst < 1e-6, `максимальная относительная невязка ${fmt(worst)} на шаге ${worstStep}`);
  return `max отн. невязка = ${fmt(worst)} (шаг ${worstStep}), клипов floor: ${s.clipCount}, |Q_clip| = ${fmt(Math.abs(s.qClip))} Кл`;
});

test('S3b химия отдельно: extent-схема сохраняет заряд машинно точно', () => {
  const s = new DBD2D({ nr: 4, nz: 24, nzDiel: 3, photoModule: false, seedSpotAmp: 1e19 });
  // насыщаем все сорта, чтобы работали все 34 реакции
  for (const id of SPECIES_IDS) {
    const a = s.n[id];
    for (let i = 0; i < s.nr; i++) for (let j = s.JG0; j <= s.JG1; j++) {
      a[i * s.nz + j] = 1e18 * (1 + 0.3 * Math.sin(i + j));
    }
  }
  for (let k = 0; k < s.ncell; k++) s.EN[k] = 150;   // сильное поле: ионизация активна
  const q0 = s.charge().volume;
  s.p.bgIonizRate = 0;
  s._chemApply(1e-11);
  const q1 = s.charge().volume;
  const rel = Math.abs(q1 - q0) / Math.abs(q0);
  ok(rel < 1e-13, `химия изменила заряд на ${fmt(rel)} относительных`);
  return `dQ/Q после одного применения всех 34 реакций = ${fmt(rel)}`;
});

// ─────────────────────────────────────────────────────────── S4: длинный прогон

test('S4  20000 шагов реального прогона: нет NaN, нет отрицательных плотностей', () => {
  const s = new DBD2D({
    nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: false, dtMax: 2e-8,
  });
  const NSTEP = Number(process.env.S4_STEPS || 20000);
  let minN = Infinity, maxN = 0, nonFinite = 0;
  const st = Date.now();
  for (let k = 0; k < NSTEP; k++) {
    s.step();
    if ((k & 255) === 0 || k === NSTEP - 1) {
      for (const id of SPECIES_IDS) {
        const a = s.n[id];
        for (let i = 0; i < s.nr; i++) for (let j = s.JG0; j <= s.JG1; j++) {
          const v = a[i * s.nz + j];
          if (!Number.isFinite(v)) nonFinite++;
          if (v < minN) minN = v;
          if (v > maxN) maxN = v;
        }
      }
      if (!Number.isFinite(s.Icond) || !Number.isFinite(s.Q) || !Number.isFinite(s.Ugap)) nonFinite++;
    }
  }
  const ms = Date.now() - st;
  ok(nonFinite === 0, `нефинитных значений: ${nonFinite}`);
  ok(minN >= 0, `отрицательная плотность: ${fmt(minN)}`);
  return `${NSTEP} шагов за ${(ms / 1e3).toFixed(1)} с (${(ms / NSTEP).toFixed(2)} мс/шаг), t = ${fmt(s.t)} с, `
       + `min n = ${fmt(minN)}, max n = ${fmt(maxN)}, отбраковок ${s.rejects}, клипов ${s.clipCount}`;
});

// ─────────────────────────────────────────────────────────── S5: U0 = 0

test('S5  U0 = 0: разряда нет, ток проводимости < 1e-12 А', () => {
  const s = new DBD2D({
    nr: 8, nz: 40, nzDiel: 4, U0kV: 0, seedSpotAmp: 0, seedBackground: 1e13,
    photoModule: false, dtMax: 1e-8,
  });
  let maxI = 0, maxE = 0, maxSig = 0;
  for (let k = 0; k < 100; k++) {
    s.step();
    maxI = Math.max(maxI, Math.abs(s.Icond));
    maxE = Math.max(maxE, s.maxEN);
    for (let i = 0; i < s.nr; i++) maxSig = Math.max(maxSig, Math.abs(s.sigLo[i]), Math.abs(s.sigHi[i]));
  }
  ok(maxI < 1e-12, `|I_cond| = ${fmt(maxI)} А`);
  ok(Math.abs(s.Idisp) < 1e-30, `I_disp = ${fmt(s.Idisp)} при U0 = 0`);
  // Поле в точности нулевым НЕ БУДЕТ и не должно: тепловой поток электронов на
  // стенки заряжает поверхности (амбиполярный эффект), и это физика, а не ошибка.
  // Порог 1e-3 Тд = 2.4e-2 В/м — на 7 порядков ниже пробойного поля.
  ok(maxE < 1e-3, `появилось значимое поле: max E/N = ${fmt(maxE)} Тд`);
  return `max|I_cond| = ${fmt(maxI)} А, I_disp = ${fmt(Math.abs(s.Idisp))} А, `
       + `max E/N = ${fmt(maxE)} Тд (амбиполярная зарядка стенок), max|sigma| = ${fmt(maxSig)} Кл/м²`;
});

// ─────────────────────────────────────────────────────────── S6: знак памяти (БЛОКЕР)

test('S6  ЗНАК ПАМЯТИ (ERRATA A1): накопленный sigma ГАСИТ поле в зазоре', () => {
  const s = new DBD2D({
    nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 0, seedBackground: 1e13,
    photoModule: false, dtMax: 1e-9, bgIonizRate: 0,
  });
  s.t = s.period / 4;                       // U = U0 = +10 кВ (максимум)
  // плазменный слой: имитируем состояние ПОСЛЕ пробоя, чтобы стенки реально зарядились
  for (let i = 0; i < s.nr; i++) {
    for (let j = s.JG0; j <= s.JG1; j++) {
      s.n.e[i * s.nz + j] = 1e18; s.n.O2p[i * s.nz + j] = 1e18;
    }
  }
  const U = s._U(s.t);
  const UgapBefore = s.gapVoltage(U, true);
  ok(Math.abs(UgapBefore - s.gapVoltage(U, false)) < 1e-9, 'на старте sigma должен быть нулевым');
  for (let k = 0; k < 300; k++) s.step();
  let sLo = 0, sHi = 0;
  for (let i = 0; i < s.nr; i++) { sLo += s.sigLo[i]; sHi += s.sigHi[i]; }
  sLo /= s.nr; sHi /= s.nr;
  // ключевая проверка: при ТОМ ЖЕ приложенном напряжении накопленный sigma
  // обязан УМЕНЬШИТЬ |U_gap|. Если увеличивает — знак sigma неверен (блокер A1).
  const withSig = s.gapVoltage(U, true);
  const noSig = s.gapVoltage(U, false);
  ok(Math.abs(sLo) > 0, 'sigma не накопился — тест бессодержателен');
  // Геометрия: активный электрод при z = 0, U > 0 => phi падает от 10 кВ к 0 =>
  // Ez = -dphi/dz > 0 => ЭЛЕКТРОНЫ дрейфуют в -z, то есть на НИЖНЮЮ поверхность
  // (она примыкает к аноду). Значит физически ожидается sigLo < 0.
  ok(sLo < 0, `ожидался sigLo < 0 (электроны на анодную поверхность), получено ${fmt(sLo)}`);
  ok(Math.abs(withSig) < Math.abs(noSig),
    `sigma УСИЛИВАЕТ поле: |U_gap| ${fmt(Math.abs(noSig))} -> ${fmt(Math.abs(withSig))} — знак неверен (ERRATA A1)`);
  const drop = (1 - Math.abs(withSig) / Math.abs(noSig)) * 100;
  return `sigLo = ${fmt(sLo)}, sigHi = ${fmt(sHi)} Кл/м²; U_gap при U = ${(U / 1e3).toFixed(1)} кВ: `
       + `${noSig.toFixed(1)} -> ${withSig.toFixed(1)} В (гашение ${drop.toFixed(2)} %)`;
});

// ─────────────────────────────────────────────────────────── S7: 2D -> 1D

test('S7  ВАЛИДАЦИЯ 2D->1D: радиально однородная затравка остаётся однородной (< 1e-8)', () => {
  const opts = {
    nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 0, seedBackground: 1e14,
    photoModule: false, dtMax: 1e-9,
  };
  const s2 = new DBD2D({ ...opts, nr: 8 });
  const s1 = new DBD2D({ ...opts, nr: 1 });
  s2.t = s1.t = s2.period / 4;
  const NS = 150;
  for (let k = 0; k < NS; k++) { s2.step(); s1.step(); }

  // (а) радиальная однородность
  let worst = 0, worstId = '';
  for (const id of [...SPECIES_IDS, 'EN']) {
    const a = id === 'EN' ? s2.EN : s2.n[id];
    for (let j = s2.JG0; j <= s2.JG1; j++) {
      let mn = Infinity, mx = -Infinity, amax = 0;
      for (let i = 0; i < s2.nr; i++) {
        const v = a[i * s2.nz + j];
        mn = Math.min(mn, v); mx = Math.max(mx, v); amax = Math.max(amax, Math.abs(v));
      }
      if (amax > 0) {
        const rel = (mx - mn) / amax;
        if (rel > worst) { worst = rel; worstId = `${id}@j=${j}`; }
      }
    }
  }
  ok(worst < 1e-8, `радиальная неоднородность ${fmt(worst)} у ${worstId}`);

  // (б) совпадение с 1D (nr=1) — таунсендовский режим
  let worst1D = 0, worst1Did = '';
  for (const id of SPECIES_IDS) {
    const a2 = s2.n[id], a1 = s1.n[id];
    for (let j = s2.JG0; j <= s2.JG1; j++) {
      const v2 = a2[0 * s2.nz + j], v1 = a1[j];
      const sc = Math.max(Math.abs(v1), Math.abs(v2));
      if (sc > 1e-30) {
        const rel = Math.abs(v2 - v1) / sc;
        if (rel > worst1D) { worst1D = rel; worst1Did = `${id}@j=${j}`; }
      }
    }
  }
  const dt2 = s2.t, dt1 = s1.t;
  ok(Math.abs(dt2 - dt1) / dt2 < 1e-12, `траектории dt разошлись: t = ${fmt(dt2)} против ${fmt(dt1)}`);
  ok(worst1D < 1e-8, `расхождение с 1D (nr=1) ${fmt(worst1D)} у ${worst1Did}`);
  return `неоднородность по r = ${fmt(worst)} (${worstId}); расхождение 2D(nr=8) vs 1D(nr=1) = ${fmt(worst1D)} (${worst1Did}); `
       + `t = ${fmt(dt2)} с, E/N = ${s2.maxEN.toFixed(2)} Тд`;
});

test('S7b фотомодуль НАРУШАЕТ радиальную однородность (краевой эффект view factor) — замер', () => {
  // Честный замер, а не ассерт «всё хорошо»: прямой view factor в КОНЕЧНОЙ области
  // даёт пристеночной ячейке меньше телесного угла, чем приосевой, поэтому при
  // радиально ОДНОРОДНОМ источнике поток фотоэмиссии радиально неоднороден.
  // Поэтому S7 гоняется с photoModule: false, и это ограничение теста, а не солвера.
  const s = new DBD2D({
    nr: 12, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 0, seedBackground: 1e14,
    photoModule: true, nPhotoSubstep: 1, dtMax: 1e-9,
  });
  s.t = s.period / 4;
  for (let k = 0; k < 5; k++) s.step();
  const F = s.photoEmitFluxL;
  let mn = Infinity, mx = -Infinity;
  for (let i = 0; i < s.nr; i++) { mn = Math.min(mn, F[i]); mx = Math.max(mx, F[i]); }
  const rel = mx > 0 ? (mx - mn) / mx : 0;
  ok(Number.isFinite(rel), 'нефинитный поток фотоэмиссии');
  return `неоднородность Gamma_pe по r при однородном источнике = ${(rel * 100).toFixed(1)} % `
       + `(поток ${fmt(mn)}..${fmt(mx)} м^-2с^-1) — ожидаемый краевой эффект, см. комментарий`;
});

// ─────────────────────────────────────────────────────────── S8: фотофлаги выключены

test('S8  photoIonization/Emission/Detachment = false: бит-в-бит как без фотомодуля', () => {
  const base = {
    nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18, dtMax: 1e-9,
  };
  const a = new DBD2D({ ...base, photoModule: true, photoIonization: false, photoEmission: false, photoDetachment: false });
  const b = new DBD2D({ ...base, photoModule: false });
  ok(a.photo !== null && b.photo === null, 'конфигурация теста неверна');
  ok(a.gammaI === b.gammaI, `gamma_i расходится: ${a.gammaI} против ${b.gammaI}`);
  const NS = 40;
  for (let k = 0; k < NS; k++) { a.step(); b.step(); }
  let diff = 0, diffId = '';
  for (const id of SPECIES_IDS) {
    const x = a.n[id], y = b.n[id];
    for (let k = 0; k < x.length; k++) {
      if (x[k] !== y[k]) { diff++; if (!diffId) diffId = `${id}[${k}]: ${x[k]} != ${y[k]}`; }
    }
  }
  for (let i = 0; i < a.nr; i++) {
    if (a.sigLo[i] !== b.sigLo[i] || a.sigHi[i] !== b.sigHi[i]) { diff++; if (!diffId) diffId = `sigma[${i}]`; }
  }
  ok(a.t === b.t, `время разошлось: ${a.t} != ${b.t}`);
  ok(a.Q === b.Q && a.Icond === b.Icond, `ток/заряд разошлись: ${a.Icond} != ${b.Icond}`);
  ok(diff === 0, `${diff} различий, первое: ${diffId}`);
  return `${NS} шагов, ${SPECIES_IDS.length} сортов x ${a.ncell} ячеек: 0 различий; t = ${fmt(a.t)}, Q = ${fmt(a.Q)} Кл`;
});

// ─────────────────────────────────────────────────────────── доп.: ток и весовое поле

test('S9  ток Сато–Морроу: E_L = 1/d_eff (ERRATA A3), I_disp = C_cell*dU/dt', () => {
  const s = new DBD2D({ nr: 8, nz: 40, nzDiel: 4, photoModule: false, seedSpotAmp: 0 });
  const dEff = 1e-3 + 1e-3 / 9;
  ok(Math.abs(s.dEff - dEff) / dEff < 1e-12, `d_eff = ${s.dEff}`);
  ok(Math.abs(s.ELz - 900) / 900 < 1e-12, `E_L = ${s.ELz} 1/м`);
  // C_cell = eps0*pi*R^2/d_eff = 6.259 фФ (NUMERICS_2D §7.2)
  const Cref = 8.8541878128e-12 * Math.PI * 0.5e-3 * 0.5e-3 / dEff;
  ok(Math.abs(s.Ccell - Cref) / Cref < 1e-12, `C_cell = ${fmt(s.Ccell)}`);
  // весовой потенциал, посчитанный численно в poisson2d, не должен зависеть от r
  const P = s.poisson;
  let dev = 0;
  for (let j = 0; j < s.nz; j++) {
    let mn = Infinity, mx = -Infinity;
    for (let i = 0; i < s.nr; i++) { const v = P.phiL[i * s.nz + j]; mn = Math.min(mn, v); mx = Math.max(mx, v); }
    dev = Math.max(dev, mx - mn);
  }
  ok(dev < 1e-12, `phi_L зависит от r на ${fmt(dev)}`);
  return `d_eff = ${fmt(s.dEff)} м, E_L = ${s.ELz.toFixed(1)} 1/м, C_cell = ${(s.Ccell * 1e15).toFixed(3)} фФ, max|dphi_L(r)| = ${fmt(dev)}`;
});

test('S10 сохранение заряда за 500 шагов подряд (интегрально, ERRATA C V12 < 1e-4)', () => {
  const s = new DBD2D({
    nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18, photoModule: false, dtMax: 1e-9,
  });
  s.t = s.period / 4;
  const q0 = s.charge(), qc0 = s.qClip;
  let scale = 0;
  for (let k = 0; k < 500; k++) {
    s.step();
    scale = Math.max(scale, Math.abs(s.charge().surface - q0.surface));
  }
  const q1 = s.charge();
  const resid = (q1.volume - q0.volume) + (q1.surface - q0.surface) - (s.qClip - qc0);
  const rel = Math.abs(resid) / Math.max(scale, 1e-300);
  ok(rel < 1e-6, `интегральная невязка заряда ${fmt(rel)}`);
  return `невязка за 500 шагов = ${fmt(rel)} от перенесённого заряда ${fmt(scale)} Кл; клипов ${s.clipCount}`;
});

test('S11 контракт state <-> recorder.mjs: кадры и ряды пишутся без потерь полей', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'dbd2d-'));
  const s = new DBD2D({ nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18, photoModule: false });
  const rec = new Recorder({ dir, runId: 'test', level: 'compact', params: { U0kV: 10 }, tEnd: 1e-8 });
  for (let k = 0; k < 25; k++) {
    s.step();
    const st = s.state;
    if (rec.shouldRecord(st)) rec.record(st);
    else rec.record(st);
  }
  const mpath = rec.close({ note: 'solver contract test' });
  const man = JSON.parse(fs.readFileSync(mpath, 'utf8'));
  const missing = [];
  for (const f of man.fields || []) if (f.present === false) missing.push(f.name);
  ok(man.stats.frames > 0, 'ни одного кадра не записано');
  ok(man.stats.seriesRecords > 0, 'ни одной записи ряда');
  ok(missing.length === 0, `поля state отсутствуют для recorder: ${missing.join(', ')}`);
  ok(rec.nonFiniteHits === 0, `recorder увидел ${rec.nonFiniteHits} нефинитных значений`);
  const sz = fs.statSync(path.join(dir, 'frames.bin')).size;
  ok(sz === 64 + man.stats.frames * man.frameLayout.frameStride, `размер frames.bin = ${sz}`);
  fs.rmSync(dir, { recursive: true, force: true });
  return `кадров ${man.stats.frames}, рядов ${man.stats.seriesRecords}, полей ${(man.fields || []).length} (все present), `
       + `frames.bin ${sz} Б сходится с раскладкой, нефинитных 0`;
});

// ═══════════════════════════════════════════════════════════
// W1..W6 — ПРИСТЕНОЧНОЕ ГУ ХАГЕЛААРА (ERRATA «STOP №2», DIVERGENCE_ANALYSIS §8.4/§9)
// ═══════════════════════════════════════════════════════════

/** Поверхностный интеграл потоков на стенку side (0 — нижняя): [с^-1]. */
function wallFlux(s, side) {
  const P = s.poisson;
  let ge = 0, gp = 0, gn = 0;
  for (let i = 0; i < s.nr; i++) {
    const A = P.Acell[i];
    ge += s.wallGe[side * s.nr + i] * A;
    gp += s.wallGpos[side * s.nr + i] * A;
    gn += s.wallGneg[side * s.nr + i] * A;
  }
  return { ge, gp, gn };
}
/** max |E_z| на гранях газ/диэлектрик. */
function maxEwall(s) {
  let m = 0;
  for (let i = 0; i < s.nr; i++) {
    const tz = i * s.nzf;
    m = Math.max(m, Math.abs(s.EzF[tz + s.JG0]), Math.abs(s.EzF[tz + s.JG1 + 1]));
  }
  return m;
}

test('W1  wallBC=legacy воспроизводит ПРЕЖНЕЕ поведение бит-в-бит (эталон снят до правки)', () => {
  // Эталон снят с версии solver2d.mjs ДО переноса ГУ Хагелаара (тогдашний дефолт
  // wallBC:'thermal'), той же конфигурацией и тем же числом шагов. Числа записаны
  // с 17 значащими цифрами, то есть round-trip точно восстанавливает double.
  const GOLD = {
    e: 1.02437784510656465e+7, O2p: 1.03780531784191951e+7, O4p: 2.89729763261792134e+5,
    Om: 2.23234842936869914e+5, O2m: 2.97965666458867236e+3, O3m: 2.33247315394539692e+3,
    O: 1.11374825503955483e+8, O3: 6.37302327985253783e+2,
    O2a: 4.65177294011503551e+6,
    t: 2.50000507896301119e-5, Q: 1.05295507875422462e-14, Icond: 5.08321003789723148e-4,
    sigLo0: -2.67679164576923044e-6, sigHi0: -8.98776892885798488e-12,
    sigLoLast: -7.64045156525141024e-11,
  };
  const s = new DBD2D({ nr: 6, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: false, dtMax: 1e-9, wallBC: 'legacy' });
  s.t = s.period / 4;
  for (let k = 0; k < 60; k++) s.step();
  const got = { t: s.t, Q: s.Q, Icond: s.Icond,
    sigLo0: s.sigLo[0], sigHi0: s.sigHi[0], sigLoLast: s.sigLo[s.nr - 1] };
  for (const id of SPECIES_IDS) got[id] = mass(s, id);
  const bad = [];
  for (const key of Object.keys(GOLD)) if (got[key] !== GOLD[key]) bad.push(`${key}: ${got[key]} != ${GOLD[key]}`);
  ok(bad.length === 0, `${bad.length} расхождений с эталоном:\n      ${bad.join('\n      ')}`);
  // 'thermal' — устаревшее имя того же режима, оно обязано давать ровно то же
  const s2 = new DBD2D({ nr: 6, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: false, dtMax: 1e-9, wallBC: 'thermal' });
  s2.t = s2.period / 4;
  for (let k = 0; k < 60; k++) s2.step();
  ok(s2.t === s.t && s2.Q === s.Q && mass(s2, 'e') === got.e, 'alias thermal != legacy');
  return `${Object.keys(GOLD).length} величин совпали бит-в-бит (n_e = ${fmt(got.e)}, Q = ${fmt(got.Q)} Кл), alias 'thermal' тождественен`;
});

test('W2  на АНОДЕ ионы отталкиваются: kw_ion = 0 ТОЧНО, на катоде kw_ion > 0', () => {
  // E_z = +5e6 В/м: положительные ионы дрейфуют в +z, то есть ОТ нижней стенки
  // (она анод для ионов) и К верхней. Электроны — наоборот.
  const s = new DBD2D({ nr: 3, nz: 32, nzDiel: 4, frozenEz: 1e7, photoModule: false,
    seedSpotAmp: 1e16, dtMax: 1e-12 });
  s.step();
  const pos = ['O2p', 'O4p'];
  let supp = 0;
  for (let i = 0; i < s.nr; i++) {
    for (const id of pos) {
      const a = s.wallCoeff(0, id, i);         // анод для положительных ионов
      const c = s.wallCoeff(1, id, i);         // катод
      ok(a.kw === 0, `kw(${id}) на аноде = ${a.kw} (обязан быть ТОЧНЫМ нулём)`);
      ok(a.lin === 0, `lin(${id}) на аноде = ${a.lin}: запертый поток не имеет отклика по полю`);
      ok(c.kw > 0, `kw(${id}) на катоде = ${c.kw} (обязан быть > 0)`);
    }
    // Электроны — зеркально, но у них ¼v_th и mu_e|E| ОДНОГО порядка, поэтому
    // точный ноль не обязан достигаться: поток лишь подавляется в десятки раз.
    // Именно эта почти-компенсация и есть плавающий потенциал (§9.2: kw_e = 0.8 %
    // от каждого слагаемого). Требуем подавления, а не нуля.
    const eRep = s.wallCoeff(1, 'e', i).kw, eCol = s.wallCoeff(0, 'e', i).kw;
    ok(eCol > 0, 'kw(e) на собирающей стенке должен быть > 0');
    ok(eRep >= 0 && eRep < 0.05 * eCol, `kw(e) на отталкивающей стенке ${eRep} не подавлен (собирающая ${eCol})`);
    supp = Math.max(supp, eRep / eCol);
  }
  // а старая форма гнала на анод полный тепловой поток ¼v_th
  const sl = new DBD2D({ nr: 3, nz: 32, nzDiel: 4, frozenEz: 5e6, photoModule: false,
    seedSpotAmp: 1e16, dtMax: 1e-12, wallBC: 'legacy' });
  sl.step();
  const kwLeg = sl.wallCoeff(0, 'O2p', 0).kw;
  ok(kwLeg > 0, 'контроль теста: legacy обязан давать НЕнулевой поток на анод');
  return `hagelaar: kw(O2+|анод) = 0 точно; legacy: ${fmt(kwLeg)} м/с (= ¼v_th — вторая ошибка старого ГУ), `
       + `kw(O2+|катод) = ${fmt(s.wallCoeff(1, 'O2p', 0).kw)} м/с; `
       + `электроны на отталкивающей стенке подавлены до ${(supp * 100).toFixed(2)} % от собирающей`;
});

test('W4  поле у стенки НЕ разгоняется: max E_wall ограничен (в 1D садится на 1.4e7 В/м)', () => {
  // Быстрый вариант: старт на пике напряжения (U ≈ 10 кВ >> U_br), зажигание за
  // доли нс, дальше фиксируется поле на пристеночной грани.
  const mk = (bc) => new DBD2D({ nr: 4, nz: 48, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: false, dtMax: 1e-9, wallBC: bc });
  const s = mk('hagelaar');
  s.t = s.period / 4;
  let mx = 0, ne = 0, en = 0;
  for (let k = 0; k < 4000; k++) {
    s.step();
    mx = Math.max(mx, maxEwall(s)); ne = Math.max(ne, s.neMax); en = Math.max(en, s.maxEN);
  }
  ok(Number.isFinite(mx) && mx < 5e7, `max E_wall = ${fmt(mx)} В/м (петля разгона §8.4 уводит его на 1.3e8)`);
  ok(en < 2e3, `max E/N = ${fmt(en)} Тд: поле всё-таки разогналось`);
  return `max E_wall = ${fmt(mx)} В/м, max E/N = ${en.toFixed(0)} Тд, max n_e = ${fmt(ne)} м^-3 за ${s.stepIndex} шагов`;
});

test('W5  зарядовый баланс d/dt[∫rho dV + ∫sigma dS] = 0, невязка < 1e-8', () => {
  const s = new DBD2D({ nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: true, nPhotoSubstep: 1, dtMax: 1e-9 });
  s.t = s.period / 4;
  const q0 = s.charge(), qc0 = s.qClip;
  let worst = 0, scale = 0;
  let q = q0, qc = qc0;
  for (let k = 0; k < 500; k++) {
    s.step();
    const q2 = s.charge(), qc2 = s.qClip;
    const dV = q2.volume - q.volume, dS = q2.surface - q.surface, dC = qc2 - qc;
    const rel = Math.abs(dV + dS - dC) / Math.max(Math.abs(dV), Math.abs(dS), 1e-300);
    if (rel > worst) worst = rel;
    q = q2; qc = qc2;
    scale = Math.max(scale, Math.abs(q2.surface - q0.surface));
  }
  const intg = Math.abs((q.volume - q0.volume) + (q.surface - q0.surface) - (s.qClip - qc0)) / Math.max(scale, 1e-300);
  ok(worst < 1e-8, `пошаговая невязка ${fmt(worst)}`);
  ok(intg < 1e-8, `интегральная невязка за 500 шагов ${fmt(intg)}`);
  return `пошаговая ${fmt(worst)}, интегральная ${fmt(intg)} от перенесённого заряда ${fmt(scale)} Кл (500 шагов)`;
});

test('W6  S7 (радиальная однородность) и S8 (бит-в-бит без фото) держатся в ОБОИХ режимах ГУ', () => {
  const out = [];
  for (const bc of ['hagelaar', 'legacy']) {
    // (a) 2D -> 1D: kw = kw(r), при радиально однородном поле обязан остаться однородным
    const opts = { nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 0, seedBackground: 1e14,
      photoModule: false, dtMax: 1e-9, wallBC: bc };
    const s2 = new DBD2D({ ...opts, nr: 8 });
    const s1 = new DBD2D({ ...opts, nr: 1 });
    s2.t = s1.t = s2.period / 4;
    for (let k = 0; k < 150; k++) { s2.step(); s1.step(); }
    let worst = 0;
    for (const id of SPECIES_IDS) {
      for (let j = s2.JG0; j <= s2.JG1; j++) {
        let mn = Infinity, mx = -Infinity, am = 0;
        for (let i = 0; i < s2.nr; i++) {
          const v = s2.n[id][i * s2.nz + j];
          mn = Math.min(mn, v); mx = Math.max(mx, v); am = Math.max(am, Math.abs(v));
        }
        if (am > 0) worst = Math.max(worst, (mx - mn) / am);
      }
    }
    ok(worst < 1e-8, `[${bc}] радиальная неоднородность ${fmt(worst)}`);
    // и совпадение kw по радиусу — прямая проверка, что kw(r) не «поехал»
    let kwSpread = 0;
    for (let side = 0; side < 2; side++) for (const id of ['e', 'O2p']) {
      let mn = Infinity, mx = -Infinity;
      for (let i = 0; i < s2.nr; i++) { const v = s2.wallCoeff(side, id, i).kw; mn = Math.min(mn, v); mx = Math.max(mx, v); }
      if (mx > 0) kwSpread = Math.max(kwSpread, (mx - mn) / mx);
    }
    ok(kwSpread < 1e-8, `[${bc}] kw(r) неоднороден на ${fmt(kwSpread)}`);
    // (б) выключенные фотофлаги = бит-в-бит как без фотомодуля
    const base = { nr: 8, nz: 40, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18, dtMax: 1e-9, wallBC: bc };
    const a = new DBD2D({ ...base, photoModule: true, photoIonization: false, photoEmission: false, photoDetachment: false });
    const b = new DBD2D({ ...base, photoModule: false });
    for (let k = 0; k < 40; k++) { a.step(); b.step(); }
    let diff = 0;
    for (const id of SPECIES_IDS) { const x = a.n[id], y = b.n[id]; for (let k = 0; k < x.length; k++) if (x[k] !== y[k]) diff++; }
    ok(diff === 0 && a.t === b.t && a.Q === b.Q, `[${bc}] ${diff} различий при выключенных фотофлагах`);
    out.push(`${bc}: неоднородность ${fmt(worst)}, разброс kw(r) ${fmt(kwSpread)}, различий S8 = 0`);
  }
  return out.join('; ');
});

test('W3  условие плавающей стенки: поле у катодного барьера садится на плавающий потенциал (5 %)', () => {
  // ЧТО ИМЕННО ПРОВЕРЯЕТСЯ И ПОЧЕМУ ИМЕННО ТАК.
  // Буквальная форма «мгновенный Gamma_e = Gamma_i на каждом барьере» в ведомом
  // синусом ДБР невыполнима по построению: барьеры ЗАРЯЖАЮТСЯ (это и есть memory
  // effect), причём в полупериод один барьер анодный (ионы отталкиваются, kw_i = 0
  // ТОЧНО — тест W2), другой катодный (электроны отталкиваются). Ноль суммарного
  // тока на поверхность наступает только в среднем за период.
  // Содержательная форма — та, в которой этот же результат прочитан в 1D
  // (DIVERGENCE_ANALYSIS §9.2): поле у стенки САМО встаёт туда, где электронный и
  // ионный пристеночные коэффициенты сравниваются, то есть где два слагаемых kw_e
  // почти компенсируют друг друга. В 1D: ¼v_th,e = 4.959e5, mu_e|E| = 4.921e5,
  // остаток kw_e = 3.733e3 м/с — 0.8 % от каждого слагаемого, и он совпал с
  // kw_O2+ = 3.500e3 м/с. Здесь проверяется ровно это: компенсация с точностью 5 %.
  const s = new DBD2D({ nr: 4, nz: 48, nzDiel: 4, U0kV: 10, seedSpotAmp: 1e18,
    photoModule: false, dtMax: 1e-9 });
  s.t = s.period / 4;
  for (let k = 0; k < 4000; k++) s.step();
  // катодный барьер — тот, куда идут ИОНЫ
  const w0 = wallFlux(s, 0), w1 = wallFlux(s, 1);
  const side = w1.gp > w0.gp ? 1 : 0;
  const w = side === 1 ? w1 : w0;
  const jf = side === 0 ? s.JG0 : s.JG1 + 1;
  const Ew = s.EzF[jf], en = s.ENzF[jf], mu = s.muEzF[jf];
  const r = s.reflE;
  const th = 0.25 * (1 - r) / (1 + r) * Math.sqrt((8 * QE * (2 / 3) * meanEnergy(en)) / (Math.PI * ME));
  const dr = (mu * Math.abs(Ew)) / (1 + r);
  const mism = Math.abs(th - dr) / th;
  ok(w.gp > 0, `на барьере ${side} нет ионного потока (${fmt(w.gp)}): разряд не зажёгся`);
  ok(mism < 0.05, `слагаемые kw_e не скомпенсированы: ¼v_th' = ${fmt(th)}, mu_e|E|' = ${fmt(dr)}, `
    + `рассогласование ${(mism * 100).toFixed(2)} % (в 1D 0.8 %)`);
  ok(s.wallCoeff(side, 'e', 0).kw <= 0.05 * th, 'электронный поток на катод не подавлен');
  return `катодный барьер ${side}: E_w = ${fmt(Ew)} В/м (${en.toFixed(0)} Тд), ¼v_th' = ${fmt(th)} против `
       + `mu_e|E|' = ${fmt(dr)} м/с — компенсация ${(mism * 100).toFixed(2)} % (1D: 0.8 %); `
       + `kw_e = ${fmt(s.wallCoeff(side, 'e', 0).kw)}, kw_O2+ = ${fmt(s.wallCoeff(side, 'O2p', 0).kw)} м/с; `
       + `Gamma_e/Gamma_i = ${(w.ge / w.gp).toFixed(4)} (отрицательное значение = чистая вторичная эмиссия, -gamma)`;
});

// ───────────────────────────────────────────────────────────

console.log(`\n${passed}/${passed + failed} тестов пройдено за ${((Date.now() - t0) / 1e3).toFixed(1)} с`);
if (failed) process.exit(1);
void QE;
