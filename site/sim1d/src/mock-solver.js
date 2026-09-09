// mock-solver.js — stand-in for src/solver.js with the IDENTICAL public API.
// Purpose: let the UI be developed/verified while the real core is written in parallel.
//
// It is NOT the fluid model. It is a lumped-circuit DBD (C_diel in series with the gas
// gap, gas conductance driven by a Townsend-like avalanche ODE) whose scalar output is
// then *painted* onto a 1D grid with analytic shape functions. That is enough to produce
// physically plausible signals: ns current pulses, self-quenching via surface charge,
// a parallelogram Q-V figure, a cathode-directed bright zone, electronegative bulk.
//
// Everything the UI is allowed to touch is in the CONTRACT; nothing else is public API.

const EPS0 = 8.8541878128e-12;
const QE = 1.602176634e-19;
const KB = 1.380649e-23;
const TORR = 133.322368;

const SPECIES_KEYS = ['e', 'O2p', 'O4p', 'Om', 'O2m', 'O', 'O3', 'O2a'];

const DEFAULTS = {
  gapMM: 1.0,
  dielMM1: 0.5,
  dielMM2: 0.5,
  epsR: 9,
  areaCM2: 1.0,
  U0kV: 10.0,
  freqKHz: 10.0,
  pressureTorr: 760,
  tempK: 300,
  gamma: 0.05,
  seedDensity: 1e13,
  nCells: 300,
  ballastOhm: 0,
  mode: 'sine',
};

const clamp = (v, lo, hi) => (v < lo ? lo : v > hi ? hi : v);
const sexp = (a) => Math.exp(clamp(a, -600, 600));

export class DBDSolver {
  constructor(params = {}) {
    this.params = Object.assign({}, DEFAULTS, params);
    this.targetSimTime = Infinity; // UI may set this; advance() respects it
    this._historyCapacity = 65536;
    this.history = {
      t: new Float64Array(this._historyCapacity),
      Uapp: new Float64Array(this._historyCapacity),
      Ugap: new Float64Array(this._historyCapacity),
      current: new Float64Array(this._historyCapacity),
      charge: new Float64Array(this._historyCapacity),
      len: 0,
      head: 0,
      capacity: this._historyCapacity,
    };
    this.periodStats = {
      energyPerPeriodJ: 0,
      powerW: 0,
      Cdiel: 0,
      Ccell: 0,
      UburnkV: 0,
    };
    this._buildGrid();
    this._allocate();
    this.reset();
  }

  // ---------------------------------------------------------------- geometry
  _derive() {
    const p = this.params;
    this.d = p.gapMM * 1e-3;
    this.d1 = p.dielMM1 * 1e-3;
    this.d2 = p.dielMM2 * 1e-3;
    this.A = p.areaCM2 * 1e-4;
    this.N = (p.pressureTorr * TORR) / (KB * p.tempK);
    this.Cd = (EPS0 * p.epsR * this.A) / Math.max(1e-9, this.d1 + this.d2);
    this.Cg = (EPS0 * this.A) / this.d;
    this.Ccell = (this.Cd * this.Cg) / (this.Cd + this.Cg);
    this.omega = 2 * Math.PI * p.freqKHz * 1e3;
    this.period = 1 / (p.freqKHz * 1e3);
  }

  _buildGrid() {
    this._derive();
    const nx = Math.max(32, Math.round(this.params.nCells));
    this.nx = nx;
    this.x = new Float64Array(nx);
    this.xFaces = new Float64Array(nx + 1);
    this.gasMask = new Uint8Array(nx).fill(1);
    const beta = 2.4;
    const th = Math.tanh(beta);
    for (let i = 0; i <= nx; i++) {
      const s = (2 * i) / nx - 1;
      this.xFaces[i] = (this.d * (1 + Math.tanh(beta * s) / th)) / 2;
    }
    for (let i = 0; i < nx; i++) this.x[i] = 0.5 * (this.xFaces[i] + this.xFaces[i + 1]);
  }

  _allocate() {
    const nx = this.nx;
    this.n = {};
    for (const k of SPECIES_KEYS) this.n[k] = new Float64Array(nx);
    this.E = new Float64Array(nx);
    this.phi = new Float64Array(nx);
    this.rho = new Float64Array(nx);
    this.ionizRate = new Float64Array(nx);
    this.EN = new Float64Array(nx);
    this._state = {
      t: 0, dt: 0,
      x: this.x, xFaces: this.xFaces, gasMask: this.gasMask,
      n: this.n, E: this.E, phi: this.phi, rho: this.rho,
      ionizRate: this.ionizRate, EN: this.EN,
      sigmaL: 0, sigmaR: 0, Uapp: 0, Ugap: 0,
      current: 0, condCurrent: 0, dispCurrent: 0, charge: 0, power: 0,
      o3ppm: 0, maxEN: 0, peakCurrent: 0, breakdownsPerPeriod: 0,
    };
  }

  // ------------------------------------------------------------------ public
  setParams(obj) {
    const p = this.params;
    const geomKeys = ['gapMM', 'dielMM1', 'dielMM2', 'nCells'];
    let regrid = false;
    for (const k of Object.keys(obj)) {
      if (obj[k] === undefined) continue;
      if (geomKeys.includes(k) && obj[k] !== p[k]) regrid = true;
      p[k] = obj[k];
    }
    if (regrid) {
      this._buildGrid();
      this._allocate();
      this.reset();
    } else {
      this._derive();
    }
  }

  reset() {
    this._derive();
    this.t = 0;
    this.dt = 1e-11;
    this.ne = Math.max(1e6, this.params.seedDensity); // gap-averaged, m^-3
    this.nOm = 1e12;
    this.nO2m = 1e12;
    this.nO = 1e15;
    this.nO3 = 0;
    this.nO2a = 1e17;
    this.Ugap = 0;
    this.Q = 0;               // charge delivered by external circuit, C
    this.Qcond = 0;           // charge that crossed the gap (surface charge), C
    this.Icond = 0;
    this.Itot = 0;
    this.Idisp = 0;
    this.frontPos = 0.5;
    this.glow = 0;
    this.peakCurrent = 0;
    this.peakCurrentHist = 1e-9;
    this.maxEN = 0;
    this._bdCount = 0;
    this._bdPrev = 0;
    this._bdLast = -1;
    this.breakdownsPerPeriod = 0;
    this.energyAccum = 0;
    this.periodIndex = 0;
    this._periodT0 = 0;
    this._segSxx = 0; this._segSxy = 0; this._segSx = 0; this._segSy = 0; this._segN = 0;
    this._segOn = false;
    this._slopeOnSum = 0; this._slopeOnN = 0;
    this._slopeOffSum = 0; this._slopeOffN = 0;
    this._uburnSum = 0; this._uburnN = 0;
    this.events = [];
    this.history.len = 0;
    this.history.head = 0;
    this._histAccum = 0;
    this._profileDirty = true;
    this.periodStats.energyPerPeriodJ = 0;
    this.periodStats.powerW = 0;
    this.periodStats.Cdiel = this.Cd;
    this.periodStats.Ccell = this.Ccell;
    this.periodStats.UburnkV = 0;
    this._refreshProfiles();
  }

  get state() {
    if (this._profileDirty) this._refreshProfiles();
    return this._state;
  }

  // -------------------------------------------------------------- integrator
  _uapp(t) {
    const U0 = this.params.U0kV * 1e3;
    const ph = this.omega * t;
    switch (this.params.mode) {
      case 'square': return U0 * (Math.sin(ph) >= 0 ? 1 : -1);
      case 'pulsed': {
        const f = ((ph / (2 * Math.PI)) % 1 + 1) % 1;
        if (f < 0.05) return U0 * Math.sin((Math.PI * f) / 0.05);
        if (f > 0.5 && f < 0.55) return -U0 * Math.sin((Math.PI * (f - 0.5)) / 0.05);
        return 0;
      }
      default: return U0 * Math.sin(ph);
    }
  }

  step() {
    const p = this.params;
    const d = this.d;
    const Eg = this.Ugap / d;                    // V/m, signed
    const absE = Math.abs(Eg);
    const EN = (absE / this.N) * 1e21;           // Td

    // Townsend ionisation / attachment (crude but monotone, LFA-flavoured)
    const mue = 0.05 * (2.45e25 / this.N);
    const alpha = absE > 1 ? 4.5e-20 * this.N * sexp(-9.0e-19 * this.N / absE) : 0; // 1/m
    const attach = 1.2e3 * (this.N / 2.45e25);   // 1/m, net attachment
    const vd = mue * absE;
    // Net growth rate. The DECAY branch is deliberately capped at ~1/(120 ns):
    // if the plasma vanished as fast as the local Townsend balance suggests, the
    // lumped model would settle into a steady glow and never pulse. Slow decay
    // makes the conductance overshoot, drives U_gap below the burning voltage and
    // turns the circuit into the relaxation oscillator a filamentary DBD actually is.
    let netNu = (alpha - attach) * vd;           // 1/s
    if (netNu < -8e6) netNu = -8e6;
    const seed = p.seedDensity * 1e3 + 1e18 * this.nO2a * 3e-22; // detachment-ish source
    const recomb = 2e-13 * this.ne * (this.ne + this.nOm + this.nO2m);

    // conduction: G = sigma*A/d  (Icond = G*Ugap)
    const sigmaCond = QE * this.ne * mue * (0.6 + 0.4 * Math.tanh(EN / 120));
    let G = (sigmaCond * this.A) / d;
    if (p.ballastOhm > 0) G = 1 / (1 / G + p.ballastOhm);
    const Ct = this.Cd + this.Cg;

    // adaptive dt: resolve BOTH the avalanche e-folding time and the circuit RC.
    // Missing the RC term is what turns the pulse into a numerical limit cycle.
    const tauAv = netNu !== 0 ? 0.15 / Math.abs(netNu) : 1;
    const tauRC = G > 0 ? Ct / G : 1;
    let dt = Math.min(2e-8, Math.max(1e-13, Math.min(tauAv, 0.5 * tauRC)));
    dt = Math.min(dt, this.period / 2000);
    this.dt = dt;

    // avalanche (semi-implicit in the linear growth term → no overshoot)
    const gr = dt * netNu;
    let ne = gr < 0.5 ? (this.ne + dt * (seed - recomb)) / (1 - gr)
      : this.ne + dt * (this.ne * netNu + seed - recomb);
    this.ne = clamp(ne, 1e4, 1e23);

    // circuit, backward-Euler on the conduction term:
    //   dUgap/dt = (Cd*dUapp/dt - G*Ugap)/Ct
    const tNext = this.t + dt;
    const Ua = this._uapp(this.t);
    const UaNext = this._uapp(tNext);
    const dUa = (UaNext - Ua) / dt;
    const Uprev = this.Ugap;
    this.Ugap = (this.Ugap + (dt * this.Cd * dUa) / Ct) / (1 + (dt * G) / Ct);
    const Icond = G * this.Ugap;
    const dUgap = (this.Ugap - Uprev) / dt;
    this.Icond = Icond;
    this.Idisp = this.Cg * dUgap;
    this.Itot = this.Cd * (dUa - dUgap);
    this.Qcond += Icond * dt;
    this.Q += this.Itot * dt;

    // negative-ion / neutral chemistry (slow track)
    const kAtt = attach * vd;
    this.nOm += dt * (0.55 * kAtt * this.ne - 1.1e-13 * this.nOm * 1e19 * 0);
    this.nO2m += dt * (0.45 * kAtt * this.ne);
    this.nOm = clamp(this.nOm, 0, 1e22);
    this.nO2m = clamp(this.nO2m, 0, 1e22);
    const Sion = this.ne * alpha * vd;
    this.nO += dt * (2.2 * Sion - 1.4e-20 * this.nO * this.N * 0.02);
    this.nO3 += dt * (1.4e-20 * this.nO * this.N * 0.02 - 2e-4 * this.nO3);
    this.nO2a += dt * (0.35 * Sion - 50 * this.nO2a);
    this.nO = clamp(this.nO, 0, 1e24);
    this.nO3 = clamp(this.nO3, 0, 1e24);
    this.nO2a = clamp(this.nO2a, 0, 1e24);

    // ionisation-wave front position (cathode-directed) + glow buffer
    const pol = this.Ugap >= 0 ? 1 : -1;
    const glowNow = clamp(Sion / 1e28, 0, 4);
    if (glowNow > this.glow) {
      this.glow = glowNow;
      this.frontPos = clamp(0.5 + pol * 0.42 * Math.tanh(glowNow * 3), 0.04, 0.96);
    } else {
      this.glow *= sexp(-dt / 2e-7);
    }
    this.Sion = Sion;

    this.t = tNext;

    // diagnostics — the DISCHARGE current is the conduction current, not the total.
    // Classifying on |I_total| would fold the ~5 mA displacement current into the
    // "discharge on" branch and wreck both the Manley slopes and the pulse count.
    const Iabs = Math.abs(this.Icond);
    if (Iabs > this.peakCurrentHist) this.peakCurrentHist = Iabs;
    if (EN > this.maxEN) this.maxEN = EN;
    this.energyAccum += Math.abs(Ua * this.Itot) * dt;

    // breakdown detection: rising crossing of 0.1*I_peak_hist with 5 ns dead time
    const thr = 0.1 * this.peakCurrentHist;
    if (Iabs > thr && this._bdPrev <= thr && this.t - this._bdLast > 5e-9) {
      this._bdLast = this.t;
      this._bdCount++;
      this.events.push({ t: this.t, polarity: pol, I: this.Itot });
      if (this.events.length > 4096) this.events.splice(0, 2048);
    }
    this._bdPrev = Iabs;

    // Manley slopes. The four sides of the parallelogram are TWO pairs of parallel
    // lines with different offsets — a single pooled least-squares fit over a pair
    // returns a slope somewhere between the line slope and the line joining the
    // pair, which is why the naive version reports C_cell > C_diel. Fit each
    // contiguous branch separately, then average the per-branch slopes.
    const on = Iabs > 0.05 * this.peakCurrentHist;
    if (on !== this._segOn) { this._closeSegment(); this._segOn = on; }
    if (on) { this._uburnSum += Math.abs(this.Ugap); this._uburnN++; }
    const ux = Ua, qy = this.Q;
    this._segSxx += ux * ux; this._segSxy += ux * qy;
    this._segSx += ux; this._segSy += qy; this._segN++;

    this._pushHistory(Ua);
    if (this.t - this._periodT0 >= this.period) this._closePeriod();

    this._profileDirty = true;
    return dt;
  }

  _pushHistory(Ua) {
    // decimate raw samples into the ring at ~4 ns granularity to keep the buffer useful
    this._histAccum += this.dt;
    const h = this.history;
    if (this._histAccum < 4e-9 && h.len > 0) return;
    this._histAccum = 0;
    const i = h.head;
    h.t[i] = this.t;
    h.Uapp[i] = Ua;
    h.Ugap[i] = this.Ugap;
    h.current[i] = this.Itot;
    h.charge[i] = this.Q;
    h.head = (i + 1) % h.capacity;
    if (h.len < h.capacity) h.len++;
  }

  _closeSegment() {
    const n = this._segN;
    if (n >= 16) {
      const den = n * this._segSxx - this._segSx * this._segSx;
      const spanU = Math.sqrt(Math.max(0, this._segSxx / n - (this._segSx / n) ** 2));
      // ignore branches that barely move in U — their slope is numerically meaningless
      if (Math.abs(den) > 1e-30 && spanU > 200) {
        const sl = (n * this._segSxy - this._segSx * this._segSy) / den;
        if (sl > 0 && isFinite(sl)) {
          if (this._segOn) { this._slopeOnSum += sl; this._slopeOnN++; }
          else { this._slopeOffSum += sl; this._slopeOffN++; }
        }
      }
    }
    this._segSxx = this._segSxy = this._segSx = this._segSy = this._segN = 0;
  }

  _closePeriod() {
    const T = this.t - this._periodT0;
    this._periodT0 = this.t;
    this.periodIndex++;
    this.periodStats.energyPerPeriodJ = this.energyAccum;
    this.periodStats.powerW = this.energyAccum / Math.max(1e-12, T);
    this._closeSegment();
    if (this._slopeOnN) this.periodStats.Cdiel = this._slopeOnSum / this._slopeOnN;
    if (this._slopeOffN) this.periodStats.Ccell = this._slopeOffSum / this._slopeOffN;
    // burning voltage = mean |U_gap| while the discharge conducts (Manley U_min)
    if (this._uburnN > 0) this.periodStats.UburnkV = this._uburnSum / this._uburnN / 1e3;
    this._uburnSum = 0; this._uburnN = 0;
    this.breakdownsPerPeriod = this._bdCount;
    this.peakCurrent = this.peakCurrentHist;
    this._bdCount = 0;
    this.energyAccum = 0;
    this._slopeOnSum = this._slopeOffSum = 0;
    this._slopeOnN = this._slopeOffN = 0;
  }

  advance(wallClockBudgetMs = 12) {
    const t0 = performance.now();
    let steps = 0;
    const tStart = this.t;
    const target = this.targetSimTime;
    while (this.t < target) {
      this.step();
      steps++;
      if ((steps & 255) === 0 && performance.now() - t0 > wallClockBudgetMs) break;
      if (steps > 4e6) break;
    }
    return { steps, simTime: this.t - tStart };
  }

  // ------------------------------------------------------------- 1D painting
  _refreshProfiles() {
    this._profileDirty = false;
    const nx = this.nx, d = this.d, x = this.x;
    const n = this.n;
    const pol = this.Ugap >= 0 ? 1 : -1;
    const xc = this.frontPos * d;
    const w = (0.06 + 0.10 / (1 + this.glow)) * d;
    const nePeak = this.ne * (2 + 6 * clamp(this.glow, 0, 1));
    const sheath = 0.05 * d;

    for (let i = 0; i < nx; i++) {
      const xi = x[i];
      const g = sexp(-((xi - xc) * (xi - xc)) / (2 * w * w));
      // electrons: depleted in the cathode sheath
      const cathDist = pol > 0 ? d - xi : xi;
      const sh = 1 - sexp(-(cathDist * cathDist) / (2 * sheath * sheath));
      n.e[i] = Math.max(1e10, nePeak * g * sh + this.ne * 0.05);
      n.O2p[i] = Math.max(1e10, nePeak * 0.35 * (g * 0.6 + 0.25) + this.ne * 0.15);
      n.O4p[i] = n.O2p[i] * 0.55;
      n.Om[i] = Math.max(1e9, this.nOm * (0.4 + 0.6 * g));
      n.O2m[i] = Math.max(1e9, this.nO2m * (0.7 + 0.3 * (1 - g)));
      n.O[i] = Math.max(1e10, this.nO * (0.7 + 0.5 * g));
      n.O3[i] = Math.max(1e9, this.nO3 * (0.85 + 0.3 * (1 - g)));
      n.O2a[i] = Math.max(1e9, this.nO2a * (0.8 + 0.4 * g));
    }

    // rho and self-consistent 1D field:  E(x) = E0 + (1/eps0) * int_0^x rho dx'
    const rho = this.rho, E = this.E, phi = this.phi;
    let acc = 0, integ = 0;
    for (let i = 0; i < nx; i++) {
      rho[i] = QE * (n.O2p[i] + n.O4p[i] - n.e[i] - n.Om[i] - n.O2m[i]);
      const dx = this.xFaces[i + 1] - this.xFaces[i];
      acc += (rho[i] / EPS0) * dx;
      E[i] = acc;
      integ += E[i] * dx;
    }
    const E0 = (this.Ugap - integ) / d;
    let ph = 0;
    for (let i = 0; i < nx; i++) {
      E[i] += E0;
      const dx = this.xFaces[i + 1] - this.xFaces[i];
      ph -= E[i] * dx;
      phi[i] = ph;
      const absE = Math.abs(E[i]);
      this.EN[i] = (absE / this.N) * 1e21;
      const mue = 0.05 * (2.45e25 / this.N);
      const alpha = absE > 1 ? 4.5e-20 * this.N * sexp(-9.0e-19 * this.N / absE) : 0;
      this.ionizRate[i] = n.e[i] * alpha * mue * absE;
    }

    const s = this._state;
    s.t = this.t; s.dt = this.dt;
    s.sigmaL = -this.Qcond / this.A;
    s.sigmaR = this.Qcond / this.A;
    s.Uapp = this._uapp(this.t);
    s.Ugap = this.Ugap;
    s.current = this.Itot;
    s.condCurrent = this.Icond;
    s.dispCurrent = this.Idisp;
    s.charge = this.Q;
    s.power = this.periodStats.powerW;
    s.o3ppm = (this.nO3 / this.N) * 1e6;
    s.maxEN = this.maxEN;
    s.peakCurrent = this.peakCurrentHist;
    s.breakdownsPerPeriod = this.breakdownsPerPeriod;
  }
}

export default DBDSolver;
