# PHYSICS.md — 1D fluid model of a dielectric barrier discharge in pure O₂ at 1 atm

Complete, numerically explicit dataset for the DBD-O₂ simulator (`/p/dbd-o2`).
All quantities are **SI** unless a table column says otherwise.
Fit variable everywhere: `x = E/N` in **townsend** (1 Td = 1e-21 V·m²).

Reference gas state (defaults): `p = 101325 Pa`, `T = 300 K`,
`N = p/(k_B T) = 2.446e25 m⁻³` (rounded to **2.45e25 m⁻³** below).
Conversion: `E [V/m] = x [Td] * 1e-21 * N` → **1 Td ↔ 2.45e4 V/m = 0.245 kV/cm** at 1 atm/300 K.

---

## 0. Geometry, circuit, and the numbers that fall out of it

| Quantity | Symbol | Value |
|---|---|---|
| Gap | `d_gap` | 1.0e-3 m |
| Barrier 1 / 2 | `d_d1`, `d_d2` | 0.5e-3 m each |
| Relative permittivity (Al₂O₃) | `eps_r` | 9 |
| Electrode area | `A` | 1.0e-4 m² (1 cm²) |
| Dielectric capacitance (both barriers in series) | `C_d = eps0*eps_r*A/(d_d1+d_d2)` | **7.97e-12 F = 7.97 pF** |
| Gas-gap capacitance | `C_g = eps0*A/d_gap` | **8.85e-13 F = 0.885 pF** |
| Cell (off) capacitance | `C_cell = C_d C_g/(C_d+C_g)` | **7.97e-13 F = 0.797 pF** |

These two capacitances are exactly the two slopes of the Q–V Lissajous parallelogram
(Manley method): the "discharge-off" branches have slope `C_cell`, the "discharge-on"
branches have slope `C_d`.

---

## 1. Species set

9 species + the background gas. Minimal but complete for O₂ DBD physics
(electronegativity, cluster ions, ozone, metastable-driven detachment).

`mu` given at **1 atm, 300 K** (N = 2.45e25 m⁻³). Reduced mobility `K0` is quoted at
`N0 = 2.687e25 m⁻³` (Loschmidt); scaling in the code must be `mu(N) = K0 * N0/N`
(and `mu(N,T)` if T is ever varied). `N0/N = 1.0967` at 300 K/1 atm.

| # | Species | m [kg] | q [C] | K0 [m²/V/s] @N0 | mu @1atm,300K [m²/V/s] | D [m²/s] | n₀ [m⁻³] |
|---|---|---|---|---|---|---|---|
| 0 | e | 9.109e-31 | −1.602e-19 | — (field dep.) | `mu_e(x)` §3 | `D_e(x)` §3 | 1.0e13 |
| 1 | O₂⁺ | 5.313e-26 | +1.602e-19 | 2.20e-4 | 2.41e-4 | 6.23e-6 | 1.0e13 |
| 2 | O₄⁺ | 1.063e-25 | +1.602e-19 | 2.16e-4 | 2.37e-4 | 6.13e-6 | 0 |
| 3 | O⁻ | 2.657e-26 | −1.602e-19 | 3.20e-4 | 3.51e-4 | 9.07e-6 | 0 |
| 4 | O₂⁻ | 5.313e-26 | −1.602e-19 | 2.16e-4 | 2.37e-4 | 6.13e-6 | 0 |
| 5 | O₃⁻ | 7.970e-26 | −1.602e-19 | 2.40e-4 | 2.63e-4 | 6.80e-6 | 0 |
| 6 | O (³P) | 2.657e-26 | 0 | — | 0 | 2.0e-5 | 1.0e15 |
| 7 | O₃ | 7.970e-26 | 0 | — | 0 | 1.5e-5 | 0 |
| 8 | O₂(a¹Δg) | 5.313e-26 | 0 | — | 0 | 2.0e-5 | 0 |
| bg | O₂ (ground) | 5.313e-26 | 0 | — | — | — | 2.45e25 (fixed) |

**Ion diffusion** by Einstein relation at gas temperature: `D_i = mu_i * k_B T / e = mu_i * 0.02585 V`
(0.02585 V = kT/e at 300 K). This is exact for ions in the low-field limit and adequate
here because ion diffusion is negligible against ion drift in a DBD.

**Neutral diffusion** — binary diffusion coefficients in O₂ at 1 atm/300 K:
D(O–O₂) ≈ 2.0e-5, D(O₃–O₂) ≈ 1.5e-5, D(O₂(a)–O₂) ≈ 2.0e-5 m²/s.
Scale as `D ∝ T^1.75 / p` if needed.

Ion mobility sources: Ellis et al., *At. Data Nucl. Data Tables* **17** (1976) & **22** (1978)
(ion mobility compilation); Viehland–Mason. O₄⁺ and O₃⁻ values are the standard
cluster-ion values used in atmospheric-pressure O₂/air models (Kossyi et al. 1992).

> **Justification for including O₄⁺, O₃⁻, O₂(a)**: at 1 atm the three-body conversion
> O₂⁺+2O₂→O₄⁺+O₂ has a time constant `1/(k*N²) = 1/(2.4e-42*(2.45e25)²) ≈ 0.7 ns` —
> O₄⁺ is the dominant positive ion in the afterglow within nanoseconds, and it recombines
> with electrons ~70× faster than O₂⁺. O⁻ converts to O₃⁻ in `1/(1.1e-42*N²) ≈ 1.5 ns`,
> so O₃⁻ carries most of the negative charge between pulses. O₂(a) is the main
> detachment partner that lets the discharge re-ignite (memory effect in the volume).

---

## 2. Reactions

`x = E/N` in Td, `T` gas temperature [K], `Te` electron temperature [K],
`eps` mean electron energy [eV]. Use `Te = (2/3)*eps*11604.5`.

### 2.1 Electron-impact (LFA, functions of E/N)

Ionization is anchored to Morrow's Townsend coefficient for O₂
(R. Morrow, *Phys. Rev. A* **32** (1985) 1799; Morrow & Lowke, *J. Phys. D* **30** (1997) 614):

```
alphaN(x) = (x < 150) ? 6.619e-21*Math.exp(-559.3/x)     // m^2
                      : 2.000e-20*Math.exp(-724.8/x);
```
The two branches are continuous at x = 149.7 Td (that is where they cross), so the
switch point is 150 Td, not the sometimes-quoted 105.

Rate coefficient from the coefficient: `k = (alpha/N) * w_e`, `w_e = muN(x)*x*1e-21`.

| ID | Reaction | ΔE [eV] | k [m³/s] (or m⁶/s) | Source |
|---|---|---|---|---|
| E1 | e + O₂ → 2e + O₂⁺ | 12.06 | `alphaN(x)*muN(x)*x*1e-21` | Morrow 1985 / Phelps |
| E2 | e + O₂ → O⁻ + O | 4.2 (thr.) | `8.0e-18*Math.exp(-160/x)/(1+x/3000)` | fit to Phelps/BOLSIG+ |
| E3 | e + O₂ + O₂ → O₂⁻ + O₂ | ~0 | `3.0e-43*Math.exp(-x/60) + 5.0e-45` **[m⁶/s]** | Chanin–Phelps–Biondi |
| E4 | e + O₂ → e + O + O | 6.0 | `4.23e-15*Math.exp(-5.56/Te)`, `Te = (2/3)<eps>` | Gudmundsson & Thorsteinsson 2007 **[ERRATA §D]** |
| E5 | e + O₂ → e + O₂(a¹Δg) | 0.98 | `1.37e-15*Math.exp(-2.14/Te)`, `Te = (2/3)<eps>` | Gudmundsson & Thorsteinsson 2007 **[ERRATA §D]** |
| E6 | e + O₂⁺ → O + O | — | `2.0e-13*Math.pow(300/Te, 0.7)` | Kossyi 1992 |
| E7 | e + O₄⁺ → O₂ + O₂ | — | `1.4e-12*Math.pow(300/Te, 0.5)` | Kossyi 1992 |
| E8* | e + O₃ → O⁻ + O₂ | — | `1.0e-15` (optional, [O₃] small during pulse) | Kossyi 1992 |

\* optional; include only if you run many periods and O₃ accumulates above ~1e22 m⁻³.

**Sanity of the fits** (values in m³/s, m⁶/s for E3):

| x [Td] | E1 (ioniz.) | E2 (diss. att.) | E3·N (eff. 2-body) | E4 (diss.) | E5 (a¹Δ) |
|---|---|---|---|---|---|
| 30 | 3.3e-24 | 3.8e-20 | 4.6e-18 | 1.1e-17 | 1.4e-16 |
| 50 | 7.6e-21 | 3.2e-19 | 3.3e-18 | 5.9e-17 | 2.6e-16 |
| 100 | 2.8e-18 | 1.6e-18 | 1.5e-18 | 2.9e-16 | 4.9e-16 |
| 150 | 2.4e-17 | 2.6e-18 | 7.3e-19 | 5.1e-16 | 6.1e-16 |
| 300 | 4.2e-16 | 4.3e-18 | 1.7e-19 | 1.0e-15 | 7.9e-16 |
| 1000 | 4.9e-15 | 5.1e-18 | 1.2e-19 | 2.0e-15 | 1.0e-15 |

E1–E3 are within a factor ≲2 of BOLSIG+/Phelps-O₂ over 20–1000 Td.

**E4/E5 were corrected after the first full run** (ERRATA §D). The original fits
(`6.0e-17*exp(-350/x)` and `2.0e-17*exp(-28/x)`) gave 5.8e-18 and 1.7e-17 m³/s at
150 Td — 1.5–2 orders of magnitude low, i.e. ~3800 eV of discharge energy per
dissociation event. Ozone came out at 1.3 ppm/period instead of the ~80 ppm/period
required by item 10 of §6, a factor ~60 miss, while every *electrical* quantity
(power, transferred charge, Lissajous slopes) was already correct — the giveaway
that the error was in the neutral-production channel, not in the discharge itself.
With the replacement the LFA now costs ~40 eV of discharge energy per dissociation
at 150 Td, which is what the empirical 100 g(O₃)/kWh yield implies, and the
simulated ozone lands within a factor 1.8 of it.

### 2.2 Ion conversion / clustering

| ID | Reaction | k | Source |
|---|---|---|---|
| C1 | O₂⁺ + O₂ + O₂ → O₄⁺ + O₂ | `2.4e-42*Math.pow(300/T,3.2)` m⁶/s | Kossyi 1992 |
| C2 | O₄⁺ + O₂ → O₂⁺ + O₂ + O₂ | `3.3e-12*Math.pow(300/T,4)*Math.exp(-5030/T)` m³/s (≈1.7e-19 @300 K) | Kossyi 1992 |
| C3 | O⁻ + O₂ + O₂ → O₃⁻ + O₂ | `1.1e-42` m⁶/s | Kossyi 1992 |
| C4 | O⁻ + O₃ → O₃⁻ + O | `5.3e-16` m³/s | Kossyi 1992 |
| C5 | O₂⁻ + O₃ → O₃⁻ + O₂ | `4.0e-16` m³/s | Kossyi 1992 |
| C6 | O₃⁻ + O → O₂⁻ + O₂ | `2.5e-16` m³/s | Kossyi 1992 |

### 2.3 Detachment (the "volume memory" that lets the discharge re-strike)

| ID | Reaction | k [m³/s] | Source |
|---|---|---|---|
| D1 | O⁻ + O₂(a¹Δg) → O₃ + e | `3.0e-16` | Kossyi 1992 |
| D2 | O⁻ + O → O₂ + e | `5.0e-16` | Kossyi 1992 |
| D3 | O₂⁻ + O₂(a¹Δg) → O₂ + O₂ + e | `2.0e-16` | Kossyi 1992 |
| D4 | O₂⁻ + O → O₃ + e | `1.5e-16` | Kossyi 1992 |
| D5 | O⁻ + O₂ → O₃ + e | `5.0e-21` | Kossyi 1992 |
| D6 | O₂⁻ + O₂ → O₂ + O₂ + e (collisional, field-driven) | `2.7e-16*Math.sqrt(Teff/300)*Math.exp(-5590/Teff)` | Kossyi 1992 |
| D7 | O₃⁻ + O → O₂ + O₂ + e | `1.0e-17` | Kossyi 1992 |

Effective ion temperature (Wannier) for D6, with `m_g = 5.313e-26` kg (O₂) and
`v_d = mu_i*E`:
```
Teff = T + (m_g * v_d*v_d) / (3*1.380649e-23);
```
At E = 3.78e6 V/m, v_d(O₂⁻) = 896 m/s → Teff = 300 + 1030 = 1330 K → k(D6) = 8.8e-18 m³/s,
i.e. detachment frequency `k*N = 2.2e8 s⁻¹`. This term is what keeps a Townsend-like
O₂ DBD from being completely electron-starved at high field — do **not** drop it.

### 2.4 Ion–ion recombination

At 1 atm the three-body channel dominates. Use for **every** (negative, positive) pair:

```
k_ii = 2.0e-13 + 2.0e-37*N;      // m^3/s ; = 5.1e-12 at N = 2.45e25
```
(2e-7 cm³/s two-body + 2e-25 cm⁶/s three-body, Kossyi 1992; the sum reproduces the
measured atmospheric-pressure value ~2e-6 cm³/s.)

| ID | Reaction | k |
|---|---|---|
| R1 | O⁻ + O₂⁺ → O + O₂ | `k_ii` |
| R2 | O₂⁻ + O₂⁺ → O₂ + O₂ | `k_ii` |
| R3 | O₃⁻ + O₂⁺ → O₃ + O₂ | `k_ii` |
| R4 | O⁻ + O₄⁺ → O + O₂ + O₂ | `k_ii` |
| R5 | O₂⁻ + O₄⁺ → 3 O₂ | `k_ii` |
| R6 | O₃⁻ + O₄⁺ → O₃ + O₂ + O₂ | `k_ii` |

### 2.5 Neutral chemistry (ozone)

| ID | Reaction | k | Source |
|---|---|---|---|
| N1 | O + O₂ + M → O₃ + M | `6.0e-46*Math.pow(300/T,2.4)` m⁶/s, M = O₂ | Atkinson/JPL |
| N2 | O + O + M → O₂ + M | `3.0e-46*Math.pow(300/T,1.0)` m⁶/s | Kossyi 1992 |
| N3 | O + O₃ → O₂ + O₂ | `8.0e-18*Math.exp(-2060/T)` m³/s (8.3e-21 @300 K) | JPL |
| N4 | O₃ + O₂(a) → O + O₂ + O₂ | `5.2e-17*Math.exp(-2840/T)` m³/s | Kossyi 1992 |
| N5 | O₂(a) + O₂ → O₂ + O₂ | `2.2e-24` m³/s (τ ≈ 19 ms) | Kossyi 1992 |
| N6 | O₂(a) + O → O₂ + O | `7.0e-22` m³/s | Kossyi 1992 |
| N7 | O₃ + O₂ → O + O₂ + O₂ | `7.3e-16*Math.exp(-11400/T)` m³/s (negligible @300 K) | Kossyi 1992 |

Ozone formation time constant: `1/(k_N1 * N²) = 1/(6.0e-46*(2.45e25)²) = 2.8 µs` — i.e. O
atoms produced in one current pulse become O₃ well before the next half period (50 µs at
10 kHz). This is the correct and expected ordering.

---

## 3. Electron transport

Primary representation: **log–log interpolated table** (BOLSIG+ / Phelps O₂ cross sections,
N = 2.45e25 m⁻³). Analytic fallbacks are given below it.

| x [Td] | mu_e·N [1/(m·V·s)] | w_e [m/s] | eps_mean [eV] | D_e·N [1/(m·s)] |
|---|---|---|---|---|
| 0.1 | 6.50e24 | 6.5e2 | 0.035 | 1.52e23 |
| 0.3 | 6.20e24 | 1.86e3 | 0.045 | 1.86e23 |
| 1 | 5.00e24 | 5.0e3 | 0.12 | 4.00e23 |
| 3 | 3.90e24 | 1.17e4 | 0.30 | 7.80e23 |
| 10 | 3.00e24 | 3.0e4 | 0.75 | 1.50e24 |
| 20 | 2.40e24 | 4.8e4 | 1.10 | 1.76e24 |
| 30 | 2.05e24 | 6.15e4 | 1.40 | 1.91e24 |
| 50 | 1.65e24 | 8.25e4 | 1.95 | 2.15e24 |
| 70 | 1.40e24 | 9.80e4 | 2.45 | 2.29e24 |
| 100 | 1.15e24 | 1.15e5 | 3.10 | 2.38e24 |
| 150 | 9.90e23 | 1.49e5 | 3.95 | 2.61e24 |
| 200 | 9.00e23 | 1.80e5 | 4.65 | 2.79e24 |
| 300 | 7.90e23 | 2.37e5 | 5.80 | 3.05e24 |
| 500 | 6.60e23 | 3.30e5 | 7.60 | 3.34e24 |
| 700 | 5.90e23 | 4.13e5 | 9.00 | 3.54e24 |
| 1000 | 5.10e23 | 5.10e5 | 10.8 | 3.67e24 |

`D_e = mu_e * (2/3) * eps_mean` (characteristic energy ≈ ⅔·mean energy; for the strongly
non-Maxwellian O₂ EEDF this is accurate to ~15 % and always errs on the safe/diffusive side).

**Analytic fallbacks** (use outside the table or if you want a table-free build):
```
muN_fit(x)  = 1.10e24 * Math.pow(100/x, 0.28);           // ±25% over 10..1000 Td
eps_fit(x)  = 0.045 + 0.115*Math.pow(x, 0.72);           // eV, ±15% over 1..1000 Td
```
Check of `eps_fit`: x=10 → 0.65; x=100 → 3.11; x=1000 → 15.4 (this one degrades above
~500 Td; prefer the table).

**Townsend coefficients** (report these in the UI — they are what makes the physics legible):
```
alphaN(x) = (x<150) ? 6.619e-21*Math.exp(-559.3/x) : 2.0e-20*Math.exp(-724.8/x);   // m^2
etaN(x)   = (kE2(x) + kE3(x)*N) / (muN(x)*x*1e-21);                                 // m^2
alpha_eff = (alphaN(x) - etaN(x)) * N;                                              // 1/m
```

| x [Td] | α/N [m²] | η/N [m²] | α_eff [1/m] |
|---|---|---|---|
| 50 | 9.2e-26 | 4.41e-23 | −1.08e3 |
| 80 | 6.1e-24 | 2.99e-23 | −5.84e2 |
| **101.3** | **2.63e-23** | **2.63e-23** | **≈ 0** |
| 120 | 6.26e-23 | 2.44e-23 | 9.37e2 |
| 150 | 1.59e-22 | 2.25e-23 | 3.35e3 |
| **154.4** | **1.81e-22** | **2.24e-23** | **3.88e3 (= ln(1+1/γ)/d)** |
| 200 | 5.34e-22 | 2.09e-23 | 1.26e4 |
| 300 | 1.79e-21 | 1.87e-23 | 4.33e4 |

---

## 4. Boundary conditions

### 4.1 Secondary electron emission
`gamma = 0.02` (default; slider 0.005–0.10) for **all positive ions** hitting either
dielectric surface. Rationale: γ for clean metals in O₂ is 0.005–0.02; oxide surfaces
(Al₂O₃, MgO) have high secondary yields and reported effective γ in atmospheric O₂/air
DBDs on alumina falls in 0.01–0.05 once surface charging and photoemission are lumped in.
0.02 sits in the middle and reproduces the expected ~3.8 kV gap breakdown for 1 mm
(see §6). Photoemission and metastable-induced emission are folded into the same γ;
do not add a separate photoemission term unless you also add photon transport.

```
Gamma_e,wall_emitted = gamma * sum_over_positive_ions( Gamma_i,wall )
```

### 4.2 Fluxes at the gas/dielectric interface
For a wall with outward normal `n`:

- Positive ions and negative ions/electrons **drifting toward** the wall: full absorption,
  flux = drift + thermal.
- Species drifting **away**: only the thermal (random-walk) half-flux escapes.

```
a = (drift velocity points into the wall) ? 1 : 0;
Gamma_e = a*|mu_e*E|*n_e + 0.25*v_th_e*n_e - gamma*Gamma_pos;
Gamma_i = a*|mu_i*E|*n_i + 0.25*v_th_i*n_i;
v_th_e  = Math.sqrt(8*1.602e-19*(2/3)*eps_mean/(Math.PI*9.109e-31));   // ~1.1e6 m/s @3 eV
v_th_i  = Math.sqrt(8*1.380649e-23*T/(Math.PI*m_i));                   // ~4.4e2 m/s
```
Reflection coefficient for charged species: 0 (full absorption/neutralisation).

### 4.3 Neutral species at the wall
Loss probabilities (sticking/recombination coefficients) on Al₂O₃:
`gamma_O = 2e-3`, `gamma_O3 = 1e-5` (≈0), `gamma_O2a = 3e-4`.
Wall flux `Gamma = 0.25*gamma_s*v_th*n`. O-atom wall loss frequency for d=1 mm:
`gamma_O*v_th_O/(2d) ≈ 2e-3*630/(2e-3) ≈ 6e2 s⁻¹` — negligible against the 3.6e5 s⁻¹
gas-phase ozone-forming loss, as it should be at atmospheric pressure.

### 4.4 Surface charge and Poisson
Surface charge density on each dielectric face:
```
dSigma/dt = e*( Gamma_pos_to_wall - Gamma_e_to_wall - Gamma_neg_to_wall )
```
(net **positive** flux increases σ; the emitted secondary flux is already subtracted in
`Gamma_e`). Poisson with the jump condition:
```
d/dz ( eps0*eps_r(z) * dV/dz ) = -rho(z)      in gas and dielectrics
eps0*eps_r*E|_diel - eps0*E|_gas = -Sigma     at each interface
```
Solve the whole stack (dielectric–gas–dielectric) in one tridiagonal system with
`eps_r = 9` in the barriers and 1 in the gas, Dirichlet `V = U(t)` at one metal and
`V = 0` at the other. With an external series resistor `R`:
`U_gas_stack = U_src(t) − R*I(t)`, `I = A*(J_cond + eps0*dE/dt)` evaluated at any plane.

Applied waveform: `U_src(t) = U0*Math.sin(2*Math.PI*f*t)`, defaults `U0 = 1.0e4 V`,
`f = 1.0e4 Hz`.

---

## 5. Seed ionization

Physical background at sea level: cosmic rays + natural radioactivity ionize at
`S_cosmic ≈ 1e7 m⁻³s⁻¹`; balanced against ion–ion recombination `k_ii = 5.1e-12`
this gives a steady-state `n ≈ sqrt(S/k_ii) ≈ 1.4e9 m⁻³`. That is too small to seed the
*first* half period in a reasonable number of time steps, and it is **not** what seeds a
real DBD — residual charge from the previous discharge does.

Recommended, numerically robust prescription:

```
n_e(0) = n_O2p(0) = 1.0e13 m^-3    (uniform in the gap)   // "residual from previous cycle"
S_bg   = 1.0e13 m^-3 s^-1          (uniform volume source of e + O2+)
n_floor= 1.0e10 m^-3               (clamp applied to e, O2+, O4+ each step)
```
The volume source contributes `S_bg * T_half = 5e8 m⁻³` per half period at 10 kHz —
five orders of magnitude below the discharge densities, so it perturbs nothing after the
first ignition but guarantees ignition happens. Expose `n_e0` and `S_bg` as UI sliders
(1e11…1e15) so the user can see that the *steady* cycle is independent of them.

---

## 6. Validation checklist (numbers the simulation must hit)

1. **Static Townsend breakdown criterion** α = η crosses at **E/N = 101.3 Td**,
   i.e. **E = 2.48e6 V/m = 24.8 kV/cm**. (Accepted value for pure O₂: 100–110 Td.)
2. **Self-sustaining criterion for the 1 mm gap**, `(α−η)d = ln(1+1/γ) = 3.93` for γ=0.02,
   satisfied at **E/N = 154.4 Td → E_br = 3.78e6 V/m = 37.8 kV/cm → U_gap,br = 3.78 kV**.
   The simulated gap voltage at the instant of each current pulse must be 3.5–4.2 kV
   and must **not** grow with U0 (that is the definition of the burning voltage).
3. **Minimum applied amplitude for ignition**
   `U_min = U_br*(1 + C_g/C_d) = 3.8*(1+0.111) = 4.2 kV`. Below ~4.2 kV amplitude the
   simulation must show a pure capacitive (no discharge) response — a *degenerate*
   Lissajous line of slope `C_cell = 0.797 pF`, no parallelogram.
4. **Lissajous slopes**: fitting the Q–V figure must return `C_cell = 0.80 ± 0.03 pF`
   (off branches) and `C_d = 7.97 ± 0.3 pF` (on branches).
5. **Manley power** at U0 = 10 kV, f = 10 kHz:
   `P = 4 f C_d U_min (U0 − U_min) = 4·1e4·7.97e-12·4220·5780 ≈ **7.8 W**` (≈ 7.8 W/cm²,
   ≈ 7.8e7 W/m³). The numerically integrated `P = f∮U dQ` must match within ~20 %.
6. **Charge transferred per half period**
   `ΔQ = 2 C_d (U0 − U_min) ≈ 2·7.97e-12·5780 ≈ **92 nC**` (per cm²). The time integral of
   the conduction current over one pulse must reproduce this.
7. **Current pulse**: one (or a few) pulse per half period, **FWHM 20–200 ns**, peak
   conduction current density **1e2–1e4 A/m² (10 mA–1 A per cm²)**. Consistency check:
   1 A/cm² × 100 ns = 100 nC ≈ item 6. *Note the 1D caveat:* a 1D model fires the whole
   1 cm² simultaneously, so it reports the *area-averaged* filament current; a real
   microdischarge is ~0.1 mm² carrying 0.1–1 nC, i.e. ~10³ filaments per cm² per half cycle.
8. **Phase shift / memory effect**: ignition must occur **before** the voltage maximum in
   steady state (typically 30–70° after the zero crossing), and the surface charge must
   reverse sign every half period with `|σ| ≈ ΔQ/(2A) ≈ 4.6e-4 C/m²`. Corresponding
   memory field `σ/eps0 ≈ 5.2e7 V/m` in the gap-side vacuum field — this is what self-quenches
   the discharge and it must appear in the gap-voltage trace as a plateau at ~U_br.
9. **Electronegativity**: at the current-pulse peak the negative-ion density must be
   comparable to or larger than n_e, and within ~1 µs after the pulse
   `n(O₃⁻) + n(O₂⁻) ≫ n_e` (electron attachment time at low field
   `1/(k_E3·N²) = 5.6 ns` at 1 Td; three-body attachment must visibly kill the
   electrons between pulses).
10. **Ozone**: after a few periods `[O₃] ≈ 1e21–1e22 m⁻³` in the gap. Cross-check with the
    empirical yield: 100 g(O₃)/kWh = 3.5e17 molecules/J; at 7.8e7 W/m³ this gives
    `d[O₃]/dt ≈ 2.7e25 m⁻³s⁻¹` → **2.7e21 m⁻³ per period** at 10 kHz. Simulated yield within
    a factor 3 of this is a pass. (No gas flow in the 1D model, so O₃ accumulates
    monotonically; it saturates near 1e23–1e24 m⁻³ where N3/N4 destruction balances.)
11. **Timescale hierarchy** (assert these in the code, they catch unit errors instantly):
    electron transit `d/w_e = 1e-3/1.2e5 = 8 ns`; ion transit `d/(mu_i E) = 1e-3/910 = 1.1 µs`;
    O₂⁺→O₄⁺ conversion 0.7 ns; O→O₃ 2.8 µs; half period 50 µs.

---

## 7. Numerics notes (not requested, but these decide whether the above is reachable)

- **Dielectric (Maxwell) relaxation time** `tau_M = eps0/(e n_e mu_e)`:
  1.2e-7 s at n_e = 1e16, **1.2e-9 s at n_e = 1e18**. An explicit Poisson-coupled scheme
  needs `dt < tau_M`, so use the **semi-implicit (Ventzek/Boeuf) Poisson** correction
  `∇·((eps + dt*sigma_p)∇V) = -rho` with `sigma_p = e(n_e mu_e + Σ n_i mu_i)`, otherwise
  dt collapses to sub-ns during the pulse.
- **Fluxes**: Scharfetter–Gummel exponential scheme for all charged species (it is the only
  cheap scheme that stays positive at the Peclet numbers here, `mu E dz/D ~ 50` for ions).
- **Grid**: ≥ 400 cells across 1 mm, refined by ~5× within 50 µm of each dielectric
  (the cathode sheath in a 1 mm atmospheric gap is 20–80 µm).
- **Time stepping**: adaptive, `dt = min(0.5*dz/w_e, 0.2/max(nu_ion), 0.2*tau_M_semi)`;
  expect ~0.1–1 ns during the pulse and 10–100 ns between pulses.
- **Neutral chemistry** (O, O₃, O₂(a)) is stiff-but-slow: sub-cycle it with its own
  implicit/exponential integrator on a 10–100 ns step, decoupled from the plasma step.

---

## 8. References

1. R. Morrow, "Theory of negative corona in oxygen", *Phys. Rev. A* **32** (1985) 1799 —
   α/N, η/N for O₂.
2. R. Morrow, J. J. Lowke, "Streamer propagation in air", *J. Phys. D* **30** (1997) 614.
3. I. A. Kossyi, A. Yu. Kostinsky, A. A. Matveyev, V. P. Silakov, "Kinetic scheme of the
   non-equilibrium discharge in nitrogen–oxygen mixtures", *Plasma Sources Sci. Technol.*
   **1** (1992) 207 — the source for essentially all heavy-particle rates above.
4. A. V. Phelps, JILA/LXCat O₂ cross-section set (www.lxcat.net) — electron transport and
   electron-impact rates; BOLSIG+ (Hagelaar & Pitchford, *PSST* **14** (2005) 722).
5. L. G. Chanin, A. V. Phelps, M. A. Biondi, "Measurement of the attachment of low-energy
   electrons to oxygen molecules", *Phys. Rev.* **128** (1962) 219 — three-body attachment.
6. H. W. Ellis et al., *At. Data Nucl. Data Tables* **17** (1976) 177 and **22** (1978) 179 —
   ion mobilities.
7. U. Kogelschatz, "Dielectric-barrier discharges: their history, discharge physics, and
   industrial applications", *Plasma Chem. Plasma Process.* **23** (2003) 1 — DBD phenomenology,
   filament charge/duration, ozone yields.
8. T. C. Manley, "The electric characteristics of the ozonator discharge",
   *Trans. Electrochem. Soc.* **84** (1943) 83 — Q–V Lissajous / power formula.
9. JPL Publication 19-5, "Chemical Kinetics and Photochemical Data" — O + O₂ + M and
   O + O₃ rate coefficients.
