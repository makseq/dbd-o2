# dbd-o2 — DBD in pure oxygen: a 2D axisymmetric fluid model with photoprocesses

Русская версия: [README.ru.md](README.ru.md)

**Live site (English / Russian):** https://barier.makseq.com — its deployed copy is in [`site/`](site/) (the precomputed run data, ~51 MB, is not in the repository; it is reproducible with `sim2d/run.mjs`). Built by Max Tkachenko together with Claude Code.


A self-consistent fluid model of a volumetric **dielectric barrier discharge (DBD)**
in pure O₂ at 1 atm. The main part is a **2D axisymmetric (r,z) solver with photoprocesses**
(offline computation in node.js, the result being frames on disk) and a **web player** that
replays those frames. The one-dimensional model everything started with is kept and described in
[`docs/README_1D.md`](docs/README_1D.md).

Pure ES2022. Zero dependencies, zero build, zero CDNs — neither in the solver nor in the player.

```
      r
  R=0.5mm ┌───────┬───────────────────────┬───────┐
          │ Al₂O₃ │      O₂ gas gap       │ Al₂O₃ │   mirror wall (Neumann)
          │ ε_r=9 │  1.0 mm, 1 atm, 300 K │ ε_r=9 │
    r=0   ├───────┼───────────────────────┼───────┤ ← axis of symmetry (Neumann)
          │0.5 mm │                       │0.5 mm │
          └───────┴───────────────────────┴───────┘
        z=0                                      z=2 mm
     metal                                    metal
   φ = U₀ sin(2πft)                              φ = 0
   U₀ = 10 kV, f = 10 kHz
```

Seeding: a 1e13 m⁻³ background **plus** a Gaussian spot on the axis at the dielectric surface
(σ_r = 40 µm). Without the spot no filament is born (`ERRATA B6`) — the model simply
yields a diffuse Townsend discharge.

---

## Quick start

```bash
./serve.sh                 # http://127.0.0.1:8777/player/  — 2D player on ready-made data
./serve.sh 9000            # a different port
```

ES modules do not load from `file://` (CORS), so an http server is required; `serve.sh`
uses `python3 -m http.server`, and if it is absent, a built-in node server.
Ready-made runs already lie in `data/`, nothing needs to be computed.

| URL | what opens |
|---|---|
| `/player/` | **2D player** (r,z): field, σ(r), waveforms, Lissajous, photoprocesses |
| `/player/plots-demo.html` | standalone demo page of the plots (plots only, without the scene) |
| `/index.html` | the old **1D** UI with the live 1D solver (see `docs/README_1D.md`) |

### Your own run

```bash
# one breakdown in the filamentary regime, 32×120 mesh, photoprocesses enabled
node sim2d/run.mjs --preset filament --nr 32 --nz 120 --nzDiel 12 \
     --betaR 1.4 --betaG 1.5 --U0kV 10 --freqKHz 10 --photo on \
     --periods 0.14 --out data/my-run

node sim2d/summarize.mjs data                 # acceptance table over all runs
node sim2d/compare-photo.mjs data/run-default data/run-nophoto
```

Progress goes to **stderr**, and stdout carries only the final JSON (it can be piped to `| jq`).
The runner understands `SIGINT/SIGTERM/SIGHUP` (a soft stop with the manifest written out),
`--maxSeconds`, `--maxSteps`, `--checkpointSec`, `--maxBytesFull/--maxBytesCompact`.
The full list is `node sim2d/run.mjs --help`.

Exactly the five runs that are in `data/`:

```bash
G="--preset filament --nr 32 --nz 120 --nzDiel 12 --betaR 1.4 --betaG 1.5"
node sim2d/run.mjs $G --out data/run-default --U0kV 10 --freqKHz 10 --photo on  --periods 0.14
node sim2d/run.mjs $G --out data/run-nophoto --U0kV 10 --freqKHz 10 --photo off --periods 0.14
node sim2d/run.mjs $G --out data/run-low     --U0kV 6  --freqKHz 10 --photo on  --periods 0.20
node sim2d/run.mjs $G --out data/run-high    --U0kV 14 --freqKHz 10 --photo on  --periods 0.10
node sim2d/run.mjs $G --out data/run-fast    --U0kV 10 --freqKHz 30 --photo on  --periods 0.40
```

### Tests

```bash
node sim2d/tests/physics.test.mjs    # 61
node sim2d/tests/poisson.test.mjs    # 14
node sim2d/tests/solver.test.mjs     # 153
node sim2d/tests/recorder.test.mjs   # 12
node sim2d/tests/photo.test.mjs      # 18
node test/solver.test.mjs            # 1D core, 55 tests (~70 s)
```

---

## What is computed

### Physics

The drift-diffusion approximation (LFA: the coefficients are functions of the local `E/N`),
self-consistent with the Poisson equation, plus surface charge on both barriers.

Species: `e`, `O₂⁺`, `O₄⁺`, `O⁻`, `O₂⁻`, `O₃⁻`, `O`, `O₃`, `O₂(a¹Δg)`.
The key channels (in full — `docs/PHYSICS.md`, the corrections — `docs/ERRATA.md §D`):

* impact ionization and **dissociative attachment** `e + O₂ → O⁻ + O`,
  three-body `e + 2O₂ → O₂⁻ + O₂` — the electronegativity of the gas determines the entire
  behaviour between pulses;
* the conversion `O⁻ → O₃⁻` (O₃⁻ being the dominant negative ion during the pause),
  ion-ion recombination `k_ii = 2e-13 + 2e-37·N` (the three-body contribution is ×25 of
  the two-body one at atmospheric pressure);
* field detachment **D6** and collisional detachment **D1/D7** (see the limitations: there
  an unresolved uncertainty of ×18–26 remains);
* ozone formation, secondary emission from the dielectric γ = 0.02.

### Photoprocesses (`docs/PHOTO_PROCESSES.md`)

Three spectral windows, three mechanisms — all of them switchable individually
(`--photo on|off` toggles them all at once):

| mechanism | how it is computed |
|---|---|
| **photoionization** | a three-term Helmholtz approximation (Bourdon), absorption lengths 238 / 90 / 14.8 µm at p(O₂) = 760 Torr |
| **photoemission** from the barriers | a view factor from the volumetric source to the surface |
| **photodetachment** from O₃⁻ | a transparent kernel `1/(4πR²)` — the gas does not absorb in this window |

What this yields quantitatively is in the section “Effect of the photoprocesses” below.

### Numerical scheme (`docs/NUMERICS_2D.md`)

* the (r,z) mesh is non-uniform: `betaR = 1.4` along the radius, `betaG = 1.5` across the gap,
  `dr₀ = 7.3 µm` on the axis (3–4 cells per half-width of the filament);
* the fluxes are Scharfetter–Gummel, explicit transport, `dt ≤ δ·dr²/(2D_e)`;
* the Poisson equation is a symmetrized `L_r`, a separable preconditioner + PCG,
  `strictGauss` (acceptance of a step by the Gauss residual);
* the time step is adaptive: CFL, diffusive, reactive, `dE/E`, the σ criterion;
  in the dark phase diffusion is the limiter (`dt ≈ 3.5e-11 s`), on the front it is the chemistry;
* the recording: `frames.bin` (quantized fields, log/asinh codecs) + `series.bin`
  (scalar series) + `manifest.json`; two levels of detail (`full`, `compact`).

---

## Player

![streamer front in flight](docs/shots/02-front.png)

On the left is the scene (r,z), on the right a column of plots, below a two-track timeline
(the whole run + a magnifier on the pulse, the mouse wheel changes the magnifier scale).

**The scene.** Orientation: `z` is the horizontal (metal on the left and on the right), `r` is the
vertical, mirrored about the axis. The raster is built in **physical coordinates**, not by
mesh indices: the mesh is non-uniform, and stretching by indices would show the filament
twice as wide as it is. The dielectrics are hatched, the surface charge σ(r) is
drawn as bands at both surfaces (a diverging scale, symmetric
normalization by the maximum |σ| over the whole run). Below the scene is the field profile along the
axis r = 0 with a label for the half-width `r½`.

**The plots** (`player/js/plots.mjs`, the contract is `player/js/PLOTS_API.md`):
U/I waveforms with min/max decimation (a 20 ns pulse in a 100 µs window is not lost),
the Q–V Lissajous figure with capacitance extraction and a cross-check against the analytics, the densities along
the axis, σ(r) at both surfaces, the photoprocess panel (volumetric event rate, log axis,
`S_ph/S_imp` at the cursor point) and the metric tiles.

**Synchronization.** One time cursor for everything: scrubbing the timeline moves the point on
the Lissajous figure and the waveform cursor, and a click on the waveform seeks the player.

**Export.** `PNG` — the current frame (in comparison mode, both scenes side by side) with a caption of
the run, the field and the time; `⏺ WebM` — a recording of the animation via `MediaRecorder` from the current
frame to the end of the run or to the end of the A–B segment; `CSV` — all the time series of the run
(`t, Uapp, Ugap, Icond, Idisp, Itot, Q, maxEN, o3ppm, sigmaMax, photoEmitTotalL/R`).

**Keys.** Space — play/pause, ←/→ — frame (Shift — ×10), Home/End, `[` `]` — loop
markers, `\` — reset the loop.

### What is visible in the screenshots (`docs/shots/`, the real `run-default`, full data)

| file | what is on it |
|---|---|
| `01-before.png` | **t = 4.01 µs, before breakdown.** The current is purely capacitive (I cond. 4.2e-9 A against I disp. 3.8e-6 A), max E/N = 92 Td. The electron cloud is broad, carried by drift towards the grounded barrier; the seed spot is already smeared out |
| `02-front.png` | **t = 8.2984 µs, the streamer head in flight.** The ionization is compressed into a spot `r½ = 43 µm` on the axis at R = 500 µm, the maximum at z = 0.825 mm — the head has detached from the cathode barrier and is heading for the anode. I cond. = 971 µA, max E/N = 261 Td |
| `03-ipeak.png` | **t = 8.2994 µs, the current peak (0.2096 A).** The space charge ρ: a positive (red) front ahead of and beside the channel, a negative (blue) tail on the axis towards the barrier z = 0.5 mm — the classic streamer structure. `r½ = 31.6 µm`, max E/N = 1806 Td |
| `04-afterglow.png` | **afterglow.** The accumulated trace shows the path of the head across the whole gap; the dark body of the channel behind the front — the ionization in the channel is suppressed, the field screened by the space charge |
| `05-compare.png` | **photo ON against photo OFF**, both scenes on a common scale (without a synchronous scale the per-frame autoscale equalizes the brightness and the difference disappears) |
| `06-plots-bottom.png` | the bottom panels: densities along the axis, σ(r) at both surfaces, photoprocesses (`S_ph/S_imp = 2.9e-4`) |

---

## Validation

The full table and the analyses are in `docs/VALIDATION.md`. In brief:

| # | Criterion | Verdict |
|---|---|---|
| V4a | `C_cell` (analytics/solver), 6.2587e-15 F | **PASSED** (9 digits) |
| V13 | `C_cell·U₀` | **PASSED** |
| V8a/b | σ after breakdown, the memory field | **PASSED** (−10…−15 %) |
| V9 | negative ions > n_e after the pulse | **PASSED** (×2.5e3 after 300 ns) |
| V10 | O₃ yield | **PASSED** by order of magnitude |
| V12 | charge conservation over a period | **PASSED**, 1e-14…3e-13 (after two defects were fixed) |
| V2c | the ignition `U_gap` does not grow with U₀ | **PASSED** (a spread of 4.6 % at U₀ = 6…20 kV) |
| V2a/b | the numerical value of the breakdown voltage | **DOES NOT MATCH** (+29.5 %), the cause established: formative lag |
| V7b/c | the pulse FWHM and the current density | **DOES NOT MATCH** (×15…×270): there is no external circuit, the regime is not a Townsend one |
| V4b, V5 | the capacitances from the Lissajous loop, the power by Manley’s method | **NOT MEASURED**: a full period is not simulated up to (see the limitations) |

Separately it has been checked that the Lissajous fit works: on a synthetic run with a full
period `C_cell = 6.26 fF` — **0.0 %** against the analytics, `C_diel = 69.1 fF` (+10.4 %),
R² = 0.997. On the real (truncated) runs the panel shows `C_cell` and
**refuses** to present `C_diel`, hanging a “loop not closed” badge — silently
issuing a slope over a curved loop is worse than not issuing one at all.

### The five runs in `data/`

| run | U₀ / f | photo | simulated up to | I_peak | n_e(0)/n_e(R) | r½ | filament |
|---|---|---|---|---|---|---|---|
| `run-default` | 10 kV / 10 kHz | on | 8.30 µs | 0.210 A | 1.2e6 | 31.0 µm | **yes** |
| `run-nophoto` | 10 kV / 10 kHz | off | 8.34 µs | 0.448 A | 1.0e9 | 20.9 µm | **yes** |
| `run-low` | 6 kV / 10 kHz | on | 14.72 µs | 0.277 A | 4.7e-4 | — | no, diffuse |
| `run-high` | 14 kV / 10 kHz | on | 5.76 µs | 0.273 A | 7.1e6 | 32.1 µm | **yes** |
| `run-fast` | 10 kV / 30 kHz | on | 2.60 µs | 0.198 A | 2.2e8 | 28.8 µm | **yes** |

All five break off at the **first** breakdown (6–15 % of a period) — see the limitations.

### Effect of the photoprocesses (`run-default` against `run-nophoto`)

| quantity | photo ON | OFF | difference |
|---|---|---|---|
| ignition instant | 8.2984 µs | 8.3355 µs | **−37.1 ns (−0.45 %)** |
| I_peak | 0.2096 A | 0.4475 A | **−53 %** |
| radius r½ | 31.0 µm | 20.9 µm | **+48 %** |
| n_e on the axis | 1.62e23 | 4.02e23 m⁻³ | −60 % |
| max E/N | 2461 Td | 3367 Td | −27 % |
| energy ∫U·I dt | 1.655e-7 J | 1.999e-7 J | −17 % |

**The photoprocesses have only a weak effect on the ignition instant, but a qualitative one on the morphology:**
without them the model yields a systematically over-compressed filament — twice as narrow, twice as dense,
twice as high in current. A caveat is mandatory: both snapshots were taken at the instant of the emergency
stop, so ratios and signs may be presented, but not absolute values.

---

## Limitations — an honest list

Not a single item here is “fixed by tuning”; these are the boundaries of the model's applicability.

### 1. Axisymmetry: azimuthal modes are absent in principle

The model describes **only the central channel**. An off-axis filament in the (r,z) geometry
becomes a **ring**, not a spot. Azimuthal filamentation modes are absent
not approximately but by construction: there is no ∂/∂φ in the equations.

The consequence: **the self-organization of a microdischarge lattice is not reproduced at all** — neither
the hexagonal packing of the channels, nor their mutual repulsion, nor the statistics of the number
of channels per unit area. Everything the player shows is one channel on the axis and its radial
profile. Comparison with photographs of a microdischarge lattice is not permissible.

### 2. The runs break off at the first breakdown

Not a single real run reaches the end of a period: the solver diverges in the
**near-wall** cell at the repeated ignition (`VALIDATION.md §8`). Mesh refinement
does not cure this: the structure remains one cell wide, `ρ ∝ 1/dz`, and the instant of the failure
converges to 22.3 µs. The hypothesis “the Hagelaar boundary condition will help” has been **tested and
disproved** (§13.2): in 2D it gives a failure at 8.4534 µs instead of 8.5071 µs and an `n_e`
two orders of magnitude worse — because in 1D it worked in combination with a semi-implicit
near-wall flux, whereas in 2D the near-wall faces do not enter the semi-implicit problem.

Hence: **V4b, V5 and V6 are not measured**, the Lissajous figure is not closed, and any
integral quantities over a period (the power, the transferred charge) are marked in the player with
a warning rather than shown as a number.

### 3. The cost of computation: the required 80×160 mesh is unrealistic

Measured, not estimated (`VALIDATION.md §13.1`):

| mesh | dr₀ | dt in the dark phase | hours per 1 µs |
|---|---|---|---|
| 80×160 (the specification's requirement) | 2.56 µm | 9.6e-12 s | **~1.5** → 2 periods = **12 days** |
| **32×120 (adopted)** | 7.3 µm | 3.5e-11 s | **0.23–0.29** |

The cause lies in the scheme: explicit transport, `dt ≤ δ·dr²/(2D_e)`, and refinement along the radius hits
twice. The requirement of “40 min per period” is faster than the coarsest admissible mesh by roughly
a factor of 600. This is cured only by implicit z-transport of the electrons (`NUMERICS_2D §4.5`).

### 4. D6 (field detachment): the ×18–26 uncertainty has not been removed

Two independent expert reviews diverged on which mass should be substituted into
the Wannier relation when computing `T_eff`: the full mass of the neutral or the reduced one.
The difference in the detachment frequency is **×26.3 at 30 kV/cm and ×18.2 at 37.8 kV/cm**, and it
depends on the field, that is, it does not reduce to a constant factor.

The **full mass of the neutral** has been adopted (`d6MassConvention: 'neutral'`) — on
methodological, not physical grounds: `k = 2.7e-16·√(T_eff/300)·exp(−5590/T_eff)`
is an empirical fit by Kossyi et al. 1992, and substituting someone else's definition of `T_eff` into
someone else's pre-exponential factor is not permissible regardless of which convention is the more physical. The convention
has been exposed as a parameter, and a sweep over it is in `VALIDATION.md §5б`. **This remains
an open question, not a settled one.**

### 5. The photoionization weights are carried over from air

The absorption lengths `λ_j` carry over legitimately: this is absorption **by oxygen**, it
scales with the partial pressure of O₂ (air → pure O₂ — all the lengths are shorter by a factor of 5.0).
The weights `A_j = [0.07, 0.26, 0.67]`, however, encode the shape of the **nitrogen** emission spectrum
in the 98–102.5 nm window, convolved with the O₂ absorption cross-section. The spectrum of the O₂ self-emission
in this window is different and **unknown**. The weights have been exposed as the parameter `photoIonWeights`,
and a sweep over them is mandatory for any quantitative interpretation of the photoionization.

### 6. γ = 0.02 is a circular calibration, not a confirmation

The secondary emission coefficient was tuned so as to obtain a breakdown voltage of 3.8 kV,
after which 3.8 kV was used as “validation”. This is a circle. The correct form of
presentation is the dependence `U_br(γ)` over the range 0.005…0.05, and it has been measured
(`VALIDATION.md §5а`); as an independent confirmation of the physics γ is unsuitable.

### 7. Other things worth knowing

* **LFA** (coefficients by the local `E/N`) is, strictly speaking, inapplicable on the streamer
  front — electron energy transport (LMEA) is needed. This is also the likely
  root cause of item 2.
* **There is no external circuit** (a ballast resistor, an inductance). Because of this the current pulse
  comes out ×15…×270 shorter and denser than the observed one (V7b/V7c).
* **The 1e13 m⁻³ seeding** is a numerical seeding of the first launch, not the “memory of the previous
  cycle”: in pure O₂ after 50 µs the memory is held in ion pairs, not in free
  electrons.
* **3 species out of 9 are written into the frame container** (`n_e`, `n_O3m`, `n_O3`) — the remaining six
  are honestly marked in the player with `×` “species not recorded in the run” and are not drawn.
  Collisional detachment is not in the container at all.
* **`run-synth`** in the list of runs is a synthetic container for debugging the format,
  not a solution of the equations; the player marks it separately.

---

## Structure

```
sim2d/            2D axisymmetric solver
  solver2d.mjs      time step, transport, BCs, surface charge
  physics2d.mjs     coefficients, chemistry, tables by E/N
  poisson2d.mjs     Poisson (r,z) + preconditioner + PCG
  photo2d.mjs       three photoprocesses (PHOTO_API.md — the contract)
  recorder.mjs      recording format: frames.bin / series.bin / manifest.json
  run.mjs           CLI runner
  summarize.mjs     acceptance: is there a filament, what was simulated up to
  compare-photo.mjs comparison of runs with photo on/off
  tests/            258 tests
player/           web player (ES modules, no build)
  index.html, css/player.css
  js/loader.mjs     reading the container
  js/scene.mjs      scene (r,z), σ(r), profile on the axis, afterglow, probe
  js/plots.mjs      waveforms, Lissajous, profiles, photo panel
  js/metrics.mjs    metric tiles, formatters, capacitance analytics
  js/app.mjs        state, transport, timeline, comparison, export
  js/PLOTS_API.md   the plots contract
  plots-demo.html   standalone demo page of the plots
data/             run results (+ `-compact` versions)
docs/             ERRATA (has priority), NUMERICS_2D, PHOTO_PROCESSES, PHYSICS,
                  VALIDATION, RATES_REVIEW, REFERENCE_TARGETS, UI_SPEC,
                  DIVERGENCE_ANALYSIS, shots/, README_1D.md
src/, test/, index.html   the 1D model and its UI (docs/README_1D.md)
```

The reading order of the documents: **`docs/ERRATA.md` has priority over all the others** —
it holds the blockers A1…A6, the traps of the axisymmetric solver B1…B6, the corrected validation
criteria (section C) and the corrected coefficients (section D). Where ERRATA contradicts
the other documents, ERRATA is right.
