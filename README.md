# GW_emission

Python (matplotlib) animation of a black-hole binary — masses close to
**GW150914** — from inspiral through merger and ringdown, played on a loop.
Three synchronized panels:

- **Top**: the quadrupolar gravitational-wave radiated field, viewed
  face-on (along the orbital angular momentum axis), built via retarded
  time — the characteristic two-armed pattern tightens and brightens as the
  two black holes spiral toward each other.
- **Middle**: the real `h+` waveform (dimensionless strain, physical
  amplitude ~10⁻²¹, typical of a LIGO/Virgo detection), drawn in gray over
  the whole cycle and "completed" in white in real time. The time axis is
  converted from `GM/c³` to seconds, minutes, hours or days (the most
  readable unit is picked automatically from the signal duration).
- **Bottom**: side by side, the `apple.jpeg` image and a drawing of the
  LISA constellation (three spacecraft linked by laser arms, with the
  undeformed triangle shown dashed), both dynamically deformed by the
  wave's two polarizations `h+` and `h_x`, with the deformation **greatly
  exaggerated** (the factor is shown above the panels) — the real value
  (~10⁻²¹) would be invisible to the naked eye.

## Project structure

```
GW_emission/
├── physics.py         # masses, Peters-equation integration (inspiral), final state
│                       # (spin, ringdown), time histories (phase, amplitudes)
├── waveform.py         # 2D radiated field grid + full h+ and h_x waveforms
│                       # (consumes the histories from physics.py)
├── apple.py             # loads apple.jpeg + the +/x deformation matrix
│                       # (pure geometry, independent of everything else)
├── apple.jpeg           # source image of the apple
├── lisa.py              # matplotlib drawing of the LISA constellation,
│                       # deformed by the same +/x matrix as the apple
├── visualization.py   # matplotlib figure, artists, animation loop,
│                       # live display or .mp4 export
├── main.py             # entry point
└── README.md
```

`physics.py` computes everything that depends only on the physical model
(no dependency on matplotlib). `waveform.py` turns these histories into a
displayable 2D field and `h+`/`h_x` waveforms. `apple.py` loads the apple
image and provides the `+`/`×` deformation matrix for a set of test masses
(independent of the binary — reusable for any source and any image).
`visualization.py` contains no physics: it builds the figure and animates
the artists from what the previous modules provide. `main.py` imports
`visualization.run()` and launches the animation — it's the only file meant
to be run directly.

## Installation

```bash
pip install numpy matplotlib
```

Tested with Python ≥ 3.9, numpy ≥ 1.20, matplotlib ≥ 3.5. Video export
additionally requires [ffmpeg](https://ffmpeg.org/) installed and available
on the `PATH`.

## Usage

```bash
python main.py
```

A matplotlib window opens and plays the animation on a loop (inspiral →
merger → ringdown → pause → new inspiral). No command-line arguments: all
parameters are set by editing the constants at the top of each file (see
[Adjustable parameters](#adjustable-parameters)).

### Saving a video instead of displaying the animation

In `main.py`, set `SAVE_MP4` to `True`:

```python
SAVE_MP4 = True             # False: live display (default)
MP4_PATH = "GW_binary_merger.mp4"
FPS = 30
```

In this mode, `visualization.run()` saves exactly **one full cycle**
(inspiral → merger → ringdown → pause) to `MP4_PATH` via `ani.save(...,
writer="ffmpeg")`, instead of opening a window and looping indefinitely.

## Physical model

The code works in **geometrized units** `G = c = 1`, with the system's
initial total mass `M = m1 + m2 = 1` as the mass unit. Lengths and times are
therefore naturally expressed in `GM/c²` and `GM/c³`. `GM/c²` is the unit
shown on the top panel's axes; times on the `h+` panel are converted to
physical units via `TIME_UNIT_S = GM/c³` (in seconds, see `physics.py`).

### 1. Inspiral — Peters & Mathews equations

The orbit-averaged evolution of the semi-major axis `a` and eccentricity
`e` under quadrupolar gravitational-wave emission follows:

```
da/dt = -(64/5) · (m1 m2 M) / (a³ (1-e²)^(7/2)) · (1 + 73/24 e² + 37/96 e⁴)
de/dt = -(304/15) · e (m1 m2 M) / (a⁴ (1-e²)^(5/2)) · (1 + 121/304 e²)
```

integrated with adaptive-step RK4 until the periastron `a(1-e)` reaches a
schematic merger radius `R_ISCO`. The instantaneous orbit (position,
angular velocity) is reconstructed at each step from the exact Keplerian
solution at fixed `(a, e)` — a standard adiabatic approximation for this
kind of model, valid as long as `a` and `e` vary slowly compared to the
orbital period.

The expected **circularization** is visible: `e` decays much faster than
`a` as the merger approaches (the orbit flattens out before merging).

### 2. Merger

Triggered when the periastron reaches `R_ISCO`. The two horizons
(Schwarzschild radii `2·m1`, `2·m2`) are animated merging into a single
final horizon over a schematic duration `T_PLUNGE` (smoothed kinematic
interpolation, not a numerical-relativity solution).

### 3. Final state and ringdown

- **Final spin** `A_F`: numerical fit for a quasi-circular merger of
  non-spinning black holes, as a function of the symmetric mass ratio
  `η = m1 m2 / M²`.
- **Ringdown frequency and damping time** (`OMEGA_RING`, `TAU_RING`):
  dominant quasi-normal mode `l = m = 2`, `n = 0`, obtained from a
  numerical fit depending only on the final spin.
- The ringdown is modeled as a damped sinusoid
  `∝ exp(-(t-t_merge)/τ) · cos(ω_R (t-t_merge))`, matched in phase and
  amplitude continuity to the end of the inspiral.

### 4. Radiated energy

Computed from a Newtonian orbital-energy balance between the start and end
of the inspiral:

```
E_orb(a) = -m1 m2 / (2a)          (independent of e, vis-viva formula)
E_rad/M  = [E_orb(a0) - E_orb(a_merge)] / M
```

The corresponding GW frequencies (`F0`, `F_MERGE`) are evaluated at
periastron passage — the only orbital frequency that depends on both `a`
**and** `e`:

```
f_GW,periastron(a, e) = (1/π) · sqrt(M/a³) · sqrt(1+e) / (1-e)^(3/2)
```

**Important note**: this method only accounts for the energy lost *during
the inspiral* (up to the ISCO). It therefore underestimates the true
radiated energy, which also includes the plunge, merger and ringdown —
strongly relativistic phases where the Newtonian balance no longer holds.
For GW150914, the measured value (~5% of the total mass, full coalescence)
is higher than the one given here.

### 5. Physical `h+` amplitude

Newtonian-order quadrupole waveform for an eccentric orbit, viewed face-on
(zero inclination):

```
h+(t) = -h0(t) · [2 cos(2ψ) + e cos(ψ) + 2e cos³(ψ) + e²]
h0(t) = 2 m1 m2 / (R · p),   p = a(1-e²)
```

where `ψ` is the true anomaly (orbital phase) and `R` the observer–source
distance, here fixed to `R_SOURCE_MPC = 410` Mpc (close to GW150914's
distance). This is what gives `h+` its realistic physical amplitude
(~10⁻²¹, shown in scientific notation on the bottom axis of the middle
panel). During the ringdown, `e → 0` and the expression naturally reduces
to `-2·h0·cos(2ψ)`, a pure damped sinusoid.

### 6. Displayed 2D field

Built via **retarded time**: for each point of the grid at distance `r`
from the origin, the state of the system is evaluated at the retarded time
`t_ret = t - r/c`. The orbital phase at that instant sets the angle at
which the field is bright, with an **m = 2 azimuthal symmetry** (two arms)
characteristic of the dominant quadrupole mode (`h_22 ∝ exp(-2iψ)`, hence
spatial-pattern frequency = 2 × orbital frequency). The displayed intensity
(colorbar, "a.u.") is an arbitrary scale built for visual clarity — see
`ARM_SCALE`, `ARM_SIGMA`, `R0_SOFT` — distinct from the real physical `h+`
amplitude plotted in the middle panel.

### 7. Deformation of the test apple

A gravitational wave crossing a distribution of test masses
(perpendicular to the plane, i.e. along the axis already used for the 2D
field and `h+`) displaces them according to first-order weak-field
geodesic deviation, for combined `+` and `×` polarizations:

```
[x']   [ 1 + h+/2      h_x/2   ] [x]
[y'] = [   h_x/2      1 - h+/2 ] [y]
```

This is the standard (linearized) response of a set of test masses to a
gravitational wave — the same physics exploited by interferometers like
LIGO/Virgo, where `h+` is precisely the relative elongation `ΔL/L` measured
between the two arms. `h_x` is the cross polarization, in phase quadrature:
`h_x(t) = -2·h0(t)·sin(2ψ)`, an exact formula for a circular, face-on orbit
(the default case, `E0 = 0`).

The code loads the `apple.jpeg` image (near-white background made
transparent by thresholding, see `apple.load_apple_rgba()`), and applies
this deformation matrix to it every frame via a matplotlib affine
transform (`Affine2D`, no manual resampling and no nonlinear function).
`h+` and `h_x` are taken directly from the already-computed waveform
(`waveform.H_FULL_HPLUS`, `waveform.H_FULL_HCROSS`), simply multiplied by a
constant factor `APPLE_EXAGGERATION = APPLE_MAX_STRAIN / max|h|` (with
`APPLE_MAX_STRAIN = 0.5`, i.e. ~10²⁰ for GW150914) before entering the matrix —
because the real value (~10⁻²¹) would produce no visible motion at all.
This factor is displayed explicitly in the panel's title.

## Adjustable parameters

| Constant | Module | Role |
|---|---|---|
| `M1_MSUN`, `M2_MSUN` | `physics.py` | Masses of the two black holes (M☉) |
| `A0`, `E0` | `physics.py` | Initial separation and eccentricity (units of M) |
| `R_ISCO` | `physics.py` | Schematic merger radius (units of M) |
| `R_SOURCE_MPC` | `physics.py` | Source distance, drives `h+`'s physical amplitude |
| `T_PLUNGE` | `physics.py` | Schematic coalescence duration |
| `N_TAU_RINGDOWN` | `physics.py` | Ringdown duration shown, in multiples of `τ` |
| `GAP` | `physics.py` | Silent pause before the loop restarts |
| `LIM`, `N_GRID` | `waveform.py` | Extent (`LIM = 5 * A0`) and resolution of the 2D field grid |
| `ARM_SIGMA`, `ARM_SCALE`, `R0_SOFT`, `BG` | `waveform.py` | Purely visual settings for the 2D field (a.u.) |
| `FRAME_DT_SIM` | `visualization.py` | Playback speed (M units advanced per frame) |
| `APPLE_MAX_STRAIN` | `visualization.py` | Peak exaggerated strain of the apple/LISA drawings; sets `APPLE_EXAGGERATION = APPLE_MAX_STRAIN / max\|h\|` (shown in the title) |
| `APPLE_HALF_EXTENT`, `APPLE_Y_OFFSET` | `visualization.py` | Apple image size and vertical position in its panel |
| `SAVE_MP4`, `MP4_PATH`, `FPS` | `main.py` | Live display or video export |

## References

| Code element | Reference |
|---|---|
| `peters_rates()` — orbit-averaged evolution `da/dt`, `de/dt` | Peters, P. C. (1964), *Gravitational Radiation and the Motion of Two Point Masses*, Phys. Rev. **136**, B1224 |
| Associated quadrupole power formula | Peters, P. C. & Mathews, J. (1963), *Gravitational Radiation from Point Masses in a Keplerian Orbit*, Phys. Rev. **131**, 435 |
| Wave amplitude, Newtonian orbital energy, eccentric `h+` formula | Maggiore, M. (2008), *Gravitational Waves, Vol. 1: Theory and Experiments*, Oxford University Press, chap. 4–5 |
| `A_F` — final spin of the merged black hole | Rezzolla, L., Barausse, E., Dorband, E. N., Pollney, D., Reisswig, C., Seiler, J. & Husa, S. (2008), *On the final spin from the coalescence of two black holes*, Phys. Rev. D **78**, 044002 |
| `OMEGA_RING`, `TAU_RING` — dominant quasi-normal mode (l=m=2, n=0) | Berti, E., Cardoso, V. & Will, C. M. (2006), *On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA*, Phys. Rev. D **73**, 064030 |
| Masses `M1_MSUN`, `M2_MSUN`, distance `R_SOURCE_MPC` | Abbott, B. P. et al. (LIGO Scientific Collaboration & Virgo Collaboration) (2016), *Observation of Gravitational Waves from a Binary Black Hole Merger*, Phys. Rev. Lett. **116**, 061102 |
| ISCO radius / Schwarzschild radius | Standard general-relativity result — see e.g. Misner, C. W., Thorne, K. S. & Wheeler, J. A. (1973), *Gravitation*, W. H. Freeman |

## Limitations and approximations (deliberately assumed)

- The inspiral uses Peters' equations at orbit-averaged Newtonian order:
  correct far from merger, but becomes approximate near the ISCO where
  strong relativistic effects dominate.
- The displayed radiated energy only accounts for the inspiral (see
  [§4](#4-radiated-energy)), not the full coalescence.
- The position during the final plunge and the inspiral → ringdown
  matching are a schematic kinematic interpolation, not a
  numerical-relativity solution.
- The displayed 2D field is a pedagogical construction via retarded time
  reproducing the two-armed spiral structure of the dominant quadrupole
  mode viewed face-on; it is not a solution of the linearized Einstein
  equations on a grid, and its intensity (colorbar) is in arbitrary
  units — only the `h+` panel carries a real physical amplitude.
- Non-precessing orbit, zero initial spins, face-on observer (zero
  inclination): no modulation from inclination or spin-precession effects.
- The apple's deformation applies the linearized weak-field law (`h+` and
  `h_x` combined) uniformly across the whole image (short-dipole
  approximation: the apple is assumed much smaller than the signal's
  wavelength, which holds here by several orders of magnitude). `h_x` uses
  the exact formula for a circular, face-on orbit (`E0 = 0`, the default
  case); if eccentricity is reintroduced, it only keeps the dominant
  harmonic, as `h+` already does for its own `2cos(2ψ)` term. The amplitude
  is multiplied by a constant factor `APPLE_EXAGGERATION` (~10²⁰, shown in
  the panel's title) to be visible — without it, no deformation would be
  perceptible.
