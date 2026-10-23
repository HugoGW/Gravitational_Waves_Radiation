"""
2D radiated field and h+ waveform, built from the time histories computed
in physics.py.
"""

import numpy as np

import physics as ph

# ============================================================
# Grid for the radiated field (face-on view, along the orbital angular
# momentum axis: the direction in which the l=m=2 quadrupole radiation is
# maximal, and where the two-armed pattern propagates outward).
# ============================================================
N_GRID = 300
LIM = 5.0 * ph.A0          # half-width of the displayed field, scaled with the initial separation
ARM_SIGMA = 0.25
ARM_SCALE = 75.0
R0_SOFT = 5.0
BG = 0.0

xs = np.linspace(-LIM, LIM, N_GRID)
ys = np.linspace(-LIM, LIM, N_GRID)
X, Y = np.meshgrid(xs, ys)
R = np.sqrt(X**2 + Y**2)
PHI = np.arctan2(Y, X)


def field_intensity(t_now):
    t_ret = t_now - R / ph.C
    valid = t_ret >= 0
    t_ret_local = np.mod(np.where(valid, t_ret, 0.0), ph.T_CYCLE)
    theta_ret = np.interp(t_ret_local, ph.T_HIST, ph.THETA_HIST)
    amp_ret = np.interp(t_ret_local, ph.T_HIST, ph.AMP_HIST)

    # quadrupole radiation: m=2 symmetric pattern (two arms), since
    # h_22 ~ exp(-2i*theta_orb) -> GW frequency = 2x orbital frequency.
    delta = PHI - theta_ret
    wrapped = ((delta + np.pi / 2) % np.pi) - np.pi / 2
    arm = ARM_SCALE * amp_ret / (R + R0_SOFT) * np.exp(-wrapped**2 / (2 * ARM_SIGMA**2))
    arm = np.where(valid, arm, 0.0)

    background_glow = BG * np.exp(-R / (LIM * 0.5))
    return background_glow + arm


# ============================================================
# Full h+ waveform (real physical strain, dimensionless, ~1e-21 for a
# source near GW150914 at physics.R_SOURCE_MPC) over one cycle,
# precomputed once: serves as a grayed-out reference, "completed" in
# white in real time by visualization.py. Densely sampled near merger
# (where the frequency is highest) via a non-uniform grid matched to the
# already-integrated time axis (physics.T_HIST).
#
# Newtonian-order quadrupole formula for an eccentric orbit, viewed face-on
# (zero inclination): see README for the reference.
# ============================================================
T_HIST_END = ph.T_HIST[ph.T_HIST <= ph.T_RING_END]
H_FULL_T = np.union1d(T_HIST_END, np.linspace(0.0, ph.T_RING_END, 6000))
H_FULL_THETA = np.interp(H_FULL_T, ph.T_HIST, ph.THETA_HIST)
H_FULL_H0 = np.interp(H_FULL_T, ph.T_HIST, ph.H0_HIST)
H_FULL_E = np.interp(H_FULL_T, ph.T_HIST, ph.E_HIST)
H_FULL_HPLUS = -H_FULL_H0 * (
    2 * np.cos(2 * H_FULL_THETA)
    + H_FULL_E * np.cos(H_FULL_THETA)
    + 2 * H_FULL_E * np.cos(H_FULL_THETA) ** 3
    + H_FULL_E**2
)

# h_x (cross polarization), in quadrature with h+: exact formula for a
# circular, face-on orbit (E0 = 0, physics.py's default case). If
# eccentricity is reintroduced, this expression only keeps the dominant
# harmonic (as h+ already does for its 2cos(2psi) term) and becomes an
# approximation.
H_FULL_HCROSS = -H_FULL_H0 * 2 * np.sin(2 * H_FULL_THETA)
