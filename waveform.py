"""
Champ rayonne 2D et forme d'onde h+, construits a partir des historiques
temporels calcules dans physics.py.
"""

import numpy as np

import physics as ph

# ============================================================
# Grille du champ rayonne (vue de face, le long de l'axe du moment
# cinetique orbital : c'est la direction ou le rayonnement quadrupolaire
# l=m=2 est maximal, et ou le motif a 2 bras se propage vers l'exterieur).
# ============================================================
N_GRID = 300
LIM = 200.0
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

    # rayonnement quadrupolaire : motif a symetrie m=2 (deux bras),
    # car h_22 ~ exp(-2i*theta_orb) -> frequence GW = 2x frequence orbitale.
    delta = PHI - theta_ret
    wrapped = ((delta + np.pi / 2) % np.pi) - np.pi / 2
    arm = ARM_SCALE * amp_ret / (R + R0_SOFT) * np.exp(-wrapped**2 / (2 * ARM_SIGMA**2))
    arm = np.where(valid, arm, 0.0)

    background_glow = BG * np.exp(-R / (LIM * 0.5))
    return background_glow + arm


# ============================================================
# Forme d'onde complete h+ (deformation physique reelle, sans dimension,
# ~1e-21 pour une source proche de GW150914 a physics.R_SOURCE_MPC) sur un
# cycle, precalculee une fois : sert de reference grisee, "completee" en
# blanc en temps reel par visualization.py. Echantillonnage dense pres de la
# fusion (ou la frequence est la plus elevee) via une grille non-uniforme
# calquee sur l'axe des temps deja integre (physics.T_HIST).
#
# Formule quadrupolaire (ordre newtonien dominant) pour une orbite
# excentrique, vue de face (inclinaison nulle) : cf. README pour la
# reference.
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

# h_x (polarisation croisee), quadrature de h+ : formule exacte pour une
# orbite circulaire face-on (E0 = 0, cas par defaut de physics.py). Si une
# excentricite est reintroduite, cette expression ne garde que l'harmonique
# dominante (comme h+ le fait deja pour son terme en 2cos(2psi)) et devient
# une approximation.
H_FULL_HCROSS = -H_FULL_H0 * 2 * np.sin(2 * H_FULL_THETA)
