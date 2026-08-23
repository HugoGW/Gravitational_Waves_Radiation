"""
Modele physique de la binaire de trous noirs.

Masses, integration de l'inspiral (equations de Peters), etat final
(spin, mode de ringdown) et historiques temporels (phase orbitale,
amplitudes) reutilises par waveform.py pour construire le champ rayonne
affiche et la forme d'onde h+. Voir README.md pour les references.
"""

import numpy as np

# ============================================================
# Masses (proches de GW150914) et unites geometrisees G = c = 1,
# masse totale M = m1 + m2 = 1.
# ============================================================
M1_MSUN, M2_MSUN = 36.0, 29.0
M_TOTAL_MSUN = M1_MSUN + M2_MSUN
m1 = M1_MSUN / M_TOTAL_MSUN
m2 = M2_MSUN / M_TOTAL_MSUN
M = m1 + m2  # == 1
ETA = m1 * m2 / M**2  # rapport de masse symetrique

A0, E0 = 22.0, 0.0      # separation et excentricite initiales (en unites de M)
R_ISCO = 3.0 * M         # rayon de fusion approx. (ISCO de Schwarzschild)
C = 1.0                  # vitesse des ondes = vitesse de la lumiere (G=c=1)

# --- Constantes physiques (SI), pour convertir h+ en amplitude de deformation
# reelle (sans dimension, ~1e-21 pour une binaire de type GW150914) ---
G_SI = 6.67430e-11
C_SI = 299792458.0
MSUN_KG = 1.98892e30
MPC_TO_M = 3.0857e22
R_SOURCE_MPC = 410.0     # distance de la source (proche de la distance de GW150914)

M_TOTAL_KG = M_TOTAL_MSUN * MSUN_KG
LENGTH_UNIT_M = G_SI * M_TOTAL_KG / C_SI**2          # 1 unite de longueur interne = GM/c^2, en metres
R_GEOM = (R_SOURCE_MPC * MPC_TO_M) / LENGTH_UNIT_M   # distance de la source, en unites de M


# ============================================================
# Inspiral : equations de Peters (orbite-moyennees)
# ============================================================
def peters_rates(a, e):
    one_me2 = 1 - e**2
    da_dt = -64.0 / 5.0 * m1 * m2 * M / (a**3 * one_me2**3.5) * (1 + 73.0 / 24 * e**2 + 37.0 / 96 * e**4)
    de_dt = -304.0 / 15.0 * e * m1 * m2 * M / (a**4 * one_me2**2.5) * (1 + 121.0 / 304 * e**2)
    return da_dt, de_dt


def theta_dot(a, e, theta):
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    h = np.sqrt(M * a * (1 - e**2))  # moment cinetique specifique (orbite relative)
    return h / r**2, r


def _derivs(state):
    a, e, theta = state
    da_dt, de_dt = peters_rates(a, e)
    dtheta_dt, _ = theta_dot(a, e, theta)
    return np.array([da_dt, de_dt, dtheta_dt])


def _rk4_step(state, dt):
    k1 = _derivs(state)
    k2 = _derivs(state + 0.5 * dt * k1)
    k3 = _derivs(state + 0.5 * dt * k2)
    k4 = _derivs(state + dt * k3)
    return state + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate_inspiral(a0, e0):
    state = np.array([a0, e0, 0.0])
    t = 0.0
    ts, aas, ees, thetas = [0.0], [a0], [e0], [0.0]
    while True:
        a, e, theta = state
        if a * (1 - e) <= R_ISCO:
            break
        da_dt, de_dt = peters_rates(a, e)
        dt_a = 0.01 * a / max(abs(da_dt), 1e-30)
        dt_theta = 0.02 * (2 * np.pi) / theta_dot(a, e, theta)[0]
        dt = min(dt_a, dt_theta, 50.0)
        state = _rk4_step(state, dt)
        state[1] = max(state[1], 0.0)
        t += dt
        ts.append(t); aas.append(state[0]); ees.append(state[1]); thetas.append(state[2])
        if len(ts) > 400000:
            break
    return np.array(ts), np.array(aas), np.array(ees), np.array(thetas)


ts, aas, ees, thetas = integrate_inspiral(A0, E0)
T_MERGE = ts[-1]
A_MERGE, E_MERGE = aas[-1], ees[-1]
THETA_MERGE = thetas[-1]


# ============================================================
# Energie rayonnee = E_orb(f_0) - E_orb(f_final), bilan d'energie newtonien
# ============================================================
def orbital_energy(a):
    # energie orbitale du systeme relatif : -G m1 m2 / (2a), independante de e (vis-viva)
    return -0.5 * m1 * m2 / a


def gw_freq_periastron(a, e):
    # frequence GW (mode dominant) au passage au periastre : depend de a ET e
    n = np.sqrt(M / a**3)
    return (n / np.pi) * np.sqrt(1 + e) / (1 - e) ** 1.5


F0 = gw_freq_periastron(A0, E0)                  # frequence GW initiale, deduite de a0 et e0
F_MERGE = gw_freq_periastron(A_MERGE, E_MERGE)   # frequence GW a la fin de l'inspiral

ERAD_FRAC = (orbital_energy(A0) - orbital_energy(A_MERGE)) / M  # E(f_0) - E(f_final), fraction de M

# --- Etat final (fit NR usuel pour une fusion quasi-circulaire non tournante) ---
A_F = 2 * np.sqrt(3) * ETA - 3.5171 * ETA**2 + 2.5763 * ETA**3     # spin final (Rezzolla et al. 2008)
M_F = M * (1 - ERAD_FRAC)

# --- Mode quasi-normal dominant l=m=2, n=0 (fit de Berti, Cardoso & Will 2006) ---
_F1, _F2, _F3 = 1.5251, -1.1568, 0.1292
_Q1, _Q2, _Q3 = 0.7000, 1.4187, -0.4990
_MF_OMEGA = _F1 + _F2 * (1 - A_F) ** _F3
_QF = _Q1 + _Q2 * (1 - A_F) ** _Q3
OMEGA_RING = _MF_OMEGA / M_F       # frequence angulaire du mode (2,2), en 1/M
TAU_RING = 2 * _QF / OMEGA_RING     # temps d'amortissement, en M

# --- Rayons de Schwarzschild (pour l'affichage des horizons) ---
RS1, RS2, RS_F = 2 * m1, 2 * m2, 2 * M_F

T_PLUNGE = 25.0                 # duree schematique de la coalescence (en M)
N_TAU_RINGDOWN = 9.0            # duree du ringdown affichee, en multiples de tau
GAP = 600.0                     # pause "silencieuse" avant que la boucle ne reprenne :
                                 # le ringdown est deja quasi eteint apres N_TAU_RINGDOWN*TAU_RING,
                                 # cette pause laisse en plus le reste du champ sortir du cadre
                                 # (silencieusement, amp=0) avant qu'un nouvel inspiral ne commence
T_RING_END = T_MERGE + T_PLUNGE + N_TAU_RINGDOWN * TAU_RING
T_CYCLE = T_RING_END + GAP

# ============================================================
# Historiques temporels (phase, amplitudes) sur inspiral + ringdown + pause,
# utilises par waveform.py (temps retardu, forme d'onde h+)
# ============================================================
# --- Amplitude "visuelle" (formule quadrupolaire, ordre dominant, sans dimension) ---
n_orb = np.sqrt(M / aas**3)
amp_insp = 4 * ETA * M * (M * n_orb) ** (2.0 / 3.0)
AMP_MERGE = amp_insp[-1]

t_ring = np.linspace(T_MERGE, T_RING_END, 4000)
theta_ring = THETA_MERGE + (OMEGA_RING / 2.0) * (t_ring - T_MERGE)
amp_ring = AMP_MERGE * np.exp(-(t_ring - T_MERGE) / TAU_RING)

t_gap = np.array([T_RING_END, T_CYCLE])
theta_gap = np.array([theta_ring[-1], theta_ring[-1]])
amp_gap = np.array([0.0, 0.0])

T_HIST = np.concatenate([ts, t_ring[1:], t_gap])
THETA_HIST = np.concatenate([thetas, theta_ring[1:], theta_gap])
AMP_HIST = np.concatenate([amp_insp, amp_ring[1:], amp_gap])

# --- Amplitude physique h0 = 2 m1 m2 / (R * p) (p = a(1-e^2), forme newtonienne
# du quadrupole en unites geometrisees) et excentricite, pour la vraie
# deformation h+ (~1e-21) ---
p_insp = aas * (1 - ees**2)
h0_insp = 2 * m1 * m2 / (R_GEOM * p_insp)
H0_MERGE = h0_insp[-1]

h0_ring = H0_MERGE * np.exp(-(t_ring - T_MERGE) / TAU_RING)
h0_gap = np.array([0.0, 0.0])
e_ring = np.zeros_like(t_ring)      # deja quasi-circularise a la fusion (cf. Peters-Mathews)
e_gap = np.array([0.0, 0.0])

H0_HIST = np.concatenate([h0_insp, h0_ring[1:], h0_gap])
E_HIST = np.concatenate([ees, e_ring[1:], e_gap])


# ============================================================
# Positions/tailles instantanees des trous noirs (non retarde), pour l'affichage
# ============================================================
def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x**2 * (3 - 2 * x)


def bodies_now(t_local):
    """Retourne les positions/tailles des trous noirs a l'instant courant (non retarde)."""
    if t_local <= T_MERGE:
        a = np.interp(t_local, ts, aas)
        e = np.interp(t_local, ts, ees)
        theta = np.interp(t_local, ts, thetas)
        r = a * (1 - e**2) / (1 + e * np.cos(theta))
        rel = np.array([r * np.cos(theta), r * np.sin(theta)])
        pos1 = (m2 / M) * rel
        pos2 = -(m1 / M) * rel
        return dict(binary=True, pos1=pos1, pos2=pos2, rs1=RS1, rs2=RS2,
                    final_alpha=0.0, final_pos=np.zeros(2), final_rs=0.0)

    a_m = A_MERGE
    e_m = E_MERGE
    r_m = a_m * (1 - e_m**2) / (1 + e_m * np.cos(THETA_MERGE))
    rel_m = np.array([r_m * np.cos(THETA_MERGE), r_m * np.sin(THETA_MERGE)])
    pos1_m = (m2 / M) * rel_m
    pos2_m = -(m1 / M) * rel_m

    merge_frac = np.clip((t_local - T_MERGE) / T_PLUNGE, 0.0, 1.0)
    shrink = 1 - _smoothstep(merge_frac)
    return dict(binary=shrink > 1e-3,
                pos1=pos1_m * shrink, pos2=pos2_m * shrink,
                rs1=RS1 * shrink, rs2=RS2 * shrink,
                final_alpha=_smoothstep(merge_frac),
                final_pos=np.zeros(2), final_rs=RS_F * _smoothstep(merge_frac))
