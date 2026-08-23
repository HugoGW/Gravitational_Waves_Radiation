# GW_emission

Animation Python (matplotlib) d'une binaire de trous noirs — masses proches
de **GW150914** — depuis l'inspiral jusqu'à la fusion et le ringdown,
rejouée en boucle. Trois panneaux synchronisés :

- **En haut** : le champ de rayonnement gravitationnel quadrupolaire, vu de
  face (le long de l'axe du moment cinétique orbital), construit par temps
  retardé — motif caractéristique à deux bras qui se resserre et s'intensifie
  à mesure que les deux trous noirs spiralent l'un vers l'autre.
- **Au milieu** : la forme d'onde `h+` réelle (déformation sans dimension,
  amplitude physique ~10⁻²¹, typique d'une détection LIGO/Virgo), tracée en
  gris sur tout le cycle et « complétée » en blanc en temps réel.
- **En bas** : l'image `apple.jpeg` déformée dynamiquement par les deux
  polarisations `h+` et `h_x` de l'onde, avec la déformation **grandement
  exagérée** (facteur affiché dans le titre du panneau) — la vraie valeur
  (~10⁻²¹) serait invisible à l'œil nu.

## Structure du projet

```
GW_emission/
├── physics.py         # masses, integration de Peters (inspiral), etat final
│                       # (spin, ringdown), historiques temporels (phase, amplitudes)
├── waveform.py         # grille du champ rayonne 2D + forme d'onde complete h+ et h_x
│                       # (consomme les historiques de physics.py)
├── apple.py             # chargement de l'image apple.jpeg + matrice de
│                       # deformation +/x (geometrie pure, independante du reste)
├── apple.jpeg           # image source de la pomme
├── visualization.py   # figure matplotlib, artistes, boucle d'animation,
│                       # affichage en direct ou export .mp4
├── main.py             # point d'entree
└── README.md
```

`physics.py` calcule tout ce qui dépend uniquement du modèle physique (aucune
dépendance à matplotlib). `waveform.py` transforme ces historiques en champ
2D et en formes d'onde `h+`/`h_x` affichables. `apple.py` charge l'image de
la pomme et fournit la matrice de déformation `+`/`×` d'un ensemble de masses
test (indépendant de la binaire — réutilisable pour n'importe quelle source
et n'importe quelle image). `visualization.py` ne contient aucune physique :
il construit la figure et anime les artistes à partir de ce que lui
fournissent les modules précédents. `main.py` importe `visualization.run()`
et lance l'animation — c'est le seul fichier à exécuter directement.

## Installation

```bash
pip install numpy matplotlib
```

Testé avec Python ≥ 3.9, numpy ≥ 1.20, matplotlib ≥ 3.5. L'export vidéo
nécessite en plus [ffmpeg](https://ffmpeg.org/) installé et accessible dans
le `PATH`.

## Utilisation

```bash
python main.py
```

Une fenêtre matplotlib s'ouvre et joue l'animation en boucle (inspiral →
fusion → ringdown → pause → nouvel inspiral). Aucun argument en ligne de
commande : tous les paramètres se règlent en modifiant les constantes en
tête de fichier (voir [Paramètres ajustables](#paramètres-ajustables)).

### Enregistrer une vidéo au lieu d'afficher l'animation

Dans `main.py`, passer `SAVE_MP4` à `True` :

```python
SAVE_MP4 = True             # False : affichage en direct (par defaut)
MP4_PATH = "GW_binary_merger.mp4"
FPS = 30
```

Dans ce mode, `visualization.run()` enregistre exactement **un cycle complet**
(inspiral → fusion → ringdown → pause) dans `MP4_PATH` via `ani.save(...,
writer="ffmpeg")`, au lieu d'ouvrir une fenêtre et de boucler indéfiniment.

## Modèle physique

Le code travaille en **unités géométrisées** `G = c = 1`, avec la masse
totale initiale du système `M = m1 + m2 = 1` comme unité de masse. Longueurs
et temps sont donc naturellement exprimés en `GM/c²` et `GM/c³` — ce sont les
unités affichées sur les axes du panneau du haut.

### 1. Inspiral — équations de Peters & Mathews

L'évolution orbite-moyennée du demi-grand axe `a` et de l'excentricité `e`
sous l'effet du rayonnement gravitationnel quadrupolaire suit :

```
da/dt = -(64/5) · (m1 m2 M) / (a³ (1-e²)^(7/2)) · (1 + 73/24 e² + 37/96 e⁴)
de/dt = -(304/15) · e (m1 m2 M) / (a⁴ (1-e²)^(5/2)) · (1 + 121/304 e²)
```

intégrées en RK4 à pas adaptatif jusqu'à ce que le périastre `a(1-e)`
atteigne un rayon de fusion schématique `R_ISCO`. L'orbite instantanée
(position, vitesse angulaire) est reconstruite à chaque pas par la solution
képlérienne exacte à `(a, e)` fixé — approximation adiabatique standard pour
ce type de modèle, valable tant que `a` et `e` varient lentement devant la
période orbitale.

On observe la **circularisation** attendue : `e` décroît beaucoup plus vite
que `a` à l'approche de la fusion (l'orbite s'aplatit avant de fusionner).

### 2. Fusion

Déclenchée quand le périastre atteint `R_ISCO`. Les deux horizons
(rayons de Schwarzschild `2·m1`, `2·m2`) sont animés en train de fusionner en
un seul horizon final sur une durée schématique `T_PLUNGE` (interpolation
cinématique lissée, pas une solution de relativité numérique).

### 3. État final et ringdown

- **Spin final** `A_F` : fit numérique pour une fusion quasi-circulaire de
  trous noirs non tournants, en fonction du rapport de masse symétrique
  `η = m1 m2 / M²`.
- **Fréquence et temps d'amortissement du ringdown** (`OMEGA_RING`,
  `TAU_RING`) : mode quasi-normal dominant `l = m = 2`, `n = 0`, obtenu par
  fit numérique en fonction du spin final seul.
- Le ringdown est modélisé comme une sinusoïde amortie
  `∝ exp(-(t-t_merge)/τ) · cos(ω_R (t-t_merge))`, raccordée en continuité de
  phase et d'amplitude à la fin de l'inspiral.

### 4. Énergie rayonnée

Calculée par bilan d'énergie orbitale newtonienne entre le début et la fin
de l'inspiral :

```
E_orb(a) = -m1 m2 / (2a)          (indépendante de e, formule vis-viva)
E_rad/M  = [E_orb(a0) - E_orb(a_merge)] / M
```

Les fréquences GW correspondantes (`F0`, `F_MERGE`) sont évaluées au
passage au périastre — la seule fréquence orbitale qui dépend à la fois de
`a` **et** de `e` :

```
f_GW,periastre(a, e) = (1/π) · sqrt(M/a³) · sqrt(1+e) / (1-e)^(3/2)
```

**Note importante** : cette méthode ne comptabilise que l'énergie perdue
*pendant l'inspiral* (jusqu'à l'ISCO). Elle sous-estime donc l'énergie
rayonnée réelle, qui inclut aussi la plongée, la fusion et le ringdown —
phases fortement relativistes où le bilan newtonien n'est plus valable. Pour
GW150914, la valeur mesurée (~5 % de la masse totale, coalescence complète)
est plus élevée que celle donnée ici.

### 5. Amplitude physique `h+`

Forme d'onde quadrupolaire newtonienne pour une orbite excentrique, vue de
face (inclinaison nulle) :

```
h+(t) = -h0(t) · [2 cos(2ψ) + e cos(ψ) + 2e cos³(ψ) + e²]
h0(t) = 2 m1 m2 / (R · p),   p = a(1-e²)
```

où `ψ` est l'anomalie vraie (phase orbitale) et `R` la distance
observateur–source, ici fixée à `R_SOURCE_MPC = 410` Mpc (distance proche de
celle de GW150914). C'est ce qui donne à `h+` son amplitude physique réaliste
(~10⁻²¹, échelle affichée en notation scientifique sur l'axe du panneau du
bas). Pendant le ringdown, `e → 0` et l'expression se réduit naturellement à
`-2·h0·cos(2ψ)`, une sinusoïde amortie pure.

### 6. Champ 2D affiché

Construit par **temps retardé** : pour chaque point de la grille à distance
`r` de l'origine, on évalue l'état du système au temps retardé
`t_ret = t - r/c`. La phase orbitale à cet instant fixe l'angle où le champ
est brillant, avec une **symétrie azimutale m = 2** (deux bras) propre au
mode quadrupolaire dominant (`h_22 ∝ exp(-2iψ)`, donc fréquence du motif
spatial = 2 × fréquence orbitale). L'intensité affichée (colorbar, "u.a.")
est une échelle arbitraire construite pour la lisibilité visuelle — voir
`ARM_SCALE`, `ARM_SIGMA`, `R0_SOFT` — distincte de la vraie amplitude
physique `h+` tracée dans le panneau du bas.

### 7. Déformation de la pomme test

Une onde gravitationnelle traversant une distribution de masses test
(perpendiculairement au plan, donc le long de l'axe déjà utilisé pour le
champ 2D et `h+`) les déplace selon la déviation géodésique au premier ordre
en champ faible, pour des polarisations `+` et `×` combinées :

```
[x']   [ 1 + h+/2      h_x/2   ] [x]
[y'] = [   h_x/2      1 - h+/2 ] [y]
```

C'est la réponse standard (linéarisée) d'un ensemble de masses test à une
onde gravitationnelle — la même physique que celle exploitée par les
interféromètres comme LIGO/Virgo, où `h+` est précisément l'allongement
relatif `ΔL/L` mesuré entre les deux bras. `h_x` est la polarisation croisée,
en quadrature de phase : `h_x(t) = -2·h0(t)·sin(2ψ)`, formule exacte pour une
orbite circulaire face-on (le cas par défaut, `E0 = 0`).

Le code charge l'image `apple.jpeg` (fond quasi blanc rendu transparent par
seuillage, voir `apple.load_apple_rgba()`), et lui applique cette matrice de
déformation à chaque frame via une transformation affine matplotlib
(`Affine2D`, pas de réechantillonnage manuel ni de fonction non linéaire).
`h+` et `h_x` sont pris directement dans la forme d'onde déjà calculée
(`waveform.H_FULL_HPLUS`, `waveform.H_FULL_HCROSS`), simplement multipliés
par un facteur constant `APPLE_EXAGGERATION = 2×10²⁰` avant d'entrer dans la
matrice — car la vraie valeur (~10⁻²¹) ne produirait aucun mouvement
visible. Ce facteur est affiché explicitement dans le titre du panneau.

## Paramètres ajustables

| Constante | Module | Rôle |
|---|---|---|
| `M1_MSUN`, `M2_MSUN` | `physics.py` | Masses des deux trous noirs (M☉) |
| `A0`, `E0` | `physics.py` | Séparation et excentricité initiales (unités de M) |
| `R_ISCO` | `physics.py` | Rayon de fusion schématique (unités de M) |
| `R_SOURCE_MPC` | `physics.py` | Distance de la source, pilote l'amplitude physique de `h+` |
| `T_PLUNGE` | `physics.py` | Durée schématique de la coalescence |
| `N_TAU_RINGDOWN` | `physics.py` | Durée du ringdown affichée, en multiples de `τ` |
| `GAP` | `physics.py` | Pause silencieuse avant que la boucle ne reprenne |
| `LIM`, `N_GRID` | `waveform.py` | Étendue et résolution de la grille du champ 2D |
| `ARM_SIGMA`, `ARM_SCALE`, `R0_SOFT`, `BG` | `waveform.py` | Réglages purement visuels du champ 2D (u.a.) |
| `FRAME_DT_SIM` | `visualization.py` | Vitesse de lecture (unités de M avancées par frame) |
| `APPLE_EXAGGERATION` | `visualization.py` | Facteur multipliant `h+`/`h_x` pour la déformation de la pomme (affiché dans le titre) |
| `SAVE_MP4`, `MP4_PATH`, `FPS` | `main.py` | Affichage en direct ou export vidéo |

## Références

| Élément du code | Référence |
|---|---|
| `peters_rates()` — évolution orbite-moyennée `da/dt`, `de/dt` | Peters, P. C. (1964), *Gravitational Radiation and the Motion of Two Point Masses*, Phys. Rev. **136**, B1224 |
| Formule de puissance quadrupolaire associée | Peters, P. C. & Mathews, J. (1963), *Gravitational Radiation from Point Masses in a Keplerian Orbit*, Phys. Rev. **131**, 435 |
| Amplitude d'onde, énergie orbitale newtonienne, formule `h+` excentrique | Maggiore, M. (2008), *Gravitational Waves, Vol. 1: Theory and Experiments*, Oxford University Press, chap. 4–5 |
| `A_F` — spin final du trou noir issu de la fusion | Rezzolla, L., Barausse, E., Dorband, E. N., Pollney, D., Reisswig, C., Seiler, J. & Husa, S. (2008), *On the final spin from the coalescence of two black holes*, Phys. Rev. D **78**, 044002 |
| `OMEGA_RING`, `TAU_RING` — mode quasi-normal dominant (l=m=2, n=0) | Berti, E., Cardoso, V. & Will, C. M. (2006), *On gravitational-wave spectroscopy of massive black holes with the space interferometer LISA*, Phys. Rev. D **73**, 064030 |
| Masses `M1_MSUN`, `M2_MSUN`, distance `R_SOURCE_MPC` | Abbott, B. P. et al. (LIGO Scientific Collaboration & Virgo Collaboration) (2016), *Observation of Gravitational Waves from a Binary Black Hole Merger*, Phys. Rev. Lett. **116**, 061102 |
| Rayon ISCO / rayon de Schwarzschild | Résultat standard de relativité générale — voir p. ex. Misner, C. W., Thorne, K. S. & Wheeler, J. A. (1973), *Gravitation*, W. H. Freeman |

## Limites et approximations (assumées volontairement)

- L'inspiral utilise les équations de Peters à l'ordre newtonien
  orbite-moyenné : correct loin de la fusion, mais devient approximatif
  près de l'ISCO où les effets relativistes forts dominent.
- L'énergie rayonnée affichée ne comptabilise que l'inspiral (voir
  [§4](#4-énergie-rayonnée)), pas la coalescence complète.
- La position pendant la plongée finale et le raccord inspiral → ringdown
  sont une interpolation cinématique schématique, pas une solution de
  relativité numérique.
- Le champ 2D affiché est une construction pédagogique par temps retardé
  reproduisant la structure en spirale à deux bras du mode quadrupolaire
  dominant vu de face ; ce n'est pas une résolution des équations d'Einstein
  linéarisées sur une grille, et son intensité (colorbar) est en unités
  arbitraires — seul le panneau `h+` porte une amplitude physique réelle.
- Orbite non précessante, spins initiaux nuls, observateur face-on
  (inclinaison nulle) : pas de modulation par inclinaison ni d'effets de
  précession de spin.
- La déformation de la pomme applique la loi linéarisée du champ faible
  (`h+` et `h_x` combinés) uniformément à toute l'image (approximation de
  dipôle court : la pomme est supposée bien plus petite que la longueur
  d'onde du signal, ce qui est vrai ici par plusieurs ordres de grandeur).
  `h_x` utilise la formule exacte pour une orbite circulaire face-on
  (`E0 = 0`, cas par défaut) ; si l'excentricité est réactivée, elle ne
  garde que l'harmonique dominante, comme `h+` le fait déjà pour son propre
  terme en `2cos(2ψ)`. L'amplitude est multipliée par un facteur constant
  `APPLE_EXAGGERATION` (~10²⁰, affiché dans le titre du panneau) pour être
  visible — sans quoi aucune déformation ne serait perceptible.
