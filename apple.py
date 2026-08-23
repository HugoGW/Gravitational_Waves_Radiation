"""
Image d'une pomme (apple.jpeg) utilisee comme "masse test", et la loi de
deformation lineaire (champ faible) qu'une onde gravitationnelle lui impose.

Reponse standard (lineaire) d'un ensemble de masses test a une onde de
polarisations h+ et h_x, dans le plan transverse a la propagation :

    [x']   [ 1 + h+/2      h_x/2  ] [x]
    [y'] = [   h_x/2      1 - h+/2] [y]

C'est exactement la matrice de deformation (deviation geodesique, ordre
dominant) utilisee par les interferometres comme LIGO/Virgo. Aucune fonction
non lineaire : juste une matrice 2x2 qui etire/cisaille l'image, appliquee
via une transformation affine matplotlib (pas de reechantillonnage manuel).
"""

import numpy as np
from matplotlib import image as mpimg

IMAGE_PATH = "apple.jpeg"


def load_apple_rgba(path=IMAGE_PATH, bg_threshold=225):
    """Charge apple.jpeg et rend transparent son fond (quasi blanc/gris uni)."""
    rgb = mpimg.imread(path).astype(float) / 255.0
    is_bg = np.all(rgb > bg_threshold / 255.0, axis=-1)
    alpha = np.where(is_bg, 0.0, 1.0)
    return np.dstack([rgb, alpha])


def strain_matrix(h_plus, h_cross, factor):
    """Matrice 3x3 (homogene, pour Affine2D) de la deformation + / x,
    exageree d'un simple facteur multiplicatif."""
    hp = factor * h_plus
    hx = factor * h_cross
    return np.array([
        [1 + hp / 2, hx / 2, 0.0],
        [hx / 2, 1 - hp / 2, 0.0],
        [0.0, 0.0, 1.0],
    ])
