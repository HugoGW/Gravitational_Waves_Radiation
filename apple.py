"""
Image of an apple (apple.jpeg) used as a "test mass", and the linear
(weak-field) deformation law a gravitational wave imposes on it.

Standard (linear) response of a set of test masses to a wave with h+ and
h_x polarizations, in the plane transverse to the propagation direction:

    [x']   [ 1 + h+/2      h_x/2  ] [x]
    [y'] = [   h_x/2      1 - h+/2] [y]

This is exactly the (dominant-order) geodesic-deviation deformation matrix
used by interferometers like LIGO/Virgo. No nonlinear function: just a 2x2
matrix that stretches/shears the image, applied through a matplotlib affine
transform (no manual resampling).
"""

import numpy as np
from matplotlib import image as mpimg

IMAGE_PATH = "apple.jpeg"


def load_apple_rgba(path=IMAGE_PATH, bg_threshold=225):
    """Loads apple.jpeg and makes its (near-uniform white/gray) background transparent."""
    rgb = mpimg.imread(path).astype(float) / 255.0
    is_bg = np.all(rgb > bg_threshold / 255.0, axis=-1)
    alpha = np.where(is_bg, 0.0, 1.0)
    return np.dstack([rgb, alpha])


def strain_matrix(h_plus, h_cross, factor):
    """3x3 (homogeneous, for Affine2D) matrix of the +/x deformation,
    exaggerated by a simple multiplicative factor."""
    hp = factor * h_plus
    hx = factor * h_cross
    return np.array([
        [1 + hp / 2, hx / 2, 0.0],
        [hx / 2, 1 - hp / 2, 0.0],
        [0.0, 0.0, 1.0],
    ])
