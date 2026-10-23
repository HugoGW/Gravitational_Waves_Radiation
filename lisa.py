"""
Schematic drawing of the LISA constellation (three spacecraft at the
vertices of an equilateral triangle, linked by laser arms), used as a
second set of "test masses" next to the apple.

The whole drawing is deformed by the same linear +/x strain matrix as the
apple (apple.strain_matrix), applied through a matplotlib affine transform:
the arm lengths change exactly as an interferometer would measure them.
A dashed, undeformed triangle is kept as a reference.
"""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon, Rectangle, RegularPolygon
from matplotlib.transforms import Affine2D

CIRCUMRADIUS = 1.15        # distance from the constellation's center to each spacecraft
SC_RADIUS = 0.13           # spacecraft body radius (drawing units)
ANGLES = np.deg2rad([90.0, 210.0, 330.0])
ARM_LENGTH_KM = 2.5e6      # real LISA arm length (shown in the panel)

LASER_COLOR = "#ff4040"
BODY_COLOR = "#c9ced6"
TRIM_COLOR = "#e8b84a"     # gold (thermal blanket)
PANEL_COLOR = "#2f5fb3"    # solar cells
TELESCOPE_COLOR = "#8a9099"


def spacecraft_positions():
    return CIRCUMRADIUS * np.column_stack([np.cos(ANGLES), np.sin(ANGLES)])


def draw_lisa(ax, y_offset=0.0):
    """Draws the constellation on `ax`. Returns the list of
    (artist, local_transform) pairs that must follow the deformation
    (see set_strain)."""
    verts = spacecraft_positions()
    parts = []

    # undeformed reference triangle (static, not part of the deformed drawing)
    ref = Polygon(verts, closed=True, fill=False, ls="--", lw=0.7, ec="dimgray", zorder=1)
    ax.add_patch(ref)
    ref.set_transform(Affine2D().translate(0, y_offset) + ax.transData)

    # laser arms: wide faint glow + thin bright core
    for i, j in ((0, 1), (1, 2), (2, 0)):
        xs, ys = [verts[i, 0], verts[j, 0]], [verts[i, 1], verts[j, 1]]
        for lw, alpha in ((6.0, 0.18), (1.3, 0.95)):
            line = Line2D(xs, ys, color=LASER_COLOR, lw=lw, alpha=alpha,
                          solid_capstyle="round", zorder=2)
            ax.add_line(line)
            parts.append((line, Affine2D()))

    # spacecraft: two telescopes (one per arm) + cylindrical body + solar panel top
    for k, v in enumerate(verts):
        for j in range(3):
            if j == k:
                continue
            d = verts[j] - v
            ang = np.arctan2(d[1], d[0])
            tel = Rectangle((0.5 * SC_RADIUS, -0.03), 0.12, 0.06,
                            fc=TELESCOPE_COLOR, ec="white", lw=0.5, zorder=3)
            ax.add_patch(tel)
            parts.append((tel, Affine2D().rotate(ang).translate(*v)))

        body = Circle(tuple(v), SC_RADIUS, fc=BODY_COLOR, ec=TRIM_COLOR, lw=1.5, zorder=4)
        cells = RegularPolygon(tuple(v), 6, radius=0.62 * SC_RADIUS,
                               fc=PANEL_COLOR, ec="white", lw=0.5, zorder=5)
        for p in (body, cells):
            ax.add_patch(p)
            parts.append((p, Affine2D()))

    return parts


def set_strain(parts, ax, mat, y_offset=0.0):
    """Applies the 3x3 strain matrix `mat` (apple.strain_matrix) to every
    deformable artist of the drawing, then shifts it vertically."""
    deform = Affine2D(mat).translate(0, y_offset) + ax.transData
    for artist, local in parts:
        artist.set_transform(local + deform)
    return [artist for artist, _ in parts]
