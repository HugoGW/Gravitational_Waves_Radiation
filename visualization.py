"""
Figure, artists and animation loop. Can either display the animation live
(matplotlib window, infinite loop) or save it as .mp4 (one full cycle:
inspiral + merger + ringdown + pause).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.transforms import Affine2D

import physics as ph
import waveform as wf
import apple
import lisa

FRAME_DT_SIM = ph.T_CYCLE / 900.0  # M units advanced per frame (pacing)
H_YMAX = 1.3 * np.max(np.abs(wf.H_FULL_HPLUS))

# exaggeration factor for the apple/LISA deformation: the real deformation
# (~1e-21) is totally invisible; h+ and h_x are simply multiplied by this
# constant factor before building the deformation matrix (no nonlinear
# function, just a scale factor). It is normalized by the peak strain, so
# the maximum exaggerated strain is always APPLE_MAX_STRAIN, whatever the
# masses or the distance.
APPLE_MAX_STRAIN = 0.5
H_MAX = max(np.max(np.abs(wf.H_FULL_HPLUS)), np.max(np.abs(wf.H_FULL_HCROSS)))
APPLE_EXAGGERATION = APPLE_MAX_STRAIN / H_MAX

# size of the apple image (half-extent, before deformation) and vertical
# offset (downward) within its panel
APPLE_HALF_EXTENT = 1.3
APPLE_Y_OFFSET = -0.35
LISA_Y_OFFSET = -0.25

# purely graphical enlargement of the black-hole disks (the physical masses
# and Schwarzschild radii are unchanged)
BH_DISPLAY_SCALE = 4.0


def _time_unit(t_max_s):
    """Picks the most readable time unit (s, min, h or days) for a
    duration t_max_s given in seconds. Returns (label, seconds per unit)."""
    for label, scale in (("days", 86400.0), ("h", 3600.0), ("min", 60.0)):
        if t_max_s >= 2 * scale:
            return label, scale
    return "s", 1.0


# h+ time axis converted from GM/c^3 to a physical unit
T_UNIT_LABEL, T_UNIT_S = _time_unit(ph.T_RING_END * ph.TIME_UNIT_S)
T_TO_DISPLAY = ph.TIME_UNIT_S / T_UNIT_S
H_T_DISPLAY = wf.H_FULL_T * T_TO_DISPLAY


def build_figure():
    fig = plt.figure(figsize=(10, 18), facecolor="black")
    gs = fig.add_gridspec(3, 2, height_ratios=[4.6, 1.4, 1.7], hspace=0.4, wspace=0.1,
                          top=0.97, bottom=0.035)
    ax = fig.add_subplot(gs[0, :])
    ax_h = fig.add_subplot(gs[1, :])
    ax_apple = fig.add_subplot(gs[2, 0])
    ax_lisa = fig.add_subplot(gs[2, 1])
    ax.set_facecolor("black")
    ax.set_xlim(-wf.LIM, wf.LIM)
    ax.set_ylim(-wf.LIM, wf.LIM)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x\ [GM/c^2]$", color="white")
    ax.set_ylabel(r"$y\ [GM/c^2]$", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("dimgray")

    im = ax.imshow(np.zeros((wf.N_GRID, wf.N_GRID)), extent=[-wf.LIM, wf.LIM, -wf.LIM, wf.LIM],
                    origin="lower", cmap="turbo", vmin=0, vmax=0.5, zorder=1)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("radiated field intensity (a.u.)", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    cbar.outline.set_edgecolor("dimgray")
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

    # nudge the 2D panel (image + colorbar) slightly to the left of the figure
    shift_left = 0.045
    for a in (ax, cbar.ax):
        p = a.get_position()
        a.set_position([p.x0 - shift_left, p.y0, p.width, p.height])

    link_line, = ax.plot([], [], "-", color="white", lw=0.8, alpha=0.6, zorder=2)
    bh1_patch = Circle((0, 0), BH_DISPLAY_SCALE * ph.RS1, facecolor="black", edgecolor="none", zorder=3)
    bh2_patch = Circle((0, 0), BH_DISPLAY_SCALE * ph.RS2, facecolor="black", edgecolor="none", zorder=3)
    bhf_patch = Circle((0, 0), BH_DISPLAY_SCALE * ph.RS_F, facecolor="black", edgecolor="none", zorder=3)
    for patch in (bh1_patch, bh2_patch, bhf_patch):
        ax.add_patch(patch)

    ax_h.set_facecolor("black")
    ax_h.set_xlim(0, ph.T_RING_END * T_TO_DISPLAY)
    ax_h.set_ylim(-H_YMAX, H_YMAX)
    ax_h.axhline(0, color="dimgray", lw=0.6, zorder=1)
    ax_h.set_xlabel(f"$t$ [{T_UNIT_LABEL}]", color="white")
    ax_h.set_ylabel(r"$h_+$ (strain, dimensionless)", color="white")
    ax_h.tick_params(colors="white")
    ax_h.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax_h.yaxis.get_offset_text().set_color("white")
    for spine in ax_h.spines.values():
        spine.set_color("dimgray")
    ax_h.plot(H_T_DISPLAY, wf.H_FULL_HPLUS, "-", color="gray", lw=0.8, alpha=0.4, zorder=1)

    h_line, = ax_h.plot([], [], "-", color="white", lw=1.2, zorder=2)
    ax_h.text(0.01, 0.85, r"$h_+$", color="white", fontsize=13, transform=ax_h.transAxes)

    # --- bottom panels: (exaggerated) deformation of the apple and of the
    # LISA constellation under the h+/h_x wave ---
    for a in (ax_apple, ax_lisa):
        a.set_facecolor("black")
        a.set_xlim(-2.0, 2.0)
        a.set_ylim(-2.0, 2.0)
        a.set_aspect("equal")
        a.axis("off")
    fig.text(
        0.5, ax_apple.get_position().y1 + 0.012,
        fr"$h_+$ and $h_\times$ polarizations, exaggerated $\times {APPLE_EXAGGERATION:.0e}$",
        color="white", fontsize=11, ha="center",
    )
    ax_apple.set_title("test apple", color="white", fontsize=11, pad=2)
    ax_lisa.set_title("LISA constellation", color="white", fontsize=11, pad=2)
    ax_lisa.text(0.5, 0.0, f"arm length {lisa.ARM_LENGTH_KM / 1e6:g} million km (not to scale)",
                 color="gray", fontsize=8, ha="center", transform=ax_lisa.transAxes)

    apple_rgba = apple.load_apple_rgba()
    e = APPLE_HALF_EXTENT
    im_apple = ax_apple.imshow(apple_rgba, extent=[-e, e, -e, e], zorder=2,
                                interpolation="bilinear")
    im_apple.set_transform(Affine2D().translate(0, APPLE_Y_OFFSET) + ax_apple.transData)

    lisa_parts = lisa.draw_lisa(ax_lisa, y_offset=LISA_Y_OFFSET)

    artists = dict(im=im, link_line=link_line, bh1_patch=bh1_patch, bh2_patch=bh2_patch,
                    bhf_patch=bhf_patch, h_line=h_line, im_apple=im_apple, ax_apple=ax_apple,
                    lisa_parts=lisa_parts, ax_lisa=ax_lisa)
    return fig, artists


def make_update_functions(artists):
    im = artists["im"]
    link_line = artists["link_line"]
    bh1_patch = artists["bh1_patch"]
    bh2_patch = artists["bh2_patch"]
    bhf_patch = artists["bhf_patch"]
    h_line = artists["h_line"]
    im_apple = artists["im_apple"]
    ax_apple = artists["ax_apple"]
    lisa_parts = artists["lisa_parts"]
    ax_lisa = artists["ax_lisa"]

    def init():
        im.set_data(np.zeros((wf.N_GRID, wf.N_GRID)))
        link_line.set_data([], [])
        for patch in (bh1_patch, bh2_patch, bhf_patch):
            patch.set_radius(1e-6)
        h_line.set_data([], [])
        im_apple.set_transform(Affine2D().translate(0, APPLE_Y_OFFSET) + ax_apple.transData)
        # matplotlib always includes images (even animated ones) in full
        # redraws, i.e. in the background cached for blitting: hide the apple
        # here (init runs before every full redraw) so its undeformed copy
        # does not stay visible behind the deformed one
        im_apple.set_visible(False)
        lisa_artists = lisa.set_strain(lisa_parts, ax_lisa, np.eye(3), LISA_Y_OFFSET)
        return (im, link_line, bh1_patch, bh2_patch, bhf_patch, h_line, im_apple, *lisa_artists)

    def update(frame):
        t = frame * FRAME_DT_SIM
        t_local = t % ph.T_CYCLE

        im.set_data(wf.field_intensity(t))

        b = ph.bodies_now(t_local)
        if b["binary"]:
            bh1_patch.center = tuple(b["pos1"])
            bh2_patch.center = tuple(b["pos2"])
            bh1_patch.set_radius(max(BH_DISPLAY_SCALE * b["rs1"], 1e-6))
            bh2_patch.set_radius(max(BH_DISPLAY_SCALE * b["rs2"], 1e-6))
            bh1_patch.set_alpha(1.0)
            bh2_patch.set_alpha(1.0)
            link_line.set_data([b["pos1"][0], b["pos2"][0]], [b["pos1"][1], b["pos2"][1]])
        else:
            bh1_patch.set_alpha(0.0)
            bh2_patch.set_alpha(0.0)
            link_line.set_data([], [])

        bhf_patch.center = tuple(b["final_pos"])
        bhf_patch.set_radius(max(BH_DISPLAY_SCALE * b["final_rs"], 1e-6))
        bhf_patch.set_alpha(b["final_alpha"])

        n_done = np.searchsorted(wf.H_FULL_T, t_local, side="right")
        h_line.set_data(H_T_DISPLAY[:n_done], wf.H_FULL_HPLUS[:n_done])

        h_plus_now = np.interp(t_local, wf.H_FULL_T, wf.H_FULL_HPLUS)
        h_cross_now = np.interp(t_local, wf.H_FULL_T, wf.H_FULL_HCROSS)
        mat = apple.strain_matrix(h_plus_now, h_cross_now, APPLE_EXAGGERATION)
        im_apple.set_transform(Affine2D(mat).translate(0, APPLE_Y_OFFSET) + ax_apple.transData)
        im_apple.set_visible(True)
        lisa_artists = lisa.set_strain(lisa_parts, ax_lisa, mat, LISA_Y_OFFSET)

        return (im, link_line, bh1_patch, bh2_patch, bhf_patch, h_line, im_apple, *lisa_artists)

    return init, update


def run(save_mp4=False, mp4_path="GW_binary_merger.mp4", fps=30):
    """Displays the animation live (infinite loop), or saves it as .mp4 if
    save_mp4=True (a single inspiral + merger + ringdown + pause cycle)."""
    fig, artists = build_figure()
    init, update = make_update_functions(artists)

    if save_mp4:
        n_frames = int(np.ceil(ph.T_CYCLE / FRAME_DT_SIM))
        ani = animation.FuncAnimation(
            fig, update, frames=n_frames, init_func=init,
            blit=True, cache_frame_data=False,
        )
        ani.save(mp4_path, writer="ffmpeg", fps=fps, dpi=150)
        print(f"Animation saved to {mp4_path}")
    else:
        # mark the moving artists as animated (and hide the apple, see init)
        # *before* the first draw, so the background cached for blitting does
        # not contain their initial state
        for a in init():
            a.set_animated(True)
        ani = animation.FuncAnimation(
            fig, update, frames=100000, init_func=init,
            interval=1, blit=True, cache_frame_data=False,
        )
        plt.show()
