"""
Entry point: animates a black-hole binary (masses close to GW150914)
through inspiral, merger and ringdown -- see README.md.

Set SAVE_MP4 below to choose between live display and video export
(requires ffmpeg installed and available on the PATH).
"""

from visualization import run

SAVE_MP4 = False           # True: saves the animation to .mp4 instead of displaying it
MP4_PATH = "GW_binary_merger.mp4"
FPS = 30

if __name__ == "__main__":
    run(save_mp4=SAVE_MP4, mp4_path=MP4_PATH, fps=FPS)
