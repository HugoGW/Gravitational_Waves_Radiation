"""
Point d'entree : anime une binaire de trous noirs (masses proches de
GW150914) en inspiral, fusion et ringdown -- voir README.md.

Regle SAVE_MP4 ci-dessous pour choisir entre affichage en direct et
enregistrement video (necessite ffmpeg installe et accessible dans le PATH).
"""

from visualization import run

SAVE_MP4 = False           # True : enregistre l'animation en .mp4 au lieu de l'afficher
MP4_PATH = "GW_binary_merger.mp4"
FPS = 30

if __name__ == "__main__":
    run(save_mp4=SAVE_MP4, mp4_path=MP4_PATH, fps=FPS)
