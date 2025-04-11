import os
import sys
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np



NUM_MAJOR_CYCLE = 5  # Nombre de colonnes (cycles majeurs)

# Affiche une image FITS sur un subplot donné avec stats et barre de couleur
def display_fits_image(ax, fits_file):
    with fits.open(fits_file) as hdulist:
        data = hdulist[0].data

    if data is not None and np.isfinite(data).any():
        # Nettoyage : on remplace les valeurs non finies par des zéros (ou on les masque)
        data_clean = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

        im = ax.imshow(data_clean, cmap='viridis', origin='lower')
        mean_val = np.mean(data_clean)
        std_val = np.std(data_clean)
        max_val = np.max(data_clean)
        max_pos = np.unravel_index(np.argmax(data_clean), data_clean.shape)
        stats_text = f"mean={mean_val:.2e}\nstd={std_val:.2e}\nmax={max_val:.2e}\nmax_pos={max_pos}"
        ax.set_title(os.path.basename(fits_file), fontsize=8)
        ax.text(0.01, 0.99, stats_text, fontsize=6, color='white',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(facecolor='black', alpha=0.5, boxstyle='round'))
        ax.axis('off')
        return im
    else:
        ax.text(0.5, 0.5, "Invalid or empty data", fontsize=6,
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return None

# Affiche une grille d'images FITS par type (ligne) et par cycle (colonne)
def display_images_by_type(base_dir, types):
    fig, axs = plt.subplots(len(types), NUM_MAJOR_CYCLE, figsize=(15, 3 * len(types)))
    axs = np.array(axs)  # S'assurer que axs est un tableau NumPy 2D
    if axs.ndim == 1:
        axs = axs.reshape(1, -1)

    for row, image_type in enumerate(types):
        files = sorted(glob.glob(os.path.join(base_dir, f"*_{image_type}.fits")))
        for col in range(NUM_MAJOR_CYCLE):
            if col < len(files):
                im = display_fits_image(axs[row][col], files[col])
                if im is not None:
                    cbar = fig.colorbar(im, ax=axs[row][col], fraction=0.046, pad=0.04)
                    cbar.ax.tick_params(labelsize=6)
            else:
                axs[row][col].axis('off')

    plt.tight_layout()
    plt.show()

# Point d'entrée principal
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python display_fits_grid.py <repertoire_fits> <type1> <type2> ...")
        print("Exemple : python display_fits_grid.py code_dft/data/fits model dirty_psf deconvolved")
        sys.exit(1)

    base_dir = sys.argv[1]
    types = sys.argv[2:]

    display_images_by_type(base_dir, types)

