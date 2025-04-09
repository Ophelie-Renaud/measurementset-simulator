import os
import sys
import glob
from astropy.io import fits
import matplotlib.pyplot as plt

NUM_MAJOR_CYCLE = 5  # Nombre de colonnes (cycles majeurs)

# Affiche une image FITS sur un subplot donné
def display_fits_image(ax, fits_file):
    with fits.open(fits_file) as hdulist:
        data = hdulist[0].data
    if data is not None:
        ax.imshow(data, cmap='viridis', origin='lower')
        ax.set_title(os.path.basename(fits_file), fontsize=8)
        ax.axis('off')

# Affiche une grille d'images FITS par type (ligne) et par cycle (colonne)
def display_images_by_type(base_dir, types):
    fig, axs = plt.subplots(len(types), NUM_MAJOR_CYCLE, figsize=(15, 3 * len(types)))
    axs = axs if len(types) > 1 else [axs]  # S'assurer que axs est bien une liste 2D

    for row, image_type in enumerate(types):
        files = sorted(glob.glob(os.path.join(base_dir, f"*_{image_type}.fits")))
        for col in range(NUM_MAJOR_CYCLE):
            if col < len(files):
                display_fits_image(axs[row][col], files[col])
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

