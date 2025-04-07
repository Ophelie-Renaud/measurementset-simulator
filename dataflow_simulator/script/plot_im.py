import numpy
import sys
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt

def write_nparr_to_fits(data, filename):
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

def display_fits_image(image_csv,image_fits):
    result = numpy.genfromtxt(image_csv, delimiter=",")[:,:-1]  # Ignorer la dernière colonne si nécessaire result = numpy.flip(result)
    write_nparr_to_fits(result, image_fits)

    # Charger les données du fichier FITS
    hdulist = fits.open(image_fits)
    data = hdulist[0].data  # Les données de l'image sont dans la première extension
    print(f"nombre de dimension du *.fits: {data.ndim}")
    hdulist.close()

    # Affichage de l'image avec matplotlib
    plt.imshow(data, cmap='viridis', origin='lower')
    plt.colorbar()  # Ajoute une barre de couleur pour l'échelle des intensités
    plt.title(f"Image {image_fits}")
    plt.show()
