import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse

def generate_psf(grid_size, sigma=20, output_csv="psf.csv"):
    # Générer les coordonnées allant de 0 à grid_size-1
    x = np.linspace(0, grid_size-1, grid_size)
    y = np.linspace(0, grid_size-1, grid_size)
    X, Y = np.meshgrid(x, y)
    
    # PSF sous forme de gaussienne
    psf = np.exp(-((X - (grid_size-1)/2)**2 + (Y - (grid_size-1)/2)**2) / (2 * sigma**2))
    psf /= np.max(psf)  # Normalisation
    
    # Affichage
    plt.figure(figsize=(6, 6))
    plt.imshow(psf, cmap='viridis', extent=[0, grid_size-1, 0, grid_size-1])
    plt.colorbar(label='Amplitude')
    plt.title('Generated PSF')
    plt.xlabel('l')
    plt.ylabel('m')
    plt.show()
    
    # Export en CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f, delimiter=' ')
        for row in psf:
            writer.writerow(row)
    print(f"PSF exported to {output_csv}")

# Exécution
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python psf.py <grid_size> <chemin_du_fichier_csv>")
        sys.exit(1)
    
    grid_size = sys.argv[1]
    csv_file = sys.argv[2]
    plot_visibilities(grid_size,csv_file)

