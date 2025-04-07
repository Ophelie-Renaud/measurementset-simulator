import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
import os
import sys

def load_visibilities_from_csv(csv_file):
    """
    Charge les visibilités à partir d'un fichier CSV.
    Le CSV est supposé avoir la structure : u, v, w, real, im, pol.
    """
    visibilities = {'u': [], 'v': [], 'w': [], 'real': [], 'im': [], 'pol': []}

    try:
        with open(csv_file, mode='r') as file:
            reader = csv.reader(file, delimiter=' ')
            next(reader)  # Sauter la première ligne (header)

            for row in reader:
                if len(row) == 6:  # Vérifier que la ligne contient bien 6 éléments
                    visibilities['u'].append(float(row[0]))
                    visibilities['v'].append(float(row[1]))
                    visibilities['w'].append(float(row[2]))
                    visibilities['real'].append(float(row[3]))
                    visibilities['im'].append(float(row[4]))
                    visibilities['pol'].append(float(row[5]))

        # Convertir les listes en numpy arrays
        return {key: np.array(value) for key, value in visibilities.items()}

    except FileNotFoundError:
        print(f"Erreur : le fichier '{csv_file}' n'existe pas.")
        exit(1)
    except Exception as e:
        print(f"Erreur lors du chargement du fichier CSV : {e}")
        exit(1)

def dft2d(grid):
    """
    Effectue une Transformée de Fourier discrète 2D brute (lente, complexité O(N^4)).
    """
    N, M = grid.shape
    image = np.zeros((N, M), dtype=np.complex64)

    for x in range(N):
        for y in range(M):
            sum_val = 0
            for k in range(N):
                for l in range(M):
                    exponent = 2j * np.pi * (x * k / N + y * l / M)
                    sum_val += grid[k, l] * np.exp(exponent)
            image[x, y] = sum_val

    return np.abs(image)  # Retourne seulement le module de l'image reconstruite

def reconstruct_image_dft(csv_file, image_size):
    """
    Reconstruit une image à partir des visibilités en utilisant une DFT explicite.
    
    :param csv_file: Chemin du fichier CSV contenant les visibilités
    :param image_size: Taille de la grille carrée pour la reconstruction
    """
    # Charger les visibilités
    df_cleaned = load_visibilities_from_csv(csv_file)

    # Extraire les coordonnées et les visibilités
    u, v = df_cleaned["u"], df_cleaned["v"]
    visibilities = df_cleaned["real"] + 1j * df_cleaned["im"]

    # Définir la grille de reconstruction
    grid = np.zeros((image_size, image_size), dtype=np.complex64)

    # Normaliser les coordonnées pour les mapper sur la grille
    u_min, u_max = u.min(), u.max()
    v_min, v_max = v.min(), v.max()

    # Éviter la division par zéro
    u_norm = np.zeros_like(u, dtype=int) if u_max == u_min else ((u - u_min) / (u_max - u_min) * (image_size - 1)).astype(int)
    v_norm = np.zeros_like(v, dtype=int) if v_max == v_min else ((v - v_min) / (v_max - v_min) * (image_size - 1)).astype(int)

    # Remplir la grille avec les visibilités (approximation par discrétisation)
    for i in range(len(u)):
        grid[u_norm[i], v_norm[i]] += visibilities[i]

    # Appliquer la transformée de Fourier discrète (DFT brute)
    image = dft2d(grid)

    # Afficher l'image reconstruite
    plt.figure(figsize=(6,6))
    plt.imshow(image, cmap="viridis", origin="lower")
    plt.colorbar(label="Intensité")
    plt.title(f"Image reconstruite avec DFT ({image_size}x{image_size})")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dft_reconstruction.py <chemin_du_fichier_vis.csv> <grid_size>")
        sys.exit(1)
    
    vis_csv = sys.argv[1]
    grid_size = int(sys.argv[2])
    
    # Exécution de la reconstruction d'image avec DFT
    reconstruct_image_dft(vis_csv, grid_size)

