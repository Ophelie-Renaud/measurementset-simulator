import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse

def load_visibilities(csv_file):
    """Charge les visibilités depuis un fichier CSV"""
    with open(csv_file, "r") as f:
        reader = csv.reader(f, delimiter=" ")
        num_vis = int(next(reader)[0])  # Première ligne = nombre de visibilités
        
        uv_list = []
        for row in reader:
            u, v, w, real, imag, pol = map(float, row)
            uv_list.append((u, v))
        
    return np.array(uv_list), num_vis
    
def generate_psf(grid_size, csv_file, output_csv="psf.csv"):
    """Génère une PSF à partir des visibilités et applique la FFT"""
    uv_data, num_vis = load_visibilities(csv_file)

    # Création d'un plan UV vierge
    uv_grid = np.zeros((grid_size, grid_size), dtype=np.complex64)

    # Transformation des coordonnées u,v en indices de la grille
    half_grid = grid_size // 2
    for u, v in uv_data:
        u_idx = int(np.round(u)) % grid_size  # Modulo pour rester dans la grille
        v_idx = int(np.round(v)) % grid_size
        uv_grid[v_idx, u_idx] = 1 + 0j  # Met toutes les visibilités à 1 pour la PSF

    # Appliquer une FFT 2D pour obtenir la PSF
    psf = np.fft.ifft2(uv_grid)
    psf = np.fft.fftshift(np.abs(psf))  # Centrage pour visualisation

    # Normalisation
    psf /= np.max(psf)
    
    # Appliquer l'échelle logarithmique (évite les valeurs nulles pour le log)
    #psf_log = np.log10(abs(psf) + 1e-6)

    # Affichage de la PSF
    plt.figure(figsize=(6, 6))
    plt.imshow(psf, cmap="viridis", extent=[-half_grid, half_grid, -half_grid, half_grid])
    plt.colorbar(label="Amplitude")
    plt.title("Point Spread Function (PSF)")
    plt.xlabel("l")
    plt.ylabel("m")
    plt.show()

    # Export en CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f, delimiter=" ")
        for row in psf:
            writer.writerow(row)

    print(f"✅ PSF générée et exportée dans {output_csv}")

# Exécution
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python psf.py <grid_size> <chemin_du_fichier_vis_csv> <chemin_du_fichier_psf_csv>")
        sys.exit(1)
    
    grid_size = sys.argv[1]
    vis_file = sys.argv[2]
    psf_file = sys.argv[3]
    generate_psf(grid_size,vis_file,psf_file)

