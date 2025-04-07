import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 3e8  # Vitesse de la lumière en m/s


def kaiser_bessel_kernel(u, beta, width):
    """ Génère un noyau Kaiser-Bessel pour le gridding et degridding."""
    arg = np.maximum(0, 1 - (2 * u / width) ** 2)  # Empêche les valeurs négatives
    kernel = np.where(arg > 0, sp.i0(beta * np.sqrt(arg)) / sp.i0(beta), 0)
    return kernel

def generate_gridding_kernel(grid_size, num_kernels, oversampling, kernel_support, beta, width, baseline_max, frequency_hz):
    """ Génère une table de noyaux de gridding en 2D avec ajustement UV selon la fréquence et la baseline."""
    # Calcul de uvw_max
    uvw_max = baseline_max * frequency_hz / SPEED_OF_LIGHT
    
    # Définir les coordonnées U et V en fonction du support et de l'échantillonnage
    u = np.linspace(-kernel_support / 2, kernel_support / 2, grid_size // oversampling)
    v = np.linspace(-kernel_support / 2, kernel_support / 2, grid_size // oversampling)
    uu, vv = np.meshgrid(u, v)
    
    # Conversion en coordonnées UV en fonction de uvw_max
    uu = uu / uvw_max
    vv = vv / uvw_max
    
    r = np.sqrt(uu**2 + vv**2)
    kernel = kaiser_bessel_kernel(r, beta, width)
    kernel_complex = kernel + 1j * np.zeros_like(kernel)  # Partie réelle et imaginaire
    return kernel_complex / np.sum(kernel)  # Normalisation


def save_kernel_to_csv(kernel, filename_real="kernel_real.csv", filename_imag="kernel_imag.csv", filename_support="kernel_support.csv"):
    """ Sauvegarde la partie réelle, imaginaire et le support du noyau dans des fichiers CSV """
    real_part = np.real(kernel)
    imag_part = np.imag(kernel)
    support = np.indices(kernel.shape)[0]
    
    with open(filename_real, "w") as f_real, open(filename_imag, "w") as f_imag, open(filename_support, "w") as f_support:
        for i in range(kernel.shape[0]):
            for j in range(kernel.shape[1]):
                f_real.write(f"{real_part[i, j]},")
                f_imag.write(f"{imag_part[i, j]},")
                f_support.write(f"{support[i, j]},")
            f_real.write("\n")
            f_imag.write("\n")
            f_support.write("\n")

def plot_real_and_phase(kernel,kernel_support, title_real="Real Part", title_phase="Imag Part"):
    """ Affiche la partie réelle et la phase du noyau côte à côte """
    real_part = np.real(kernel)
    phase = np.angle(kernel)  # Phase en radians

    # Créer la figure avec deux sous-graphes côte à côte
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    # Afficher la partie réelle
    im_real = ax[0].imshow(real_part, extent=[-kernel_support / 2, kernel_support / 2, -kernel_support / 2, kernel_support / 2], cmap='viridis')
    ax[0].set_title(title_real)
    ax[0].set_xlabel("U")
    ax[0].set_ylabel("V")
    fig.colorbar(im_real, ax=ax[0], label="Amplitude")

    # Afficher la phase
    im_phase = ax[1].imshow(phase, extent=[-kernel_support / 2, kernel_support / 2, -kernel_support / 2, kernel_support / 2], cmap='viridis')
    ax[1].set_title(title_phase)
    ax[1].set_xlabel("U")
    ax[1].set_ylabel("V")
    fig.colorbar(im_phase, ax=ax[1], label="Phase (radians)")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python fits_to_csv.py <grid_size> <num_kernels> <oversampling_factor> <kernel_support> <baseline_max> <chemin_du_fichier_csv>")
        sys.exit(1)
    
    grid_size = int(sys.argv[1])
    num_kernels = int(sys.argv[2])
    oversampling_factor = int(sys.argv[3])
    kernel_support = int(sys.argv[4])
    baseline_max = int(sys.argv[5])

    csv_file = sys.argv[6]
    BETA = 2.0  # Paramètre de la fenêtre Kaiser-Bessel
    WIDTH = 4.0  # Largeur du noyau
    
    FREQUENCY_HZ = SPEED_OF_LIGHT/0.21
    
    kernel = generate_gridding_kernel(grid_size, num_kernels, oversampling_factor, kernel_support, BETA, WIDTH, baseline_max, FREQUENCY_HZ)
    plot_real_and_phase(kernel, kernel_support)
    save_kernel_to_csv(kernel)



