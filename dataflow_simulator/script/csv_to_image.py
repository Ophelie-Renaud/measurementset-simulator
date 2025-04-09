import numpy as np
import matplotlib.pyplot as plt
import sys

def load_image_from_csv(csv_file, grid_size, space):
    # Charger les données à partir du fichier CSV
    data = np.loadtxt(csv_file, delimiter=',')
    
    if data.shape[1] != 2:
        raise ValueError("Le fichier CSV doit contenir deux colonnes : partie réelle et imaginaire.")
    
    # Reconstituer les données complexes
    complex_data = data[:, 0] + 1j * data[:, 1]
    total_points = complex_data.shape[0]
    expected_points = grid_size * grid_size

    # Ajustement automatique en fonction des dimensions de la grille
    if total_points < expected_points:
        print(f"Seulement {total_points} points, pas assez pour {grid_size}x{grid_size}")
        grid_size_final = int(np.floor(np.sqrt(total_points)))
        print(f"Ajustement automatique à {grid_size_final}x{grid_size_final}")
        complex_data = complex_data[:grid_size_final * grid_size_final]
    elif total_points > expected_points:
        print(f"{total_points} points : trop pour {grid_size}x{grid_size}, on tronque.")
        complex_data = complex_data[:expected_points]
        grid_size_final = grid_size
    else:
        grid_size_final = grid_size

    # Affichage des résultats
    display_image(complex_data.reshape((grid_size_final, grid_size_final)), csv_file, space)

def display_image(image_data, csv_file, space):
    abs_image = np.abs(image_data)
    phase_image = np.angle(image_data)
    
    # stat de base
    max_amplitude = np.max(abs_image)
    max_phase = np.max(np.abs(phase_image))
    mean_amp = np.mean(abs_image)
    std_amp = np.std(abs_image)
    mean_phase = np.mean(phase_image)
    std_phase = np.std(phase_image)

    # Création du subplot
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))  # 1 ligne, 2 colonnes
    
    # Si le plan est 'uv', appliquer l'échelle logarithmique
    if space == "uv":
# Affichage de l'amplitude
        im1 = axs[0].imshow(abs_image, cmap='viridis', origin='lower')
        axs[0].set_title(f"Amplitude : {csv_file}\n max: {max_amplitude:.2f},\n mean: {mean_amp:.2f}, std: {std_amp:.2f}")
        axs[0].set_xlabel("u")
        axs[0].set_ylabel("v")
        fig.colorbar(im1, ax=axs[0], label='Amplitude')

        # Affichage de la phase
        im2 = axs[1].imshow(phase_image, cmap='viridis', origin='lower')
        axs[1].set_title(f"Phase : {csv_file}\n max: {max_phase:.2f},\n mean: {mean_phase:.2f}, std: {std_phase:.2f}")
        axs[1].set_xlabel("u")
        axs[1].set_ylabel("v")
        fig.colorbar(im2, ax=axs[1], label='Phase')
    else:
# Affichage de l'amplitude
        im1 = axs[0].imshow(abs_image, cmap='viridis', origin='lower')
        axs[0].set_title(f"Amplitude : {csv_file}\n max: {max_amplitude:.2f},\n mean: {mean_amp:.2f}, std: {std_amp:.2f}")
        axs[0].set_xlabel("l")
        axs[0].set_ylabel("m")
        fig.colorbar(im1, ax=axs[0], label='Amplitude')

        # Affichage de la phase
        im2 = axs[1].imshow(phase_image, cmap='viridis', origin='lower')
        axs[1].set_title(f"Phase : {csv_file}\n max: {max_phase:.2f},\n mean: {mean_phase:.2f}, std: {std_phase:.2f}")
        axs[1].set_xlabel("l")
        axs[1].set_ylabel("m")
        fig.colorbar(im2, ax=axs[1], label='Phase')

    plt.xlabel("u" if space == "uv" else "l")
    plt.ylabel("v" if space == "uv" else "m")
    
    # Histogrammes
    axs[2].hist(abs_image.ravel(), bins=50, color='blue', alpha=0.7)
    axs[2].set_xlabel('Amplitude')
    axs[2].set_ylabel('Nombre de pixels')

    axs[3].hist(phase_image.ravel(), bins=50, color='orange', alpha=0.7)
    axs[3].set_xlabel('Phase (radians)')
    axs[3].set_ylabel('Nombre de pixels')
    

    plt.tight_layout()
    plt.show()
    
    center_y = image_data.shape[0] // 2
    amp_profile = abs_image[center_y, :]
    phase_profile = phase_image[center_y, :]

    fig_profile, axs_profile = plt.subplots(1, 2, figsize=(10, 1))
    axs_profile[0].plot(amp_profile, label='Amplitude')
    axs_profile[0].set_title('Profil amplitude (ligne centrale)')
    axs_profile[0].set_ylabel('Amplitude')

    axs_profile[1].plot(phase_profile, label='Phase', color='green')
    axs_profile[1].set_title('Profil phase (ligne centrale)')
    axs_profile[1].set_ylabel('Phase (rad)')
    axs_profile[1].set_xlabel('Position horizontale (pixels)')



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python csv_to_image.py <grid_size> <chemin_du_fichier_csv> <lm|uv>")
        sys.exit(1)

    try:
        grid_size = int(sys.argv[1])
    except ValueError:
        print("Erreur : grid_size doit être un entier.")
        sys.exit(1)

    csv_file = sys.argv[2]
    space = sys.argv[3].lower()

    if space not in ["lm", "uv"]:
        print("Erreur : Le troisième paramètre doit être 'lm' ou 'uv'.")
        sys.exit(1)

    load_image_from_csv(csv_file, grid_size, space)

