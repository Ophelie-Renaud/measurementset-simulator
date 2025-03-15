#include <stdio.h>
#include "common.h"
#include "string.h"
#include "top.h"

float2 make_float2(float x, float y) {
	float2 result;
	result.x = x;
	result.y = y;
	return result;
}

PRECISION2 complex_mult(const PRECISION2 z1, const PRECISION2 z2)
{
	return make_float2(z1.x * z2.x - z1.y * z2.y, z1.y * z2.x + z1.x * z2.y);
}

void std_degridding(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR, PRECISION2* kernels,
		int2* kernel_supports, PRECISION2* input_grid, PRECISION3* corrected_vis_uvw_coords, int* num_corrected_visibilities, Config* config,
		PRECISION2* output_visibilities){
	printf("Degridding visibilities using FFT degridder\n");

	memset(output_visibilities, 0, sizeof(PRECISION2) * NUM_VISIBILITIES);

	int grid_center = GRID_SIZE / 2;

	for(int i = 0; i < *num_corrected_visibilities; ++i){
		// Calcul de l'indice w basé sur la coordonnée z des visibilités corrigées
		int w_idx = (int)(SQRT(ABS(corrected_vis_uvw_coords[i].z * config->w_scale)) + 0.5);
		w_idx = 0;

		// Récupère l'épaisseur du noyau et le décalage pour le plan w
		int half_support = kernel_supports[w_idx].x;
		int w_offset = kernel_supports[w_idx].y;

		// Calcul de la position sur la grille à partir des coordonnées UV
		PRECISION2 grid_pos = {.x = corrected_vis_uvw_coords[i].x * config->uv_scale, .y = corrected_vis_uvw_coords[i].y * config->uv_scale};
		printf("i: %d, grid_pos: (%f, %f)\n", i, grid_pos.x, grid_pos.y);

		// Détermine si la visibilité doit être conjuguée selon la coordonnée z
		PRECISION conjugate = (corrected_vis_uvw_coords[i].z < 0.0) ? -1.0 : 1.0;

		// Initialisation de la norme des coefficients du noyau pour la normalisation
		float comm_norm = 0.f;

        // Parcours des positions sur l'axe y dans la fenêtre du noyau
        for(int v = CEIL(grid_pos.y - half_support); v < CEIL(grid_pos.y + half_support); ++v)
		{
        	// Correction de l'indice v si nécessaire pour respecter les bornes de la grille
			int corrected_v = v;
			if(v < -grid_center){
				corrected_v = (-grid_center - v) - grid_center;
			}
			else if(v >= grid_center){
				corrected_v = (grid_center - v) + grid_center;
			}

        	// Détermination de l'indice du noyau sur l'axe y
			int2 kernel_idx;
			kernel_idx.y = abs((int)ROUND((corrected_v - grid_pos.y) * OVERSAMPLING_FACTOR));

        	// Parcours des positions sur l'axe x dans la fenêtre du noyau
			for(int u = CEIL(grid_pos.x - half_support); u < CEIL(grid_pos.x + half_support); ++u)
			{
				int corrected_u = u;

				if(u < -grid_center){
					corrected_u = (-grid_center - u) - grid_center;
				}
				else if(u >= grid_center){
					corrected_u = (grid_center - u) + grid_center;
				}

				// Détermination de l'indice du noyau sur l'axe x
				kernel_idx.x = abs((int)ROUND((corrected_u - grid_pos.x) * OVERSAMPLING_FACTOR));

				// Calcul de l'indice du noyau dans le tableau des noyaux
				int k_idx = w_offset + kernel_idx.y * (half_support + 1) * OVERSAMPLING_FACTOR + kernel_idx.x;

				// Récupère l'échantillon du noyau et applique la conjugaison si nécessaire
				PRECISION2 kernel_sample = kernels[k_idx];



				kernel_sample.y *= conjugate;

				// Calcul de l'indice de la grille d'entrée pour récupérer la valeur
				int grid_idx = (corrected_v + grid_center) * GRID_SIZE + (corrected_u + grid_center);
				//printf("grid_idx: %d, input_grid[grid_idx]: (%f, %f)\n", grid_idx, input_grid[grid_idx].x, input_grid[grid_idx].y);

				// Effectue le produit complexe entre la valeur de la grille et l'échantillon du noyau
				PRECISION2 prod = complex_mult(input_grid[grid_idx], kernel_sample);

				// Accumule la norme des échantillons du noyau pour la normalisation
				comm_norm += kernel_sample.x * kernel_sample.x + kernel_sample.y * kernel_sample.y;

				// Ajoute le produit au résultat final de la visibilité
				output_visibilities[i].x += prod.x;
				output_visibilities[i].y += prod.y;
			}
		}

		// Calcul de la norme des noyaux pour la normalisation
		comm_norm = sqrt(comm_norm);

		// Normalisation des visibilités pour éviter les valeurs trop petites
		output_visibilities[i].x = comm_norm < 1e-5f ? output_visibilities[i].x / 1e-5f : output_visibilities[i].x / comm_norm;
		output_visibilities[i].y = comm_norm < 1e-5f ? output_visibilities[i].x / 1e-5f : output_visibilities[i].y / comm_norm;

        // Portion temporaire pour copier directement la valeur de la grille dans la sortie (à commenter une fois que des noyaux corrects sont trouvés)
        // int x = (int)(grid_pos.x + (PRECISION)grid_center + 0.5);
		//int y = (int)(grid_pos.y + (PRECISION)grid_center + 0.5);

		//int idx = x + y * GRID_SIZE;

		//output_visibilities[i].x = input_grid[idx].x;
		//output_visibilities[i].y = input_grid[idx].y;
	}
	printf("UPDATE >>> Image degridded successfully\n");
}


void load_image_from_file(PRECISION2 *input_grid, const char *file_name) {
	FILE *f = fopen(file_name, "r");
	if (f == NULL) {
		perror(">>> ERROR");
		printf(">>> ERROR: Unable to open file %s...\n\n", file_name);
		return;
	}

	// Calcul du nombre d'éléments dans le fichier
	fseek(f, 0, SEEK_END);
	long file_size = ftell(f);
	rewind(f);

	// On suppose qu'il y a un nombre pair de valeurs dans le fichier (x, y)
	long num_elements = file_size / sizeof(float);  // Nombre total d'éléments (x et y compris)

	if (num_elements % 2 != 0) {
		printf("Le fichier doit contenir un nombre pair de valeurs (x, y).\n");
		fclose(f);
		return;
	}

	// Lecture des valeurs dans le fichier
	int i = 0;
	while (fscanf(f, "%f %f", &input_grid[i].x, &input_grid[i].y) == 2) {
		i++;
	}
	printf("UPDATE >>> Image loaded from %s\n", file_name);
}




void convert_vis_to_csv(int NUM_VISIBILITIES,PRECISION2* output_visibilities, PRECISION3* corrected_vis_uvw_coords,Config *config) {
	FILE* file = fopen(config->visibility_source_file, "w");
	if (file == NULL) {
		perror("Erreur lors de l'ouverture du fichier");
		exit(EXIT_FAILURE);
	}

	// Écrire le nombre de visibilités dans la première ligne
	fprintf(file, "%d\n", NUM_VISIBILITIES);

	// Écrire les visibilités u, v, real, imag dans le fichier CSV
	for (int i = 0; i < NUM_VISIBILITIES; i++) {
		// Assurez-vous que les champs sont correctement initialisés
		if (isnan(corrected_vis_uvw_coords[i].x) || isnan(corrected_vis_uvw_coords[i].y) || isnan(corrected_vis_uvw_coords[i].z) ||
			isnan(output_visibilities[i].x) || isnan(output_visibilities[i].y)) {
			fprintf(stderr, "Erreur : des données invalides (NaN) rencontrées à l'index %d.\n", i);
			continue; // Ou gérer l'erreur selon le besoin
			}
		// Écrire les données dans le fichier
		fprintf(file, "%.6f %.6f %.6f %.6f %.6f 1\n",
				corrected_vis_uvw_coords[i].x,
				corrected_vis_uvw_coords[i].y,
				corrected_vis_uvw_coords[i].z,
				output_visibilities[i].x,
				output_visibilities[i].y);
	}

	// Fermer le fichier
	fclose(file);
}


void vis_coord_set_up(int NUM_VISIBILITIES, int GRID_SIZE, PRECISION3* vis_uvw_coords, Config *config) {
	// Paramètres pour définir la plage UV
	const float MAX_UVW =300.0;

	// Option 1: Distribution uniforme (aléatoire) dans la plage
	/*for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		// Générer des coordonnées aléatoires dans un carré de taille MAX_UVW
		vis_uvw_coords[i].x = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
		vis_uvw_coords[i].y = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
		vis_uvw_coords[i].z = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
	}*/

	// Option 2: Distribution gaussienne autour du centre (phase center) - Exemple simple

	for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		vis_uvw_coords[i].x = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;
		vis_uvw_coords[i].y = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;
		vis_uvw_coords[i].z = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;
		float norm = sqrt(vis_uvw_coords[i].x * vis_uvw_coords[i].x +
						  vis_uvw_coords[i].y * vis_uvw_coords[i].y +
						  vis_uvw_coords[i].z * vis_uvw_coords[i].z);
		vis_uvw_coords[i].x /= norm;  // Normaliser
		vis_uvw_coords[i].y /= norm;
		vis_uvw_coords[i].z /= norm;
	}

}

int main(void) {

	// dimension du csv (depuis le .fits)
	int RANG_X = 512;
	int RANG_Y = 512;

	int NUM_RECEIVERS = 5;
	int NUM_BASELINE = NUM_RECEIVERS*(NUM_RECEIVERS-1)/2;
	int OVERSAMPLING_FACTOR = 16;
	int GRID_SIZE = 512;

	// si y'a plus de vis que de taille d'image d'entré ça couvre tout le plan uv
	int NUM_VISIBILITIES = NUM_BASELINE*OVERSAMPLING_FACTOR;
	int NUM_KERNEL = 17;
	int TOTAL_KERNEL_SAMPLES = NUM_KERNEL*OVERSAMPLING_FACTOR;

	int UV_LEN = 1000;




	PRECISION* final_image = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE*GRID_SIZE);
	PRECISION2* kernels = (PRECISION2*)malloc(sizeof(PRECISION2) * TOTAL_KERNEL_SAMPLES);
	int2* kernel_supports = (int2*)malloc(sizeof(int2) * NUM_KERNEL);
	PRECISION2* input_grid = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

	PRECISION3* vis_uvw_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	/*for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		vis_uvw_coords[i].x = (float)(rand() % 601) - 300;  // Aléatoire entre -300 et 300
		vis_uvw_coords[i].y = (float)(rand() % 601) - 300;  // Aléatoire entre -300 et 300
		vis_uvw_coords[i].z = (float)(rand() % 601) - 300;  // Aléatoire entre -300 et 300
	}*/
	PRECISION3* corrected_vis_uvw_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		corrected_vis_uvw_coords[i].x = (float)(rand() % 601) - 300;  // Aléatoire entre -30000 et 30000
		corrected_vis_uvw_coords[i].y = (float)(rand() % 601) - 300;   // Aléatoire entre -30000 et 30000
		corrected_vis_uvw_coords[i].z = (float)(rand() % 601) - 300;   // Aléatoire entre -30000 et 30000
	}

	int* num_corrected_visibilities = (int*)malloc(sizeof(int) * 1);
	for(int i = 0; i < 1; ++i){
		num_corrected_visibilities[i] = TOTAL_KERNEL_SAMPLES;
	}
	PRECISION2* output_visibilities = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	if (output_visibilities == NULL) {
		fprintf(stderr, "Erreur d'allocation mémoire\n");
		exit(EXIT_FAILURE);
	}
	//PRECISION2* measured_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);


	Config config;
	config.max_w = 1895.410847844;//osef = baseline_max * freq obs /celerite
	config.w_scale = pow(NUM_KERNEL - 1, 2.0) / config.max_w;
	config.cell_size = 1.953125; // info depuis le .fits: norme du vecteur (CDELT1, CDELT2)
	config.uv_scale =  GRID_SIZE/(2*300);

	config.visibility_source_file = "vis.csv";
	config.output_path = "";
	config.final_image_output = "image.csv";
	config.degridding_kernel_support_file = "config/wproj_manualconj_degridding_kernel_supports_x16.csv";
	config.degridding_kernel_imag_file = "config/wproj_manualconj_degridding_kernels_imag_x16.csv";
	config.degridding_kernel_real_file= "config/wproj_manualconj_degridding_kernels_real_x16.csv";


	load_image_from_file( input_grid, config.final_image_output);

	// Affichage des résultats
	for (int i = 0; i < 10; i++) {
		printf("Grille d'entrée %d: %.6f + %.6fi\n", i, input_grid[i].x, input_grid[i].y);
	}

	vis_coord_set_up(NUM_VISIBILITIES, GRID_SIZE,vis_uvw_coords, &config);


	degridding_kernel_host_set_up( NUM_KERNEL, TOTAL_KERNEL_SAMPLES, &config, kernel_supports,  kernels);

	/*for (int i=0;i<NUM_VISIBILITIES;i++) {
		output_visibilities[i].x = input_grid[i].x;
		output_visibilities[i].y = input_grid[i].y;
	}*/



	std_degridding(GRID_SIZE, NUM_VISIBILITIES, NUM_KERNEL, TOTAL_KERNEL_SAMPLES,OVERSAMPLING_FACTOR, kernels, kernel_supports, input_grid,vis_uvw_coords, num_corrected_visibilities, &config, output_visibilities);

	// Affichage des résultats
	for (int i = 0; i < 10; i++) {
		printf("Visibilité %d: %.6f + %.6fi\n", i, output_visibilities[i].x, output_visibilities[i].y);
	}

	convert_vis_to_csv(NUM_VISIBILITIES,output_visibilities, vis_uvw_coords, &config);

	//correction_set_up(GRID_SIZE, input_grid);

	//visibility_host_set_up(NUM_VISIBILITIES, &config, vis_uvw_coords, measured_vis);

	free(kernel_supports);
	free(kernels);
	free(input_grid);
	free(corrected_vis_uvw_coords);
	free(output_visibilities);
	free(num_corrected_visibilities);
	//free(measured_vis);

	return 0;
}
