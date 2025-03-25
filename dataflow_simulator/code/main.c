#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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

// Génération d'un nombre aléatoire suivant une distribution gaussienne (normale) avec moyenne 0 et écart-type 1
double randn(double mean, double stddev) {
	double u1 = ((double) rand() / RAND_MAX);  // Génère un nombre entre 0 et 1
	double u2 = ((double) rand() / RAND_MAX);  // Génère un autre nombre entre 0 et 1
	double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);  // Box-Muller transform
	return z0 * stddev + mean;  // Transformation pour avoir une moyenne et un écart-type donnés
}
double rayleigh(double sigma) {
	double u = (double) rand() / RAND_MAX;
	return sigma * sqrt(-2.0 * log(1.0 - u));
}

void std_degridding(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR, PRECISION2* kernels,
		int2* kernel_supports, PRECISION2* input_grid, PRECISION3* vis_uvw_coords, int* num_corrected_visibilities, Config* config,
		PRECISION2* output_visibilities){
	printf("Degridding visibilities using FFT degridder\n");

	memset(output_visibilities, 0, sizeof(PRECISION2) * NUM_VISIBILITIES);

	int grid_center = GRID_SIZE / 2;

	for(int i = 0; i < *num_corrected_visibilities; ++i){
		// Calcul de l'indice w basé sur la coordonnée z des visibilités corrigées
		int w_idx = (int)(SQRT(ABS(vis_uvw_coords[i].z * config->w_scale)) + 0.5);
		// we don't care first on w
		w_idx = 0;

		// Récupère l'épaisseur du noyau et le décalage pour le plan w
		int half_support = kernel_supports[w_idx].x;
		int w_offset = kernel_supports[w_idx].y;

		// Calcul de la position sur la grille à partir des coordonnées UV
		PRECISION2 grid_pos = {.x = vis_uvw_coords[i].x * config->uv_scale, .y = vis_uvw_coords[i].y * config->uv_scale};
		//printf("i: %d, grid_pos: (%f, %f)\n", i, grid_pos.x, grid_pos.y);

		// Détermine si la visibilité doit être conjuguée selon la coordonnée z
		PRECISION conjugate = (vis_uvw_coords[i].z < 0.0) ? -1.0 : 1.0;

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


void load_image_from_file(PRECISION2 *input_grid, Config* config) {
	const char *file_name = config->final_image_output;
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
		printf("Nombre impair de valeurs détecté, suppression du dernier élément.\n");
		num_elements--;
	}

	// Lecture des valeurs dans le fichier
	int i = 0;
	while (fscanf(f, "%f, %f,", &input_grid[i].x, &input_grid[i].y) == 2) {
		i++;
	}
	printf("UPDATE >>> Image loaded from %s\n", file_name);
}




void handle_file_error(FILE *file, const char *message) {
	if (file) {
		fclose(file);  // Ferme le fichier si ouvert
	}
	perror(message);
	exit(EXIT_FAILURE);
}

void convert_vis_to_csv(int NUM_VISIBILITIES, PRECISION2* output_visibilities, PRECISION3* corrected_vis_uvw_coords, Config *config) {
	FILE* file = fopen(config->visibility_source_file, "w");
	if (file == NULL) {
		handle_file_error(file, "Erreur lors de l'ouverture du fichier");
	}

	// Écrire le nombre de visibilités dans la première ligne
	if (fprintf(file, "%d\n", NUM_VISIBILITIES) < 0) {
		handle_file_error(file, "Erreur lors de l'écriture du nombre de visibilités");
	}

	for (int i = 0; i < NUM_VISIBILITIES; i++) {
		// Vérification des valeurs NaN
		if (isnan(corrected_vis_uvw_coords[i].x) || isnan(corrected_vis_uvw_coords[i].y) || isnan(corrected_vis_uvw_coords[i].z) ||
			isnan(output_visibilities[i].x) || isnan(output_visibilities[i].y)) {
			fprintf(stderr, "Erreur : des données invalides (NaN) rencontrées à l'index %d.\n", i);
			handle_file_error(file, "Erreur lors de la rencontre de données invalides");
			}

		// Écriture des données dans le fichier CSV
		if (fprintf(file, "%.6f %.6f %.6f %.6f %.6f 1\n",
					corrected_vis_uvw_coords[i].x,
					corrected_vis_uvw_coords[i].y,
					corrected_vis_uvw_coords[i].z,
					output_visibilities[i].x,
					output_visibilities[i].y) < 0) {
			handle_file_error(file, "Erreur lors de l'écriture dans le fichier");
					}
	}

	// Fermer le fichier
	if (fclose(file) != 0) {
		handle_file_error(file, "Erreur lors de la fermeture du fichier");
	}

	printf("UPDATE >>> write csv to %s\n", config->visibility_source_file);
}




void vis_coord_set_up(int NUM_VISIBILITIES, int GRID_SIZE, PRECISION3* vis_uvw_coords, Config *config) {
	// Paramètres pour définir la plage UV
	const float MAX_UVW =config->max_w;
	double meters_to_wavelengths = SPEED_OF_LIGHT/config->frequency_hz ; // lambda = c/f

	const double sigma = MAX_UVW / 30.0;  // Ajuster selon vos besoins

	// Option 1: Distribution uniforme (aléatoire) dans la plage
	/*for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		// Générer des coordonnées aléatoires dans un carré de taille MAX_UVW
		vis_uvw_coords[i].x = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
		vis_uvw_coords[i].y = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
		vis_uvw_coords[i].z = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;  // Aléatoire entre -MAX_UVW et +MAX_UVW
	}*/

	// Option 2: Distribution gaussienne autour du centre (phase center) - Box-Muller simple

	for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		// Générer une coordonnée x, y, z suivant une distribution normale
		vis_uvw_coords[i].x = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;
		vis_uvw_coords[i].y = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;
		vis_uvw_coords[i].z = (float)(rand() % (int)(2 * MAX_UVW + 1)) - MAX_UVW;

		// Appliquer une distribution gaussienne avec une plus grande densité autour du centre
		vis_uvw_coords[i].x = (float)randn(0.0, MAX_UVW);  // Génération d'un nombre aléatoire suivant une loi normale
		vis_uvw_coords[i].y = (float)randn(0.0, MAX_UVW);
		vis_uvw_coords[i].z = (float)randn(0.0, MAX_UVW);


	}

	// Option 3: Distribution gaussienne autour du centre (phase center) - Rayleigh simple

	/*for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		// Générer une coordonnée x, y, z suivant une distribution normale
		vis_uvw_coords[i].x = (float)rayleigh(sigma);
		vis_uvw_coords[i].y = (float)rayleigh(sigma);
		vis_uvw_coords[i].z = (float)rayleigh(sigma);

	}*/
	for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		vis_uvw_coords[i].x *= meters_to_wavelengths;  // Génération d'un nombre aléatoire suivant une loi normale
		vis_uvw_coords[i].y *= meters_to_wavelengths;
		vis_uvw_coords[i].z *= meters_to_wavelengths;
	}

}

int main(void) {

	int FOV_DEGREES = 1; //champs de vue
	int NUM_CHANNEL = 1;//nombre de canaux de frequence
	int NUM_POL = 2;// nombre de polarisation
	int TIMING_SAMPLE = 10000;

	int NUM_RECEIVERS = 5;//nombre d'antennes
	int NUM_BASELINE = NUM_RECEIVERS*(NUM_RECEIVERS-1)/2; // nombre de paire d'antennes
	int OVERSAMPLING_FACTOR = 16; // nombre de fois que les données sont surechantillonée
	int GRID_SIZE = 512; // taille de l'image fits "true sky"
	int NUM_VISIBILITIES = NUM_BASELINE*TIMING_SAMPLE*NUM_CHANNEL*NUM_POL; // greater than (GRID_SIZE/cell_size)²
	int NUM_KERNEL = 17;//nombre de noyaux de convolution
	int k = 100;//GRID_SIZE doit etr au moins 2 fois superieur à la taille du noyau kernel_size = 2*half+1
	int half_support =(int)( (GRID_SIZE/2 -1)/k);
	int TOTAL_KERNEL_SAMPLES = (int)(pow(2*half_support*OVERSAMPLING_FACTOR,2)*NUM_KERNEL); //c'est recalculé dans le processus mais faut mettre une taille suffisante



	PRECISION2* kernels = (PRECISION2*)malloc(sizeof(PRECISION2) * TOTAL_KERNEL_SAMPLES);
	int2* kernel_supports = (int2*)malloc(sizeof(int2) * NUM_KERNEL);
	PRECISION2* input_grid = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
	PRECISION3* vis_uvw_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	int* num_corrected_visibilities = (int*)malloc(sizeof(int) * 1);
	PRECISION2* output_visibilities = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	Config config;


	//initialize values (pas necessaire mais au cas ou)

	memset(kernels, 0, sizeof(PRECISION2)*TOTAL_KERNEL_SAMPLES);
	memset(kernel_supports, 0, sizeof(int2) * NUM_KERNEL);
	memset(input_grid, 0, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
	memset(vis_uvw_coords, 0, sizeof(PRECISION3) * NUM_VISIBILITIES);
	memset(output_visibilities, 0, sizeof(PRECISION2) * NUM_VISIBILITIES);



	for(int i = 0; i < 1; ++i){
		num_corrected_visibilities[i] = NUM_VISIBILITIES;
	}

	if (kernels == NULL || kernel_supports == NULL || input_grid == NULL ||
	vis_uvw_coords == NULL  ||
	num_corrected_visibilities == NULL || output_visibilities == NULL) {
		fprintf(stderr, "Erreur d'allocation mémoire\n");
		exit(EXIT_FAILURE);
	}

	config.frequency_hz = 29979.0; // observed frequency --> osef
	int baseline_max = 100;
	config.max_w = baseline_max*config.frequency_hz/SPEED_OF_LIGHT;// baseline_max * freq obs /celerite mais osef en fait
	config.w_scale = pow(NUM_KERNEL - 1, 2.0) / config.max_w;
	config.cell_size = (FOV_DEGREES * M_PI) / (180.0 * GRID_SIZE);// lower than 1/2.f_max
	config.uv_scale =  config.cell_size*GRID_SIZE;

	config.num_baselines = NUM_BASELINE;
	config.num_kernels = NUM_KERNEL;
	config.total_kernel_samples = TOTAL_KERNEL_SAMPLES;

	config.visibility_source_file = "vis.csv";
	config.output_path = "";
	config.final_image_output = "image.csv";
	config.degridding_kernel_support_file = "config/wproj_manualconj_degridding_kernel_supports_x16.csv";
	config.degridding_kernel_imag_file = "config/wproj_manualconj_degridding_kernels_imag_x16.csv";
	config.degridding_kernel_real_file= "config/wproj_manualconj_degridding_kernels_real_x16.csv";
	config.oversampling = OVERSAMPLING_FACTOR;


	load_image_from_file( input_grid, &config);

	// Affichage des résultats
	for (int i = 0; i < 5; i++) {
		printf("Grille d'entrée %d: %.6f + %.6fi\n", i, input_grid[i].x, input_grid[i].y);
	}

	vis_coord_set_up(NUM_VISIBILITIES, GRID_SIZE,vis_uvw_coords, &config);


	degridding_kernel_host_set_up( NUM_KERNEL, TOTAL_KERNEL_SAMPLES, &config, kernel_supports,  kernels);
	// Affichage des résultats
	for (int i = 0; i < 5; i++) {
		printf("kernel_supports %d: %d + %di\n", i, kernel_supports[i].x, kernel_supports[i].y);
	}
	for (int i = 0; i < 5; i++) {
		printf("kernels %d: %.6f + %.6fi\n", i, kernels[i].x, kernels[i].y);
	}

	/*for (int i=0;i<NUM_VISIBILITIES;i++) {
		output_visibilities[i].x = input_grid[i].x;
		output_visibilities[i].y = input_grid[i].y;
	}*/



	std_degridding(GRID_SIZE, NUM_VISIBILITIES, NUM_KERNEL, TOTAL_KERNEL_SAMPLES,OVERSAMPLING_FACTOR, kernels, kernel_supports, input_grid,vis_uvw_coords, num_corrected_visibilities, &config, output_visibilities);

	// Affichage des résultats
	for (int i = 0; i < 5; i++) {
		printf("Visibilité %d: %.6f + %.6fi\n", i, output_visibilities[i].x, output_visibilities[i].y);
	}

	convert_vis_to_csv(NUM_VISIBILITIES,output_visibilities, vis_uvw_coords, &config);


	free(kernel_supports);
	kernel_supports = NULL;
	free(kernels);
	kernels = NULL;
	free(input_grid);
	input_grid = NULL;
	free(output_visibilities);
	output_visibilities = NULL;
	free(num_corrected_visibilities);
	num_corrected_visibilities = NULL;
	free(vis_uvw_coords);
	vis_uvw_coords = NULL;

	return 0;
}
