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
		int w_idx = (int)(SQRT(ABS(corrected_vis_uvw_coords[i].z * config->w_scale)) + 0.5);
		int half_support = kernel_supports[w_idx].x;
		int w_offset = kernel_supports[w_idx].y;

		PRECISION2 grid_pos = {.x = corrected_vis_uvw_coords[i].x * config->uv_scale, .y = corrected_vis_uvw_coords[i].y * config->uv_scale};

		PRECISION conjugate = (corrected_vis_uvw_coords[i].z < 0.0) ? -1.0 : 1.0;

		float comm_norm = 0.f;

		for(int v = CEIL(grid_pos.y - half_support); v < CEIL(grid_pos.y + half_support); ++v)
		{
			int corrected_v = v;
			if(v < -grid_center){
				corrected_v = (-grid_center - v) - grid_center;
			}
			else if(v >= grid_center){
				corrected_v = (grid_center - v) + grid_center;
			}

			int2 kernel_idx;
			kernel_idx.y = abs((int)ROUND((corrected_v - grid_pos.y) * OVERSAMPLING_FACTOR));

			for(int u = CEIL(grid_pos.x - half_support); u < CEIL(grid_pos.x + half_support); ++u)
			{
				int corrected_u = u;

				if(u < -grid_center){
					corrected_u = (-grid_center - u) - grid_center;
				}
				else if(u >= grid_center){
					corrected_u = (grid_center - u) + grid_center;
				}

				kernel_idx.x = abs((int)ROUND((corrected_u - grid_pos.x) * OVERSAMPLING_FACTOR));
				int k_idx = w_offset + kernel_idx.y * (half_support + 1) * OVERSAMPLING_FACTOR + kernel_idx.x;

				PRECISION2 kernel_sample = kernels[k_idx];
				kernel_sample.y *= conjugate;

				int grid_idx = (corrected_v + grid_center) * GRID_SIZE + (corrected_u + grid_center);

				PRECISION2 prod = complex_mult(input_grid[grid_idx], kernel_sample);

				comm_norm += kernel_sample.x * kernel_sample.x + kernel_sample.y * kernel_sample.y;

				output_visibilities[i].x += prod.x;
				output_visibilities[i].y += prod.y;
			}
		}

		comm_norm = sqrt(comm_norm);

		output_visibilities[i].x = comm_norm < 1e-5f ? output_visibilities[i].x / 1e-5f : output_visibilities[i].x / comm_norm;
		output_visibilities[i].y = comm_norm < 1e-5f ? output_visibilities[i].x / 1e-5f : output_visibilities[i].y / comm_norm;

		//comment this portion out once i find some proper degridding kernels
		int x = (int)(grid_pos.x + (PRECISION)grid_center + 0.5);
		int y = (int)(grid_pos.y + (PRECISION)grid_center + 0.5);

		int idx = x + y * GRID_SIZE;

		output_visibilities[i].x = input_grid[idx].x;
		output_visibilities[i].y = input_grid[idx].y;
	}
}


void load_image_from_file(PRECISION *image, const char *file_name, int start_x, int start_y, int range_x, int range_y) {
	FILE *f = fopen(file_name, "r");
	if (f == NULL) {
		printf(">>> ERROR: Unable to open file %s...\n\n", file_name);
		return;
	}

	char line[1024];
	int row = 0;

	while (fgets(line, sizeof(line), f) && row < range_y) {
		char *token = strtok(line, ", ");
		int col = 0;
		while (token != NULL && col < range_x) {
#if SINGLE_PRECISION
			image[(start_y + row) * range_x + (start_x + col)] = (PRECISION)strtof(token, NULL);
#else
			image[(start_y + row) * range_x + (start_x + col)] = (PRECISION)strtod(token, NULL);
#endif
			token = strtok(NULL, ", ");
			col++;
		}
		row++;
	}

	fclose(f);
	printf("UPDATE >>> Image loaded from %s\n", file_name);
}

void convert_image_to_grid(PRECISION* final_image, PRECISION2* input_grid, int image_size, int grid_size) {
	float scale_factor = (float)grid_size / (float)image_size;  // Facteur de mise à l'échelle

	for (int y = 0; y < grid_size; ++y) {
		for (int x = 0; x < grid_size; ++x) {
			// Calculer les indices correspondants dans l'image finale
			int img_x = (int)(x / scale_factor);
			int img_y = (int)(y / scale_factor);

			// Initialisation d'un 'float2' à partir de la valeur de final_image
			input_grid[y * grid_size + x].x = final_image[img_y * image_size + img_x];  // La première composante
			input_grid[y * grid_size + x].y = 0;  // La deuxième composante
		}
	}

	printf("UPDATE >>> Image converted to grid (input_grid) successfully\n");
}
void save_image_to_file(Config *config, PRECISION *image, const char *file_name, int start_x, int start_y, int range_x, int range_y, int cycle)
{
	char buffer[256];
	snprintf(buffer,255,"%scycle_%d_%s",config->output_path,cycle,file_name);
	printf("UPDATE >>> Attempting to save image to %s... \n\n", buffer);

	//MD5_Update(sizeof(PRECISION) * range_x * range_y, image);

	FILE *f = fopen(buffer, "w");

	if(f == NULL)
	{
		printf(">>> ERROR: Unable to save image to file %s, check file/folder structure exists...\n\n", buffer);
		return;
	}

	for(int row = start_y; row < start_y + range_y; ++row)
	{
		for(int col = start_x; col < start_x + range_x; ++col)
		{
			PRECISION pixel = image[row * range_x + col];

#if SINGLE_PRECISION
			fprintf(f, "%f, ", pixel);
#else
			fprintf(f, "%lf, ", pixel);
#endif
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <csv_path>\n", argv[0]);
        return 1;
    }
	const char *image = argv[1];

	int GRID_SIZE = 2560;
	int NUM_VISIBILITIES = 3924480;
	int NUM_KERNEL = 17;
	int TOTAL_KERNEL_SAMPLES = 108800;
	int OVERSAMPLING_FACTOR = 16;

	PRECISION* final_image = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE*GRID_SIZE);
	PRECISION2* kernels = (PRECISION2*)malloc(sizeof(PRECISION2) * TOTAL_KERNEL_SAMPLES);
	int2* kernel_supports = (int2*)malloc(sizeof(int2) * NUM_KERNEL);
	PRECISION2* input_grid = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
	PRECISION3* vis_uvw_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	PRECISION3* corrected_vis_uvw_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	for (int i = 0; i < NUM_VISIBILITIES; ++i) {
		corrected_vis_uvw_coords[i].x = 0;
		corrected_vis_uvw_coords[i].y = 0;
		corrected_vis_uvw_coords[i].z = 0;
	}

	int* num_corrected_visibilities = (int*)malloc(sizeof(int) * 1);
	for(int i = 0; i < 1; ++i){
		num_corrected_visibilities[i] = 0;
	}
	PRECISION2* output_visibilities = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);

	Config config;
	config.max_w = 1895.410847844;
	config.w_scale = pow(NUM_KERNEL - 1, 2.0) / config.max_w;
	config.cell_size = 8.52211548825356E-06;
	config.uv_scale = config.cell_size * GRID_SIZE;
	config.weak_source_percent_gc = 0;
	config.weak_source_percent_img = 0;
	//config.psf_max_value = 1.f;
	config.visibility_source_file = "vis.csv";
	config.output_path = "/path";
	config.final_image_output = image;



	load_image_from_file( final_image, config.final_image_output, 0, 0, GRID_SIZE, GRID_SIZE);

	convert_image_to_grid(final_image, input_grid, GRID_SIZE, GRID_SIZE);

	degridding_kernel_host_set_up( NUM_KERNEL, TOTAL_KERNEL_SAMPLES, &config, kernel_supports,  kernels);

	//correction_set_up(GRID_SIZE, input_grid);

	// visibility_host_set_up(NUM_VISIBILITIES, &config, vis_uvw_coords, measured_vis);



	std_degridding(GRID_SIZE, NUM_VISIBILITIES, NUM_KERNEL, TOTAL_KERNEL_SAMPLES,OVERSAMPLING_FACTOR, kernels, kernel_supports, input_grid,vis_uvw_coords, num_corrected_visibilities, &config, output_visibilities);

	// faut mettre output dans config.visibility_source_file

	save_image_to_file(&config, final_image, config.visibility_source_file, 0, 0, GRID_SIZE, GRID_SIZE, 0);

	return 0;
}