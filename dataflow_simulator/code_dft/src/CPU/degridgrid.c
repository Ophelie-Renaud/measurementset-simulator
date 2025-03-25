#include "degridgrid.h"

#include "map.h"

PRECISION2 complex_mult(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x * z2.x - z1.y * z2.y, z1.y * z2.x + z1.x * z2.y);
}

//takes in original visibilities and their continuous positions and corrects them to the oversampled finegrid positions. This function also performs reduction
//on visibilities that fall within the same cell and w-layer by adding together their values
void correct_to_finegrid(int NUM_VISIBILITIES, int GRID_SIZE, int OVERSAMPLING_FACTOR, int PERFORM_SIMPLIFICATION, PRECISION3* vis_uvw_coords, PRECISION2* input_visibilities,
		Config* config, PRECISION2* output_visibilities, PRECISION3* output_finegrid_vis_coords, int* num_output_visibilities){
	printf("Simplifying and correcting measured visibilities to fine grid\n");

	memset(output_visibilities, 0, sizeof(PRECISION2) * NUM_VISIBILITIES);
	*num_output_visibilities = 0;
	hashmap* finegrid_hash = hashmap_create();

	//assumes constant weighting, will need to change to accomodate for other weighting schemes (can just multiply with weight when adding)
	for(int i = 0; i < NUM_VISIBILITIES; ++i){
		//need 64 bit precision due to how many possible cells in a fine grid (eg. an image of 5000x5000 with an oversampling of 16 is already too much for normal 32bit
		//integers to handle
		unsigned long finegrid_x = (vis_uvw_coords[i].x * config->uv_scale + GRID_SIZE / 2) * OVERSAMPLING_FACTOR + 0.5f;
		unsigned long finegrid_y = (vis_uvw_coords[i].y * config->uv_scale + GRID_SIZE / 2) * OVERSAMPLING_FACTOR + 0.5f;

		unsigned long w_idx = (unsigned long)(SQRT(ABS(vis_uvw_coords[i].z * config->w_scale)) + 0.5);

		unsigned long idx = (finegrid_x + finegrid_y * GRID_SIZE * OVERSAMPLING_FACTOR) + w_idx * (GRID_SIZE * GRID_SIZE * OVERSAMPLING_FACTOR * OVERSAMPLING_FACTOR);

		uintptr_t result;

		if(PERFORM_SIMPLIFICATION == 1){
			if (hashmap_get(finegrid_hash, &idx, sizeof(unsigned long), &result))
			{
				output_visibilities[result].x += input_visibilities[i].x;
				output_visibilities[result].y += input_visibilities[i].y;
			}
			else
			{
				output_visibilities[*num_output_visibilities] = input_visibilities[i];
				output_finegrid_vis_coords[*num_output_visibilities].x = ((double)finegrid_x / OVERSAMPLING_FACTOR - GRID_SIZE / 2) / config->uv_scale;
				output_finegrid_vis_coords[*num_output_visibilities].y = ((double)finegrid_y / OVERSAMPLING_FACTOR - GRID_SIZE / 2) / config->uv_scale;
				output_finegrid_vis_coords[*num_output_visibilities].z = vis_uvw_coords[i].z;

				hashmap_set(finegrid_hash, &idx, sizeof(unsigned long), *num_output_visibilities);
				*num_output_visibilities = *num_output_visibilities + 1;
			}
		}
		else{
			output_visibilities[*num_output_visibilities] = input_visibilities[i];
			output_finegrid_vis_coords[*num_output_visibilities].x = ((double)finegrid_x / OVERSAMPLING_FACTOR - GRID_SIZE / 2) / config->uv_scale;
			output_finegrid_vis_coords[*num_output_visibilities].y = ((double)finegrid_y / OVERSAMPLING_FACTOR - GRID_SIZE / 2) / config->uv_scale;
			output_finegrid_vis_coords[*num_output_visibilities].z = vis_uvw_coords[i].z;

			hashmap_set(finegrid_hash, &idx, sizeof(unsigned long), *num_output_visibilities);
			*num_output_visibilities = *num_output_visibilities + 1;
		}

	}

	printf("Corrected %d visibilities to %d\n", NUM_VISIBILITIES, *num_output_visibilities);

	hashmap_free(finegrid_hash);
}

//gridding function assumes that the uvw coords are already corrected to the fine grid (either through calling std_degrid, or through calling correct_to_finegrid in the
//case where we are dealing with the original visibilities)
void std_gridding(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR, int BYPASS, int* maj_iter, int* num_corrected_visibilities,
			PRECISION2* kernels, int2* kernel_supports, PRECISION3* corrected_vis_uvw_coords, PRECISION2* visibilities, Config* config, PRECISION2* prev_grid, PRECISION2* output_grid){
	printf("Gridding visibilities\n");

	memset(output_grid, 0, GRID_SIZE * GRID_SIZE * sizeof(PRECISION2));

	if(BYPASS && *maj_iter > 0){
		memcpy(output_grid, prev_grid, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

		return;
	}

	int grid_center = GRID_SIZE / 2;

	for(int i = 0; i < *num_corrected_visibilities; ++i){
		int w_idx = (int)(SQRT(ABS(corrected_vis_uvw_coords[i].z * config->w_scale)) + 0.5);
		PRECISION2 grid_pos = {.x = corrected_vis_uvw_coords[i].x * config->uv_scale, .y = corrected_vis_uvw_coords[i].y * config->uv_scale};

		int half_support = kernel_supports[w_idx].x;

		PRECISION conjugate = (corrected_vis_uvw_coords[i].z < 0.0) ? -1.0 : 1.0;

		int w_offset = kernel_supports[w_idx].y;

		for(int v = CEIL(grid_pos.y - half_support); v < CEIL(grid_pos.y + half_support); ++v)
		{
			//out of bounds, assume zero padding
			if(v < -grid_center || v >= grid_center){
				continue;
			}

			int2 kernel_idx;
			kernel_idx.y = abs((int)ROUND((v - grid_pos.y) * OVERSAMPLING_FACTOR));

			for(int u = CEIL(grid_pos.x - half_support); u < CEIL(grid_pos.x + half_support); ++u)
			{
				//out of bounds, assume zero padding
				if(u < -grid_center || u >= grid_center){
					continue;
				}

				kernel_idx.x = abs((int)ROUND((u - grid_pos.x) * OVERSAMPLING_FACTOR));

				int k_idx = w_offset + kernel_idx.y * (half_support + 1) * OVERSAMPLING_FACTOR + kernel_idx.x;
				PRECISION2 kernel_sample = kernels[k_idx];
				kernel_sample.y *= conjugate;

				int grid_idx = (v + grid_center) * GRID_SIZE + (u + grid_center);

				PRECISION2 prod = complex_mult(visibilities[i], kernel_sample);

				output_grid[grid_idx].x += prod.x;
				output_grid[grid_idx].y += prod.y;
			}
		}
	}
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

void g2g_degridgrid(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_GRIDDING_KERNELS, int NUM_DEGRIDDING_KERNELS,
			int TOTAL_GRIDDING_KERNEL_SAMPLES, int TOTAL_DEGRIDDING_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR,
    		PRECISION2* gridding_kernels, int2* gridding_kernel_supports, PRECISION2* degridding_kernels, int2* degridding_kernel_supports, PRECISION2* input_grid,
			PRECISION3* corrected_vis_uvw_coords, int* num_corrected_visibilities, Config* config, PRECISION2* output_grid){
	printf("Performing Grid 2 Grid degridding gridding\n");

	memset(output_grid, 0, GRID_SIZE * GRID_SIZE * sizeof(PRECISION2));

	int grid_center = GRID_SIZE / 2;

	for(int i = 0; i < *num_corrected_visibilities; ++i){
		int w_idx = 0;//(int)(SQRT(ABS(corrected_vis_uvw_coords[i].z * config->w_scale)) + 0.5);
		int half_support_degridding = degridding_kernel_supports[w_idx].x;
		int degridding_w_offset = degridding_kernel_supports[w_idx].y;

		int degrid_half_support_size = half_support_degridding;
		PRECISION2 grid_pos = {.x = corrected_vis_uvw_coords[i].x * config->uv_scale, .y = corrected_vis_uvw_coords[i].y * config->uv_scale};

		PRECISION conjugate = (corrected_vis_uvw_coords[i].z < 0.0) ? -1.0 : 1.0;

		//degridding
		PRECISION2 degridded_visibility;
		degridded_visibility.x = 0;
		degridded_visibility.y = 0;

		float comm_norm = 0.f;

		for(int v = CEIL(grid_pos.y - degrid_half_support_size); v < CEIL(grid_pos.y + degrid_half_support_size); ++v)
		{
			//out of bounds, assume zero padding
			if(v < -grid_center || v >= grid_center){
				continue;
			}

			int2 kernel_idx;
			kernel_idx.y = abs((int)ROUND((v - grid_pos.y) * OVERSAMPLING_FACTOR));

			for(int u = CEIL(grid_pos.x - degrid_half_support_size); u < CEIL(grid_pos.x + degrid_half_support_size); ++u)
			{
				//out of bounds, assume zero padding
				if(u < -grid_center || u >= grid_center){
					continue;
				}

				kernel_idx.x = abs((int)ROUND((u - grid_pos.x) * OVERSAMPLING_FACTOR));
				int k_idx = degridding_w_offset + kernel_idx.y * (degrid_half_support_size + 1) * OVERSAMPLING_FACTOR + kernel_idx.x;

				PRECISION2 kernel_sample = degridding_kernels[k_idx];
				kernel_sample.y *= conjugate;

				int grid_idx = (v + grid_center) * GRID_SIZE + (u + grid_center);

				PRECISION2 prod = complex_mult(input_grid[grid_idx], kernel_sample);

				degridded_visibility.x += prod.x;
				degridded_visibility.y += prod.y;

				comm_norm += kernel_sample.x * kernel_sample.x + kernel_sample.y * kernel_sample.y;
			}
		}

		comm_norm = sqrt(comm_norm);

		//degridded_visibility.x = comm_norm < 1e-5f ? 0.f : degridded_visibility.x / comm_norm;
		//degridded_visibility.y = comm_norm < 1e-5f ? 0.f : degridded_visibility.y / comm_norm;



		//comment out when i have proper degridding kernels
		int x = (int)(grid_pos.x + (PRECISION)grid_center + 0.5);
		int y = (int)(grid_pos.y + (PRECISION)grid_center + 0.5);

		int idx = x + y * GRID_SIZE;

		if (input_grid && ((void *)&input_grid[idx] >= (void *)input_grid)) {
			degridded_visibility.x = input_grid[idx].x;
			degridded_visibility.y = input_grid[idx].y;
		} else {
			printf("Warning: input_grid is NULL or idx is out of bounds, setting default values.\n");
			degridded_visibility.x = 0;
			degridded_visibility.y = 0;
		}






		int half_support_gridding = gridding_kernel_supports[w_idx].x;
		int gridding_w_offset = gridding_kernel_supports[w_idx].y;

		//gridding
		for(int v = CEIL(grid_pos.y - half_support_gridding); v < CEIL(grid_pos.y + half_support_gridding); ++v)
		{
			if(v < -grid_center || v >= grid_center){
				continue;
			}

			int2 kernel_idx;
			kernel_idx.y = abs((int)ROUND((v - grid_pos.y) * OVERSAMPLING_FACTOR));

			for(int u = CEIL(grid_pos.x - half_support_gridding); u < CEIL(grid_pos.x + half_support_gridding); ++u)
			{
				if(u < -grid_center || u >= grid_center){
					continue;
				}

				kernel_idx.x = abs((int)ROUND((u - grid_pos.x) * OVERSAMPLING_FACTOR));

				int k_idx = gridding_w_offset + kernel_idx.y * (half_support_gridding + 1) * OVERSAMPLING_FACTOR + kernel_idx.x;
				PRECISION2 kernel_sample = gridding_kernels[k_idx];
				kernel_sample.y *= conjugate;

				int grid_idx = (v + grid_center) * GRID_SIZE + (u + grid_center);

				PRECISION2 prod = complex_mult(degridded_visibility, kernel_sample);

				output_grid[grid_idx].x += prod.x;
				output_grid[grid_idx].y += prod.y;
			}
		}
	}
}

void subtract_image_space(int GRID_SIZE, PRECISION* measurements, PRECISION* estimate, PRECISION* result){
	for(int y = 0; y < GRID_SIZE; ++y){
		for(int x = 0; x < GRID_SIZE; ++x){
			int idx = y * GRID_SIZE + x;
			result[idx] = measurements[idx] - estimate[idx];
		}
	}
}

void degridding_kernel_sink(int NUM_DEGRIDDING_KERNELS, int TOTAL_DEGRIDDING_KERNEL_SAMPLES, int2* supports, PRECISION2* kernels){
	return;
}
