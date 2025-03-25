#include "splitmerge.h"

void split_visibilities(IN Config* config, IN int NUM_VISIBILITIES, IN int SLICES, IN PRECISION2* measured_visibilities, IN PRECISION3* vis_coords,
		OUT PRECISION2* visibilities, OUT PRECISION3* out_vis_coords, OUT Config* output_configs, OUT Config* output_config){
	memcpy(visibilities, measured_visibilities, NUM_VISIBILITIES * sizeof(PRECISION2));
	memcpy(out_vis_coords, vis_coords, NUM_VISIBILITIES * sizeof(PRECISION3));

	for(int i = 0; i < SLICES; ++i){
		output_configs[i] = *config;
	}
	*output_config = *config;
}

void merge_gridded_visibilities(IN int GRID_SIZE, IN int SLICES, IN PRECISION2* grids, OUT PRECISION2* output_grid){
	for(unsigned int i = 0; i < GRID_SIZE * GRID_SIZE; ++i){
		PRECISION2 curr_vis = {.x = 0, .y = 0};
		for(unsigned int j = 0; j < SLICES; ++j){
			unsigned int idx = j * GRID_SIZE * GRID_SIZE + i;
			curr_vis.x += grids[idx].x;
			curr_vis.y += grids[idx].y;
		}
		output_grid[i] = curr_vis;
	}
}

void merge_configs(IN Config* in_configs, OUT Config* out_config){
	*out_config = in_configs[0];
}

void split_fftdegrid_visibilities(int NUM_VISIBILITIES, int SLICES, int GRID_SIZE, IN Config* config, IN PRECISION3* vis_coords, IN PRECISION2* measured_visibilities,
		IN int* num_corrected_vis, IN int* maj_iter, OUT PRECISION2* output_visibilities, OUT PRECISION3* out_vis_coords, OUT Config* output_configs, OUT int* num_vis, OUT int* maj_iter_out){
	memcpy(out_vis_coords, vis_coords, NUM_VISIBILITIES * sizeof(PRECISION3));
	memcpy(output_visibilities, measured_visibilities, NUM_VISIBILITIES * sizeof(PRECISION2));

	for(int i = 0; i < SLICES; ++i){
		output_configs[i] = *config;
		maj_iter_out[i] = *maj_iter;
	}

	int vis_left = *num_corrected_vis;
	int ideal_vis_per_slice = NUM_VISIBILITIES/SLICES;

	for(int i = 0; i < SLICES; ++i){
		int curr_vis = vis_left > ideal_vis_per_slice ? ideal_vis_per_slice : vis_left;
		vis_left -= curr_vis;
		num_vis[i] = curr_vis;
	}
}

void split_uvgrid(int GRID_SIZE, int SLICES, IN PRECISION2* input_grid, OUT PRECISION2* output_grid){
	for(int i = 0; i < SLICES; ++i){
		unsigned int offset = i * GRID_SIZE * GRID_SIZE;
		PRECISION2* dest = output_grid + offset;
		memcpy(dest, input_grid, GRID_SIZE * GRID_SIZE * sizeof(PRECISION2));
	}
}
