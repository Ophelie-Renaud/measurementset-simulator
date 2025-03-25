#ifndef SPLITMERGE_H
#define SPLITMERGE_H

#ifndef CPU_VERSION
#define CPU_VERSION 1
#endif


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"

void split_visibilities(IN Config* config, IN int NUM_VISIBILITIES, IN int SLICES, IN PRECISION2* measured_visibilities, IN PRECISION3* vis_coords,
		OUT PRECISION2* visibilities, OUT PRECISION3* out_vis_coords, OUT Config* output_configs, OUT Config* output_config);
void merge_gridded_visibilities(IN int GRID_SIZE, IN int SLICES, IN PRECISION2* grids, OUT PRECISION2* output_grid);
void merge_configs(IN Config* in_configs, OUT Config* out_config);
void split_fftdegrid_visibilities(int NUM_VISIBILITIES, int SLICES, int GRID_SIZE, IN Config* config, IN PRECISION3* vis_coords, IN PRECISION2* measured_visibilities,
		IN int* num_corrected_vis, IN int* maj_iter, OUT PRECISION2* output_visibilities, OUT PRECISION3* out_vis_coords, OUT Config* output_configs, OUT int* num_vis, OUT int* maj_iter_out);
void split_uvgrid(int GRID_SIZE, int SLICES, IN PRECISION2* input_grid, OUT PRECISION2* output_grid);

#ifdef __cplusplus
}
#endif

#endif
