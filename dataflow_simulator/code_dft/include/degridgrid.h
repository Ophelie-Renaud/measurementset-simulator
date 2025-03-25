#ifndef DEGRIDGRID_H
#define DEGRIDGRID_H


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"
	void correct_to_finegrid(int NUM_VISIBILITIES, int GRID_SIZE, int OVERSAMPLING_FACTOR, int PERFORM_SIMPLIFICATION, PRECISION3* vis_uvw_coords, PRECISION2* input_visibilities,
			Config* config, PRECISION2* output_visibilities, PRECISION3* output_finegrid_vis_coords, int* num_output_visibilities);

	void std_gridding(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR, int BYPASS, int* maj_iter, int* num_corrected_visibilities,
			PRECISION2* kernels, int2* kernel_supports, PRECISION3* corrected_vis_uvw_coords, PRECISION2* visibilities, Config* config, PRECISION2* prev_grid, PRECISION2* output_grid);

    void std_degridding(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR, PRECISION2* kernels,
    		int2* kernel_supports, PRECISION2* input_grid, PRECISION3* corrected_vis_uvw_coords, int* num_corrected_visibilities, Config* config,
    		PRECISION2* output_visibilities);

    void g2g_degridgrid(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_GRIDDING_KERNELS, int NUM_DEGRIDDING_KERNELS,
			int TOTAL_GRIDDING_KERNEL_SAMPLES, int TOTAL_DEGRIDDING_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR,
    		PRECISION2* gridding_kernels, int2* gridding_kernel_supports, PRECISION2* degridding_kernels, int2* degridding_kernel_supports, PRECISION2* input_grid,
			PRECISION3* corrected_vis_uvw_coords, int* num_corrected_visibilities, Config* config, PRECISION2* output_grid);

    void subtract_image_space(int GRID_SIZE, PRECISION* measurements, PRECISION* estimate, PRECISION* result);

    void degridding_kernel_sink(int NUM_DEGRIDDING_KERNELS, int TOTAL_DEGRIDDING_KERNEL_SAMPLES, int2* supports, PRECISION2* kernels);

#ifdef __cplusplus
}
#endif


#endif
