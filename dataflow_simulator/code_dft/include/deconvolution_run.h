#ifndef DECONVOLUTION_RUN_H
#define DECONVOLUTION_RUN_H


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"

void deconvolution_run(int GRID_SIZE, int CALIBRATION, int NUMBER_MINOR_CYCLES_CAL, int NUMBER_MINOR_CYCLES_IMG, int NUM_MAX_SOURCES, int NUM_MAJOR_CYCLES,
			PRECISION* dirty_image_in, PRECISION* psf_in, PRECISION3* sources_in, int* num_sources_in, Config* config, int* num_sources_out, PRECISION3* sources_out, PRECISION* residual_image);

void hogbom_clean(int GRID_SIZE, int NUM_MINOR_CYCLES, int MAX_SOURCES, IN PRECISION* residual, IN PRECISION* partial_psf, int2* partial_psf_halfdims,
		IN Config* config, IN PRECISION* current_model, OUT int* num_sources_out, OUT PRECISION3* sources_out, OUT PRECISION* output_model);

void scale_dirty_image_by_psf_actor(int GRID_SIZE,
			IN PRECISION* dirty_image_in, IN PRECISION* psf, IN Config *config, OUT PRECISION* dirty_image_out);

void image_to_grid_coords_conversion_actor(int GRID_SIZE, int NUM_MAX_SOURCES,
			IN PRECISION3 *sources_in, IN int *num_sources, IN Config *config, OUT PRECISION3 *sources_out);

void reset_num_sources_actor(int CALIBRATION, IN int* num_sources_in, OUT int* num_sources_out);

void loop_iterator(int NUMBER_MINOR_CYCLES, OUT int *cycle_count);

void find_max_source_row_reduction_actor(int GRID_SIZE,
			IN PRECISION *dirty_image_in, IN int *loop_token, IN Config *config, OUT PRECISION3 *max_locals);

void find_max_source_col_reduction_actor(int GRID_SIZE, int NUM_MAX_SOURCES, int CALIBRATION,
			IN PRECISION3 *sources_in, IN int *num_sources_in, IN PRECISION3 *max_locals, IN int *cycle_number, IN Config *config, OUT PRECISION3 *sources_out, OUT int *num_sources_out);

void grid_to_image_coords_conversion_actor(int GRID_SIZE, int NUM_MAX_SOURCES,
			IN PRECISION3 *sources_in, IN int *num_sources, IN Config *config, OUT PRECISION3 *sources_out);

void subtract_psf_from_image_actor(int GRID_SIZE, int NUM_MAX_SOURCES,
			IN PRECISION *dirty_image_in, IN PRECISION3 *sources, IN PRECISION* psf, IN int *cycle_number, IN int *num_sources, IN Config *config, OUT PRECISION *dirty_image_out);

void compress_sources_actor(int NUM_MAX_SOURCES, IN PRECISION3 *sources_in, IN int *num_sources_in, OUT PRECISION3 *sources_out, OUT int *num_sources_out);

#ifdef __NVCC__
	__global__ void scale_dirty_image_by_psf(PRECISION *image, PRECISION *psf, PRECISION psf_max, const int grid_size);

	__global__ void image_to_grid_coords_conversion(PRECISION3 *sources, const PRECISION cell_size, const int half_grid_size, const int source_count);

	__global__ void find_max_source_row_reduction(const PRECISION *image, PRECISION3 *local_max, const int grid_size);

	__global__ void find_max_source_col_reduction(PRECISION3 *sources, const PRECISION3 *local_max, const int cycle_number,
			const int grid_size, const double loop_gain, const double weak_source_percent, const double noise_factor);

	__global__ void grid_to_image_coords_conversion(PRECISION3 *sources, const PRECISION cell_size, const int half_grid_size,
			const int source_count);

	__global__ void subtract_psf_from_image(PRECISION *image, PRECISION3 *sources, const PRECISION *psf,
			const int cycle_number, const int grid_size, const PRECISION loop_gain);

	__global__ void compress_sources(PRECISION3 *sources);
#endif

#ifdef __cplusplus
}
#endif


#endif
