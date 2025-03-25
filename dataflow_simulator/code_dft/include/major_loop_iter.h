#ifndef MAJOR_LOOP_ITER_H
#define MAJOR_LOOP_ITER_H


#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"

	void save_dirty_image_actor(int SAVE_DIRTY_IMAGE, int GRID_SIZE, IN PRECISION *dirty_image, IN int *cycle_count, IN Config *config);

	void save_estimated_gain_actor(int SAVE_ESTIMATED_GAIN, int NUM_RECEIVERS, IN PRECISION2 *estimated_gains, IN int *cycle_count, IN Config *config);

	void save_predicted_visibilities_actor(int SAVE_PREDICTED_VISIBILITIES, int NUM_VISIBILITIES, IN PRECISION2 *predicted_visibilities, IN PRECISION3 *vis_uvw_coords, IN int *cycle_count, IN Config *config);

	void save_residual_image_actor(int SAVE_RESIDUAL_IMAGE, int GRID_SIZE, IN PRECISION *residual_image, IN int *cycle_count, IN Config *config);

	void residual_image_sink(int SAVE_RESIDUAL_IMAGE, int GRID_SIZE, IN PRECISION *residual_image);

	void token_sink(IN int *token_in);

	void save_image_to_file(Config *config, PRECISION *image, const char *file_name, int start_x, int start_y, int range_x, int range_y, int cycle);

	void save_extracted_sources_actor(int GRID_SIZE, int SAVE_EXTRACTED_SOURCES, int NUM_MAX_SOURCES, IN int *num_sources, IN PRECISION3 *sources, IN int *cycle_count,
			IN Config *config, PRECISION* image_in, PRECISION* clean_psf, int2* clean_psf_halfdims, PRECISION* image_out);

	void save_extracted_sources(PRECISION3 *sources, int number_of_sources, const char *path, const char *output_file, int cycle);

	void save_output(int GRID_SIZE, IN PRECISION* residual, IN PRECISION* model, IN PRECISION* clean_psf, int2* clean_psf_halfdims, IN PRECISION* psf, IN Config *config, IN int *cycle);

	void source_list_sink(int MAX_SOURCES, IN PRECISION3* source_list, IN int* num_sources);

	void delta_visibility_sink(int NUM_VISIBILITIES, IN PRECISION2* visibilities);

	void psf_sink(int GRID_SIZE, IN PRECISION *psf);

	void pass_through_image(int GRID_SIZE, IN PRECISION* image, OUT PRECISION* output_image);

	void save_partial_psf(int GRID_SIZE, PRECISION* psf, int2* extents, int* cycle, Config* config);

#ifdef __cplusplus
}
#endif


#endif
