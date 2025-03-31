#ifndef TOP_H
#define TOP_H

#ifndef CPU_VERSION
#define CPU_VERSION 1
#endif


#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"

	void end_sink(int NUM_RECEIVERS, IN PRECISION2 *gains);

	void config_struct_set_up(int GRID_SIZE, int NUM_KERNELS, OUT Config *config);
	void config_struct_set_up_v2(int GRID_SIZE, int NUM_KERNELS, int NUM_BASELINES, int TOTAL_KERNEL_SAMPLES, int OVERSAMPLING_FACTOR,OUT Config *config_struct) ;

	void gains_host_set_up(int NUM_RECEIVERS, int NUM_BASELINES, IN Config *config, OUT PRECISION2 *gains, OUT int2 *receiver_pairs);

	//void calculate_receiver_pairs(Config *config, int2 *receiver_pairs);
    void calculate_receiver_pairs(int NUM_BASELINES, int NUM_RECEIVERS, int2 *receiver_pairs);

	void visibility_host_set_up(int NUM_VISIBILITIES, IN Config *config, OUT PRECISION3 *vis_uvw_coords, OUT PRECISION2 *measured_vis);

	void kernel_host_set_up(int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, IN Config *config, OUT int2 *kernel_supports, OUT PRECISION2 *kernels);

	void degridding_kernel_host_set_up(int NUM_DEGRIDDING_KERNELS, int TOTAL_DEGRIDDING_KERNEL_SAMPLES, IN Config *config,
			OUT int2 *degridding_kernel_supports, OUT PRECISION2 *degridding_kernels);

	void correction_set_up(int GRID_SIZE, OUT PRECISION *prolate);

	void create_1D_half_prolate(PRECISION *prolate, int grid_size);

	double prolate_spheroidal(double nu);

	void psf_host_set_up(int GRID_SIZE, int PSF_GRID_SIZE, IN Config *config, OUT PRECISION *psf, OUT double *psf_max_value);

	void config_struct_set_up_sequel(IN Config *config_in, IN double *psf_max_value, OUT Config *config_out);

	void clean_psf_host_set_up(int GRID_SIZE, int GAUSSIAN_CLEAN_PSF, IN Config *config, IN PRECISION *dirty_psf, OUT PRECISION *clean_psf, OUT int2* partial_psf_halfdims);

	void model_set_up(int GRID_SIZE, IN Config *config, OUT PRECISION *initial_image);

	void sources_set_up(int MAX_SOURCES, OUT int *num_sources, OUT PRECISION3 *sources);

	void visibility_sink(PRECISION2* delta_vis);

	//TESTING FUNCTIONS
	void initsink(Config* config_delta, Config* config_psi, Config* config_save, PRECISION* prolate, PRECISION* psf_delta, PRECISION* psf_psi, PRECISION* psf_save,
				PRECISION* psf_clean, int2* receiver_pairs, PRECISION2* gains, int2* kernel_supports, PRECISION2* kernels, PRECISION2* measured_vis,
				PRECISION3* vis_uwv_coords, int2* partial_psf_halfdims, int* cycle);

	void initdeltasink(int GRID_SIZE, Config* config_psi, Config* config_save, PRECISION* psf_psi, PRECISION* psf_save, PRECISION* clean_psf, int2* partial_psf_halfdims
			, PRECISION* delta_image, PRECISION* image_out
			);

	void init_sourcelist(int MAX_SOURCES, int* num_sources, PRECISION3* source_list);

	void init_image(int GRID_SIZE, PRECISION* image);

	void delta_pseudo(int NUM_VISIBILITIES, int NUM_RECEIVERS, int GRID_SIZE, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES, int NUM_BASELINES, int MAX_SOURCES,
			Config* config, int* cycle, PRECISION2* delta_vis, PRECISION2* gains, PRECISION* image_estimate, int2* kernel_supports,
			PRECISION2* kernels, int* num_sources_in, PRECISION* prolate, PRECISION* psf, int2* receiver_pairs, PRECISION3* source_list, PRECISION3* vis_coords,
			PRECISION* delta_image, PRECISION2* delta_vis_out, PRECISION2* gains_out, PRECISION* image_out
			);

	void extra_delta_sinks(int NUM_VISIBILITIES, int NUM_RECEIVERS, PRECISION2* delta_vis_out, PRECISION2* gains);

	void cycle_num(OUT int* cycle);

	void additional_sinks(PRECISION3* source_list, int* num_sources, PRECISION* image_sink);

	void pseudo_psi(Config* config, int* cycle, PRECISION* delta_image, PRECISION* input_model, int2* partial_psf_halfdims, PRECISION* psf,
			PRECISION* image_estimate, int* num_sources_out, PRECISION3* source_list);

#ifdef __cplusplus
}
#endif


#endif
