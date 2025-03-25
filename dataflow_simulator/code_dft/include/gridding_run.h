#ifndef GRIDDING_RUN_H
#define GRIDDING_RUN_H


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"


    void gridding_actor(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES,
                IN PRECISION2* kernels, IN int2* kernel_supports, IN PRECISION3* vis_uvw_coords, IN PRECISION2* visibilities, IN Config *config, OUT PRECISION2* uv_grid);

    void gridding_CPU(PRECISION2 *grid, const PRECISION2 *kernel, const int2 *supports,
                const PRECISION3 *vis_uvw, const PRECISION2 *vis, const int num_vis, const int oversampling,
                const int grid_size, const double uv_scale, const double w_scale);

    PRECISION2 complex_mult_CPU(const PRECISION2 z1, const PRECISION2 z2);

    void add_visibilities(int NUM_VISIBILITIES, IN PRECISION2* v1, IN PRECISION2* v2, OUT PRECISION2* output);


#ifdef __NVCC__
	__global__ void gridding(PRECISION2 *grid, const PRECISION2 *kernel, const int2 *supports,
    	        const PRECISION3 *vis_uvw, const PRECISION2 *vis, const int num_vis, const int oversampling,
    	        const int grid_size, const double uv_scale, const double w_scale);

    __device__ PRECISION2 complex_mult(const PRECISION2 z1, const PRECISION2 z2);
#endif




#ifdef __cplusplus
}
#endif


#endif
