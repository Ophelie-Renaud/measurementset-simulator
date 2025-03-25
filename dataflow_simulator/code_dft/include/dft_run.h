#ifndef DFT_RUN_H
#define DFT_RUN_H

#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"
#include "common.h"


void dft_actor(int NUM_VISIBILITIES, int NUM_MAX_SOURCES,
			IN PRECISION3 *sources, IN PRECISION3 *vis_uvw_coords, IN int *num_sources, IN Config *config, OUT PRECISION2 *visibilities);

#ifdef __NVCC__

__global__ void direct_fourier_transform(const PRECISION3 *vis_uvw, PRECISION2 *predicted_vis,
		const int vis_count, const PRECISION3 *sources, const int source_count);

#endif

#ifdef __cplusplus
}
#endif


#endif