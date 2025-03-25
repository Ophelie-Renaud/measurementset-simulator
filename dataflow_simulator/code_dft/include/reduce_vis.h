#ifndef REDUCEVIS_H
#define REDUCEVIS_H


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"
	void reduce_vis(int NUM_ACTUAL_VISIBILITIES, int NUM_VISIBILITIES, PRECISION2* vis_in, PRECISION3* vis_coords_in, PRECISION2* vis_out, PRECISION3* vis_coords_out);
#ifdef __cplusplus
}
#endif


#endif
