#include "reduce_vis.h"

void reduce_vis(int NUM_ACTUAL_VISIBILITIES, int NUM_VISIBILITIES, PRECISION2* vis_in, PRECISION3* vis_coords_in, PRECISION2* vis_out, PRECISION3* vis_coords_out){
	memcpy(vis_out, vis_in, sizeof(PRECISION2) * NUM_VISIBILITIES);
	memcpy(vis_coords_out, vis_coords_in, sizeof(PRECISION3) * NUM_VISIBILITIES);
}
