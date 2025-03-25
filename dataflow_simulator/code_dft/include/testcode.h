#ifndef TESTCODE_H
#define TESTCODE_H


#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"

	void initsink(Config* config_delta, Config* config_psi, Config* config_save, PRECISION* prolate, PRECISION* psf_delta, PRECISION* psf_psi, PRECISION* psf_save,
			PRECISION* psf_clean, int2* receiver_pairs, PRECISION2* gains, int2* kernel_supports, PRECISION2* kernels, PRECISION2* measured_vis,
			PRECISION3* vis_uwv_coords, int2* partial_psf_halfdims);

#ifdef __cplusplus
}
#endif


#endif
