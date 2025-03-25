#ifndef CONVOLUTION_CORRECTION_RUN_H
#define CONVOLUTION_CORRECTION_RUN_H


#ifdef __cplusplus
extern "C" {
#endif

#include "preesm.h"

void execute_convolution_correction_actor(int GRID_SIZE, IN PRECISION* dirty_image_in, IN PRECISION* prolate, IN Config *config, OUT PRECISION* dirty_image_out);
    void execute_inv_convolution_correction_actor(int GRID_SIZE, IN PRECISION* dirty_image_in, IN PRECISION* prolate, IN Config *config, OUT PRECISION* dirty_image_out);

#ifdef __NVCC__
__global__ void execute_convolution_correction(PRECISION *image, const PRECISION *prolate, const int image_size);
#endif

#ifdef __cplusplus
}
#endif


#endif
