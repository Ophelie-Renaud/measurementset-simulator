#ifndef FFT_RUN_H
#define FFT_RUN_H




#ifdef __cplusplus
extern "C" {
#endif
	#include "preesm.h"
    #include <fftw3.h>

	void fft_shift_complex_to_complex_actor(int GRID_SIZE, IN PRECISION2* uv_grid_in, IN Config *config, OUT PRECISION2* uv_grid_out);

	void CUFFT_EXECUTE_FORWARD_C2C_actor(int GRID_SIZE, IN PRECISION2* uv_grid_in, OUT PRECISION2* uv_grid_out);

	void CUFFT_EXECUTE_INVERSE_C2C_actor(int GRID_SIZE, IN PRECISION2* uv_grid_in, OUT PRECISION2* uv_grid_out);

	void fft_shift_complex_to_real_actor(int GRID_SIZE, IN PRECISION2* uv_grid, IN Config *config, OUT PRECISION* dirty_image);

	void fft_shift_real_to_complex_actor(int GRID_SIZE, IN PRECISION *image, IN Config *config, OUT PRECISION2 *fourier);


#ifdef __NVCC__
	__global__ void fft_shift_complex_to_complex(PRECISION2 *grid, const int width);

	__global__ void fft_shift_complex_to_real(PRECISION2 *grid, PRECISION *image, const int width);
#endif


#ifdef __cplusplus
}
#endif


#endif
