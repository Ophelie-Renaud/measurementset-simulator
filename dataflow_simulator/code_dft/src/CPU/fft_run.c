
#include "fft_run.h"


void fft_shift_complex_to_complex_actor(int GRID_SIZE, PRECISION2 *uv_grid_in, Config *config, PRECISION2 *uv_grid_out) {
    int row_index,col_index;

    printf("UPDATE >>> Shifting grid data for FFT...\n\n");
    // Perform 2D FFT shift
    for (row_index = 0;row_index < GRID_SIZE; row_index++)
	{
        for (col_index = 0; col_index < GRID_SIZE; col_index++)
        {
            int a = 1 - 2 * ((row_index + col_index) & 1);
            uv_grid_out[row_index * GRID_SIZE + col_index].x = uv_grid_in[row_index * GRID_SIZE + col_index].x * a;
            uv_grid_out[row_index * GRID_SIZE + col_index].y = uv_grid_in[row_index * GRID_SIZE + col_index].y * a;
        }
	}
	MD5_Update(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE, uv_grid_out);
}


void CUFFT_EXECUTE_INVERSE_C2C_actor(int GRID_SIZE, PRECISION2 *uv_grid_in, PRECISION2 *uv_grid_out)
{
    printf("UPDATE >>> Performing iFFT...\n\n");

    //memcpy(uv_grid_out, uv_grid_in, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
    //return;

    //memset(uv_grid_out, 0, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

	#if SINGLE_PRECISION

    // Need to link with -lfftw3f instead of or in addition to -lfftw3
	fftwf_plan fft_plan;

    fft_plan = fftwf_plan_dft_2d(GRID_SIZE, GRID_SIZE, (float (*)[2]) uv_grid_in, (float (*)[2]) uv_grid_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(fft_plan);
    fftwf_destroy_plan(fft_plan);

	#else

    fftw_plan fft_plan;

    fft_plan = fftw_plan_dft_2d(GRID_SIZE, GRID_SIZE, (double (*)[2]) uv_grid_in, (double (*)[2]) uv_grid_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);
	#endif

    //memcpy(uv_grid_out, uv_grid_in, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

    MD5_Update(sizeof(PRECISION2) * GRID_SIZE*GRID_SIZE, uv_grid_out);
}

void CUFFT_EXECUTE_FORWARD_C2C_actor(int GRID_SIZE, PRECISION2 *uv_grid_in, PRECISION2 *uv_grid_out)
{
    printf("UPDATE >>> Performing FFT...\n\n");

    //memcpy(uv_grid_out, uv_grid_in, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
    //return;

	#if SINGLE_PRECISION

    // Need to link with -lfftw3f instead of or in addition to -lfftw3
	fftwf_plan fft_plan;

    fft_plan = fftwf_plan_dft_2d(GRID_SIZE, GRID_SIZE, (float (*)[2]) uv_grid_in, (float (*)[2]) uv_grid_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(fft_plan);
    fftwf_destroy_plan(fft_plan);

	#else

    fftw_plan fft_plan;

    fft_plan = fftw_plan_dft_2d(GRID_SIZE, GRID_SIZE, (double (*)[2]) uv_grid_in, (double (*)[2]) uv_grid_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);
	#endif

    //memcpy(uv_grid_out, uv_grid_in, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

    MD5_Update(sizeof(PRECISION2) * GRID_SIZE*GRID_SIZE, uv_grid_out);
}


void fft_shift_complex_to_real_actor(int GRID_SIZE, PRECISION2 *uv_grid, Config *config, PRECISION *dirty_image) {
    int grid_square = GRID_SIZE * GRID_SIZE;
    int row_index,col_index;

    printf("UPDATE >>> Shifting grid data for FFT...\n\n");
    // Perform 2D FFT shift back
    for (row_index = 0; row_index < GRID_SIZE; row_index++)
	{
        for (col_index = 0; col_index < GRID_SIZE; col_index++)
        {
            int a = 1 - 2 * ((row_index + col_index) & 1);
            dirty_image[row_index * GRID_SIZE + col_index] = uv_grid[row_index * GRID_SIZE + col_index].x * a;
        }
	}

    MD5_Update(sizeof(PRECISION) * grid_square, dirty_image);
}

void fft_shift_real_to_complex_actor(int GRID_SIZE, PRECISION *image, Config *config, PRECISION2 *fourier) {
    int grid_square = GRID_SIZE * GRID_SIZE;
    int row_index,col_index;

    printf("UPDATE >>> Shifting grid data for FFT...\n\n");
    // Perform 2D FFT shift back
    for (row_index = 0; row_index < GRID_SIZE; row_index++)
	{
        for (col_index = 0; col_index < GRID_SIZE; col_index++)
        {
            int a = 1 - 2 * ((row_index + col_index) & 1);
            fourier[row_index * GRID_SIZE + col_index].x = image[row_index * GRID_SIZE + col_index] * a;
            fourier[row_index * GRID_SIZE + col_index].y = 0;
        }
	}

    MD5_Update(sizeof(PRECISION2) * grid_square, fourier);
}

