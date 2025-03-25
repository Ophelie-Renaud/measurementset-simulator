#ifndef TIMINGS_H
#define TIMINGS_H

#ifndef CPU_VERSION
#define CPU_VERSION 1
#endif


#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"
	#include "top.h"
	#include "dft_run.h"
	#include "gains_apply_run.h"
	#include "major_loop_iter.h"
	#include "gridding_run.h"
	#include "degridgrid.h"
	#include "fft_run.h"
	#include "convolution_correction_run.h"
	#include "deconvolution_run.h"

	void time_constant_setups(int NUM_SAMPLES);
	void time_gridsize_setups(int NUM_SAMPLES, int GRID_SIZE);
	void time_visibility_setups(int NUM_SAMPLES, int NUM_VISIBILITIES);
	void time_save_output(int NUM_SAMPLES, int GRID_SIZE);
	void time_dft(int NUM_SAMPLES, int NUM_MINOR_CYCLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES);
	void time_gains_application(int NUM_SAMPLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES);
	void time_add_visibilities(int NUM_SAMPLES, int NUM_VISIBILITIES);
	void time_prolate(int NUM_SAMPLES, int GRID_SIZE);
	void time_finegrid(int NUM_SAMPLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES);
	void time_subtract_ispace(int NUM_SAMPLES, int GRID_SIZE);
	void time_fftshift(int NUM_SAMPLES, int GRID_SIZE);
	void time_fft(int NUM_SAMPLES, int GRID_SIZE);
	void time_hogbom(int NUM_SAMPLES, int GRID_SIZE, int NUM_MINOR_CYCLES);

#ifdef __cplusplus
}
#endif


#endif
