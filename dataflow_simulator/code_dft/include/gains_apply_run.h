#ifndef GAINS_APPLY_RUN_H
#define GAINS_APPLY_RUN_H


#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"

//	 #include "common.h"
	// #include "timer.h"

    void reciprocal_transform_actor(int NUM_RECEIVERS, PRECISION2* gains_in, Config *config, PRECISION2* gains_out);

   // void reciprocal_transform_actor(int NUM_RECEIVERS,
	//			IN PRECISION2* gains_in, IN Config *config, OUT PRECISION2* gains_out);

    void apply_gains_actor(int CALIBRATION, int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES,
                           IN PRECISION2* measured_vis, IN PRECISION2* visibilities_in, IN PRECISION2* gains, IN int2* receiver_pairs, IN Config *config, OUT PRECISION2* visibilities_out);

    //this is exactly the same as the apply_gains_actor function minus the calibration variable
	void subtract_from_measurements(int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES,
					IN PRECISION2* measured_vis, IN PRECISION2* visibilities_in, IN PRECISION2* gains, IN int2* receiver_pairs, IN Config *config, OUT PRECISION2* visibilities_out);

	//applies the gains without subtracting. Used in the oversampled finegrid versions of visibilities as the gains need to be applied before
	//discretization of the measured visibilities, which are then in turn used for subtraction
	void apply_gains_only(int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES, IN PRECISION2* measured_vis, IN PRECISION2* gains,
			IN int2* receiver_pairs, IN Config *config, OUT PRECISION2* visibilities_out);

	void subtract_from_measurements_finegrid(int NUM_VISIBILITIES, IN PRECISION2* gridded_measured_vis, IN PRECISION2* input_vis, int* num_finegrid_vis, OUT PRECISION2* delta_vis);

#ifndef __NVCC__
    PRECISION2 complex_multiply_apply(const PRECISION2 z1, const PRECISION2 z2);
    PRECISION2 complex_conjugate_apply(const PRECISION2 z1);
    PRECISION2 complex_subtract_apply(const PRECISION2 z1, const PRECISION2 z2);
#endif

#ifdef __NVCC__
	__global__ void reciprocal_transform(PRECISION2 *gains_out, PRECISION2 *gains_in, const int num_recievers);

	__global__ void apply_gains(PRECISION2 *measured_vis, PRECISION2 *predicted_vis, const int num_vis,
				const PRECISION2 *gains_recip, const int2 *receiver_pairs, const int num_recievers, const int num_baselines);

	__global__ void apply_gains_subtraction(PRECISION2 *measured_vis, PRECISION2 *predicted_vis, const int num_vis,
				const PRECISION2 *gains_recip, const int2 *receiver_pairs, const int num_recievers, const int num_baselines);

	__device__ PRECISION2 complex_reciprocal_apply(const PRECISION2 z);

	__device__ PRECISION2 complex_conjugate_apply(const PRECISION2 z1);

	__device__ PRECISION2 complex_multiply_apply(const PRECISION2 z1, const PRECISION2 z2);

	__device__ PRECISION2 complex_divide(const PRECISION2 z1, const PRECISION2 z2);

	__device__ PRECISION2 complex_subtract_apply(const PRECISION2 z1, const PRECISION2 z2);
#endif



#ifdef __cplusplus
}
#endif


#endif
