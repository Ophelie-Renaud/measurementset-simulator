#ifndef GAINS_CALIBRATION_RUN_H
#define GAINS_CALIBRATION_RUN_H

#ifdef __NVCC__
#include <cusolverDn.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

   // #include <iostream>
    #include <stdio.h>
    #include <lapacke.h>

    #include "preesm.h"

	 #include "common.h"
	// #include "timer.h"

	#define CUDA_SOLVER_CHECK_RETURN(value) check_cuda_solver_error_aux(__FILE__,__LINE__, #value, value)


	void update_gain_calibration_actor(int NUM_BASELINES, int NUM_RECEIVERS, int NUM_VISIBILITIES,
				IN int *loop_cycle_token, IN PRECISION2 *measured_vis, IN PRECISION2 *visibilities, IN PRECISION2 *gains, IN int2 *receiver_pairs, IN Config *config, OUT PRECISION *A_array, OUT PRECISION *Q_array);

	void execute_calibration_SVD_actor(int NUM_COLS, int NUM_ROWS,
				IN PRECISION *A_array, OUT PRECISION *S_array, OUT PRECISION *U_array, OUT PRECISION *V_array);

	void calculate_SUQ_product_actor(int NUM_COLS, int NUM_ROWS,
				IN PRECISION *S_array, IN PRECISION *U_array, IN PRECISION *Q_array, IN Config *config, OUT PRECISION *SUQ_array);

	void calculate_delta_update_gains_actor(int NUM_COLS, int NUM_RECEIVERS,
				IN PRECISION *SUQ_array, IN PRECISION *V_array, IN Config *config, IN PRECISION2 *gains_in, OUT PRECISION2* gains_out);



    void update_gain_calibration_CPU(const PRECISION2 *vis_measured_array, const PRECISION2 *vis_predicted_array,
                                     const PRECISION2 *gains_array, const int2 *receiver_pairs, PRECISION *A_array, PRECISION *Q_array, const int num_recievers, const int num_baselines);

    PRECISION2 complex_multiply_CPU(const PRECISION2 z1, const PRECISION2 z2);

    PRECISION2 flip_for_i_CPU(const PRECISION2 z);

    PRECISION2 flip_for_neg_i_CPU(const PRECISION2 z);

    PRECISION2 complex_conjugate_CPU(const PRECISION2 z1);

    PRECISION2 complex_subtract_CPU(const PRECISION2 z1, const PRECISION2 z2);


#ifdef __NVCC__

	#define CUDA_SOLVER_CHECK_RETURN(value) check_cuda_solver_error_aux(__FILE__,__LINE__, #value, value)

	void check_cuda_solver_error_aux(const char *file, unsigned line, const char *statement, cusolverStatus_t err);

	__global__ void update_gain_calibration(const PRECISION2 *vis_measured_array, const PRECISION2 *vis_predicted_array,
				const PRECISION2 *gains_array, const int2 *receiver_pairs, PRECISION *A_array, PRECISION *Q_array, const int num_recievers, const int num_baselines);

	__global__ void calculate_suq_product(const PRECISION *d_S, const PRECISION *d_U, const PRECISION *d_Q, PRECISION *d_SUQ, const int num_entries);

	__global__ void calculate_delta_update_gains(const PRECISION *d_V, const PRECISION *d_SUQ, PRECISION2 *d_gains, const int num_recievers, const int num_cols);

	__device__ PRECISION2 complex_multiply(const PRECISION2 z1, const PRECISION2 z2);

	__device__ PRECISION2 complex_conjugate(const PRECISION2 z1);

	__device__ PRECISION2 complex_subtract(const PRECISION2 z1, const PRECISION2 z2);

	__device__ PRECISION2 complex_reciprocal(const PRECISION2 z);

	__device__ PRECISION2 flip_for_i(const PRECISION2 z);

	__device__ PRECISION2 flip_for_neg_i(const PRECISION2 z);
#endif


#ifdef __cplusplus
}
#endif


#endif