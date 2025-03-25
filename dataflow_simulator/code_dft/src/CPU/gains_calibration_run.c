
#include "gains_calibration_run.h"



//#if SINGLE_PRECISION
//extern void sgesvd_(char* jobu, char* jobvt, int* m, int* n, float* a,	int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork, int* info);
//#else
//extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a,	int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info);
//#endif









void update_gain_calibration_actor(int NUM_BASELINES, int NUM_RECEIVERS, int NUM_VISIBILITIES,
			int *loop_cycle_token, PRECISION2 *measured_vis, PRECISION2 *visibilities, PRECISION2 *gains, int2 *receiver_pairs, Config *config, PRECISION *A_array, PRECISION *Q_array)
{

//    printf("gains MD5 \t: ");
//    MD5_Update(sizeof(PRECISION2) * NUM_RECEIVERS, gains);

//    printf("measured_vis MD5 \t: ");
//    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, measured_vis);
//
//    printf("visibilities MD5 \t: ");
//    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities);
//
//    printf("receiver_pairs MD5 \t: ");
//    MD5_Update(sizeof(int2) * NUM_BASELINES, receiver_pairs);

    memset(Q_array, 0, sizeof(PRECISION) * 2 * NUM_RECEIVERS);
    memset(A_array, 0, sizeof(PRECISION) * 2 * NUM_RECEIVERS * 2 * NUM_RECEIVERS);

	printf("UPDATE >>> Performing Gain Calibration cycle: %d ...\n", *loop_cycle_token);
		//EXECUTE CUDA KERNEL UPDATING Q AND A array (NEED ATOMIC ACCESS!)
	update_gain_calibration_CPU(
			measured_vis,
			visibilities,
			gains,
			receiver_pairs,
			A_array,
			Q_array,
			NUM_RECEIVERS,
			NUM_BASELINES
	);
	// cudaDeviceSynchronize();

//    printf("A_array MD5 \t: ");
//    MD5_Update(sizeof(PRECISION) * 2 * NUM_RECEIVERS * 2 * NUM_RECEIVERS, A_array);
//    printf("Q_array MD5 \t: ");
//    MD5_Update(sizeof(PRECISION) * 2 * NUM_RECEIVERS, Q_array);

}



void update_gain_calibration_CPU(const PRECISION2 *vis_measured_array, const PRECISION2 *vis_predicted_array,
                                 const PRECISION2 *gains_array, const int2 *receiver_pairs, PRECISION *A_array, PRECISION *Q_array,
                                 const int num_recievers, const int num_baselines)
{
    // const int index = threadIdx.x + blockDim.x * blockIdx.x;
    // if(index >= num_baselines)
    // 	return;

    for (int index = 0; index < num_baselines; index++)
    {

        PRECISION2 vis_measured = vis_measured_array[index];
        PRECISION2 vis_predicted = vis_predicted_array[index];

        int2 antennae = receiver_pairs[index];

        PRECISION2 gainA = gains_array[antennae.x];
        //only need gainB as conjugate??
        PRECISION2 gainB_conjugate = complex_conjugate_CPU(gains_array[antennae.y]);

        //NOTE do not treat residual as a COMPLEX!!!!! (2 reals)
        PRECISION2 residual = complex_subtract_CPU(vis_measured, complex_multiply_CPU(vis_predicted,complex_multiply_CPU(gainA, gainB_conjugate)));


        //CALCULATE Partial Derivatives

        PRECISION2 part_respect_to_real_gain_a = complex_multiply_CPU(vis_predicted, gainB_conjugate);

        PRECISION2 part_respect_to_imag_gain_a = flip_for_i_CPU(complex_multiply_CPU(vis_predicted, gainB_conjugate));

        PRECISION2 part_respect_to_real_gain_b = complex_multiply_CPU(vis_predicted,gainA);

        PRECISION2 part_respect_to_imag_gain_b = flip_for_neg_i_CPU(complex_multiply_CPU(vis_predicted, gainA));

        //Calculate Q[2a],Q[2a+1],Q[2b],Q[2b+1] arrays - In this order... NEED ATOMIC UPDATE
        double qValue = part_respect_to_real_gain_a.x * residual.x
                        + part_respect_to_real_gain_a.y * residual.y;
        // atomicAdd(&(Q_array[2*antennae.x]), qValue);
        Q_array[2*antennae.x] += qValue;

        qValue = part_respect_to_imag_gain_a.x * residual.x
                 + part_respect_to_imag_gain_a.y * residual.y;
        // atomicAdd(&(Q_array[2*antennae.x+1]), qValue);
        Q_array[2*antennae.x+1] += qValue;

        qValue = part_respect_to_real_gain_b.x * residual.x
                 + part_respect_to_real_gain_b.y * residual.y;
        // atomicAdd(&(Q_array[2*antennae.y]), qValue);
        Q_array[2*antennae.y] += qValue;

        qValue = part_respect_to_imag_gain_b.x * residual.x
                 + part_respect_to_imag_gain_b.y * residual.y;
        // atomicAdd(&(Q_array[2*antennae.y+1]), qValue);
        Q_array[2*antennae.y+1] += qValue;


        int num_cols = 2 * num_recievers;
        //CALCULATE JAcobian product on A matrix... 2a2a, 2a2a+1, 2a2b, 2a2b+1
        //2a2a
        double aValue = part_respect_to_real_gain_a.x * part_respect_to_real_gain_a.x +
                        part_respect_to_real_gain_a.y * part_respect_to_real_gain_a.y;

        int aIndex = (2 *  antennae.x * num_cols) + (2 * antennae.x);

        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a2a+1,
        aValue = part_respect_to_real_gain_a.x * part_respect_to_imag_gain_a.x +
                 part_respect_to_real_gain_a.y * part_respect_to_imag_gain_a.y;

        aIndex = (2 *  antennae.x * num_cols) + (2 * antennae.x + 1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a2b
        aValue = part_respect_to_real_gain_a.x * part_respect_to_real_gain_b.x +
                 part_respect_to_real_gain_a.y * part_respect_to_real_gain_b.y;

        aIndex = (2 *  antennae.x * num_cols) + (2 * antennae.y);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a2b+1
        aValue = part_respect_to_real_gain_a.x * part_respect_to_imag_gain_b.x +
                 part_respect_to_real_gain_a.y * part_respect_to_imag_gain_b.y;

        aIndex = (2 *  antennae.x * num_cols) + (2 * antennae.y+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //CACLUATE JAcobian product on A matrix... [2a+1,2a], [2a+1,2a+1], [2a+1,2b], [2a+1,2b+1]
        //2a+1,2a
        aValue = part_respect_to_imag_gain_a.x * part_respect_to_real_gain_a.x +
                 part_respect_to_imag_gain_a.y * part_respect_to_real_gain_a.y;

        aIndex = ((2 *  antennae.x+1) * num_cols) + (2 * antennae.x);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a+1,2a+1
        aValue = part_respect_to_imag_gain_a.x * part_respect_to_imag_gain_a.x +
                 part_respect_to_imag_gain_a.y * part_respect_to_imag_gain_a.y;

        aIndex = ((2 *  antennae.x+1) * num_cols) + (2 * antennae.x+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a+1,2b
        aValue = part_respect_to_imag_gain_a.x * part_respect_to_real_gain_b.x +
                 part_respect_to_imag_gain_a.y * part_respect_to_real_gain_b.y;

        aIndex = ((2 *  antennae.x+1) * num_cols) + (2 * antennae.y);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2a+1,2b+1
        aValue = part_respect_to_imag_gain_a.x * part_respect_to_imag_gain_b.x +
                 part_respect_to_imag_gain_a.y * part_respect_to_imag_gain_b.y;

        aIndex = ((2 *  antennae.x+1) * num_cols) + (2 * antennae.y+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //CACLUATE JAcobian product on A matrix... [2b,2a], [2b,2a+1], [2b,2b], [2b,2b+1]
        //2b,2a
        aValue = part_respect_to_real_gain_b.x * part_respect_to_real_gain_a.x +
                 part_respect_to_real_gain_b.y * part_respect_to_real_gain_a.y;

        aIndex = (2 *  antennae.y * num_cols) + (2 * antennae.x);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b,2a+1
        aValue = part_respect_to_real_gain_b.x * part_respect_to_imag_gain_a.x +
                 part_respect_to_real_gain_b.y * part_respect_to_imag_gain_a.y;

        aIndex = (2 *  antennae.y * num_cols) + (2 * antennae.x+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b,2b
        aValue = part_respect_to_real_gain_b.x * part_respect_to_real_gain_b.x +
                 part_respect_to_real_gain_b.y * part_respect_to_real_gain_b.y;

        aIndex = (2 *  antennae.y * num_cols) + (2 * antennae.y);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b, 2b+1
        aValue = part_respect_to_real_gain_b.x * part_respect_to_imag_gain_b.x +
                 part_respect_to_real_gain_b.y * part_respect_to_imag_gain_b.y;

        aIndex = (2 *  antennae.y * num_cols) + (2 * antennae.y + 1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //CALCULATE JAcobian product on A matrix... [2b+1,2a], [2b+1,2a+1], [2b+1,2b], [2b+1,2b+1]
        //2b+1,2a
        aValue =  part_respect_to_imag_gain_b.x * part_respect_to_real_gain_a.x +
                  part_respect_to_imag_gain_b.y * part_respect_to_real_gain_a.y;

        aIndex = ((2 *  antennae.y+1) * num_cols) + (2 * antennae.x);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b+1,2a+1
        aValue =  part_respect_to_imag_gain_b.x * part_respect_to_imag_gain_a.x+
                  part_respect_to_imag_gain_b.y * part_respect_to_imag_gain_a.y;

        aIndex = ((2 *  antennae.y+1) * num_cols) + (2 * antennae.x+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b+1,2b
        aValue =  part_respect_to_imag_gain_b.x * part_respect_to_real_gain_b.x+
                  part_respect_to_imag_gain_b.y * part_respect_to_real_gain_b.y;

        aIndex = ((2 *  antennae.y+1) * num_cols) + (2 * antennae.y);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;

        //2b+1, 2b+1
        aValue = part_respect_to_imag_gain_b.x * part_respect_to_imag_gain_b.x +
                 part_respect_to_imag_gain_b.y * part_respect_to_imag_gain_b.y;

        aIndex = ((2 *  antennae.y+1) * num_cols) + (2 * antennae.y+1);
        // atomicAdd(&(A_array[aIndex]), aValue);
        A_array[aIndex] += aValue;
    }
}

void execute_calibration_SVD_actor(int NUM_COLS, int NUM_ROWS,
                                   PRECISION *A_array, PRECISION *S_array, PRECISION *U_array, PRECISION *V_array)
{
    int lwork,info;
    PRECISION wkopt;
//	double wkopt;

//    double* a_local = (double*) malloc(NUM_COLS*NUM_ROWS * sizeof(double));
//	double* s_local = (double*) malloc(NUM_COLS * sizeof(double));
//	double* u_local = (double*) malloc(NUM_COLS*NUM_ROWS * sizeof(double));
//	double* v_local = (double*) malloc(NUM_COLS*NUM_COLS * sizeof(double));
//
//	for (int i = 0; i < NUM_COLS*NUM_ROWS; i++)
//		a_local[i] = (double) A_array[i];

	/* Query and allocate the optimal workspace */
    lwork = -1;
	#if SINGLE_PRECISION
		LAPACK_sgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, A_array, &NUM_ROWS, S_array, U_array, &NUM_ROWS, V_array, &NUM_COLS, &wkopt, &lwork, &info);
//		LAPACK_dgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, a_local, &NUM_ROWS, s_local, u_local, &NUM_ROWS, v_local, &NUM_COLS, &wkopt, &lwork, &info);
	#else
		LAPACK_dgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, A_array, &NUM_ROWS, S_array, U_array, &NUM_ROWS, V_array, &NUM_COLS, &wkopt, &lwork, &info);
	#endif

    lwork = (int)wkopt;

    PRECISION* work = (PRECISION*) malloc(lwork * sizeof(PRECISION));
    memset(work, 0, sizeof(PRECISION) * lwork);

//	double* work = (double*) malloc(lwork * sizeof(double));
//	memset(work, 0, sizeof(double) * lwork);

    info = 0;

    #if SINGLE_PRECISION
		LAPACK_sgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, A_array, &NUM_ROWS, S_array, U_array, &NUM_ROWS, V_array, &NUM_COLS, work, &lwork, &info );
//		LAPACK_dgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, a_local, &NUM_ROWS, s_local, u_local, &NUM_ROWS, v_local, &NUM_COLS, work, &lwork, &info );
    #else
        LAPACK_dgesvd( "All", "All", &NUM_ROWS, &NUM_COLS, A_array, &NUM_ROWS, S_array, U_array, &NUM_ROWS, V_array, &NUM_COLS, work, &lwork, &info );
	#endif


//	for (int i = 0; i < NUM_COLS; i++)
//		S_array[i] = (float) s_local[i];
//
//	for (int i = 0; i < NUM_COLS*NUM_ROWS; i++)
//		U_array[i] = (float) u_local[i];
//
//	for (int i = 0; i < NUM_COLS*NUM_COLS; i++)
//		V_array[i] = (float) v_local[i];
//
//	free(a_local);
//	free(s_local);
//	free(u_local);
//	free(v_local);

    /* Check for convergence */
    if( info > 0 ) {
        printf("The algorithm computing SVD failed to converge.\n");
        exit(1);
    }

    free(work);
}









//void execute_calibration_SVD_actor(int NUM_COLS, int NUM_ROWS, PRECISION *A_array, PRECISION *S_array, PRECISION *U_array, PRECISION *V_array)
//{
////    PRECISION *d_A_array = NULL;
////    PRECISION *d_U_array = NULL;
////    PRECISION *d_V_array = NULL;
////    PRECISION *d_S_array = NULL;
//
////    PRECISION *work;
////
////    cusolverDnHandle_t solver_handle;
////    cusolverDnCreate(&solver_handle);
////
////    int work_size = 0;
////    CUDA_SOLVER_CHECK_RETURN(SVD_BUFFER_SIZE(solver_handle, NUM_ROWS, NUM_COLS, &work_size));
////   // LAPACKE_dgesvd
////
////    CUDA_CHECK_RETURN(cudaMalloc(&work, work_size * sizeof(PRECISION)));
////    CUDA_CHECK_RETURN(cudaMemset(work, 0, sizeof(PRECISION) * work_size));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMalloc(&d_A_array,  sizeof(PRECISION) * NUM_ROWS * NUM_COLS));
////    CUDA_CHECK_RETURN(cudaMemcpy(d_A_array, A_array, sizeof(PRECISION) * NUM_ROWS * NUM_COLS, cudaMemcpyHostToDevice));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMalloc(&d_S_array, sizeof(PRECISION) * NUM_ROWS));
////    CUDA_CHECK_RETURN(cudaMemset(d_S_array, 0, sizeof(PRECISION) * NUM_ROWS));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMalloc(&d_U_array, sizeof(PRECISION) * NUM_ROWS * NUM_COLS));
////    CUDA_CHECK_RETURN(cudaMemset(d_U_array, 0, sizeof(PRECISION) * NUM_ROWS * NUM_COLS));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMalloc(&d_V_array, sizeof(PRECISION) * NUM_ROWS * NUM_COLS));
////    CUDA_CHECK_RETURN(cudaMemset(d_V_array, 0, sizeof(PRECISION) * NUM_ROWS * NUM_COLS));
////    cudaDeviceSynchronize();
////
////
////    int *devInfo;
////    CUDA_CHECK_RETURN(cudaMalloc(&devInfo, sizeof(int)));
//
//    //ALL OUR MATRICES ARE "ROW" MAJOR HOWEVER AS A IS SYMMETRIC DOES NOT NEED TO BE TRANSPOSED FOR FOR SVD ROUTINE
//    //SO NEED TO NOT TRANSPOSE U AND TRANSPOSE VSTAR
////    CUDA_SOLVER_CHECK_RETURN(SVD(solver_handle, 'A', 'A', NUM_ROWS, NUM_COLS, d_A_array, NUM_ROWS, d_S_array, d_U_array, NUM_ROWS, d_V_array, NUM_COLS, work, work_size, NULL, devInfo));
//
//    int m = NUM_ROWS, n = NUM_COLS, lda = NUM_ROWS, ldu = NUM_ROWS, ldvt = NUM_COLS, info, lwork;
//    double wkopt;
//    double* work;
//
//    lwork = -1;
//    dgesvd( "All", "All", &m, &n, A_array, &lda, S_array, U_array, &ldu, V_array, &ldvt, &wkopt, &lwork, &info );
//    lwork = (int)wkopt;
//    work = (double*)malloc( lwork*sizeof(double) );
//    dgesvd( "All", "All", &m, &n, A_array, &lda, S_array, U_array, &ldu, V_array, &ldvt, work, &lwork, &info );
//
//
//
//
//
//
//    /* Check for convergence */
//    if( info > 0 ) {
//        printf( "The algorithm computing SVD failed to converge.\n" );
//        exit( 1 );
//    }
////    cudaDeviceSynchronize();
//
////    int devInfo_h = 0;
////    CUDA_CHECK_RETURN(cudaMemcpy(&devInfo_h, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
////    cudaDeviceSynchronize();
////
////
////
////    CUDA_CHECK_RETURN(cudaMemcpy(S_array, d_S_array, sizeof(PRECISION) * NUM_ROWS, cudaMemcpyDeviceToHost));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMemcpy(U_array, d_U_array, sizeof(PRECISION) * NUM_ROWS * NUM_COLS, cudaMemcpyDeviceToHost));
////    cudaDeviceSynchronize();
////
////    CUDA_CHECK_RETURN(cudaMemcpy(V_array, d_V_array, sizeof(PRECISION) * NUM_ROWS * NUM_COLS, cudaMemcpyDeviceToHost));
////    cudaDeviceSynchronize();
////
////
////    CUDA_CHECK_RETURN(cudaFree(d_S_array));
////    CUDA_CHECK_RETURN(cudaFree(d_U_array));
////    CUDA_CHECK_RETURN(cudaFree(d_V_array));
////    CUDA_CHECK_RETURN(cudaFree(d_A_array));
////    CUDA_CHECK_RETURN(cudaFree(work));
//
//
////    printf("h_V_array MD5 \t: ");
////    MD5_Update(sizeof(PRECISION) * NUM_ROWS * NUM_COLS, V_array);
////
////    printf("h_U_array MD5 \t: ");
////    MD5_Update(sizeof(PRECISION) * NUM_ROWS * NUM_COLS, U_array);
////
////    printf("h_S_array MD5 \t: ");
////    MD5_Update(sizeof(PRECISION) * NUM_ROWS, S_array);
//
//
//
//    //CHECKING S PRODUCT!!
//    //bool success = (devInfo_h == 0);
//    //printf("UPDATE >>> SVD complete...\n\n");
////    if (devInfo) CUDA_CHECK_RETURN(cudaFree(devInfo));
//    //return success;
//}






void calculate_SUQ_product_actor(int NUM_COLS, int NUM_ROWS,
                                 PRECISION *S_array, PRECISION *U_array, PRECISION *Q_array, Config *config, PRECISION *SUQ_array)
{
    int index;

    for (index =0; index < NUM_COLS; index ++) {
        PRECISION sinv = (ABS(S_array[index]) > 1E-6) ? 1.0 / S_array[index] : 0.0;


        PRECISION product = 0;
        for (int i = 0; i < NUM_COLS; ++i) {
            product += Q_array[i] * U_array[index * NUM_COLS + i];
        }

        SUQ_array[index] = product * sinv;
    }
}

void calculate_delta_update_gains_actor(int NUM_COLS, int NUM_RECEIVERS,
                                        PRECISION *SUQ_array, PRECISION *V_array, Config *config, PRECISION2 *gains_in, PRECISION2* gains_out)
{
    int index;

    for (index=0; index<NUM_RECEIVERS; index++)
    {
        PRECISION delta_top = 0;
        PRECISION delta_bottom = 0;

        int vindex = index * 2;
        for (int i = 0; i < NUM_COLS; ++i) {
            delta_top += SUQ_array[i] * V_array[vindex * NUM_COLS + i];//[i*num_cols + vindex];
            delta_bottom += SUQ_array[i] * V_array[(vindex + 1) * NUM_COLS + i];//[i*num_cols + vindex+1];
        }

        gains_out[index].x = gains_in[index].x + delta_top;
        gains_out[index].y = gains_in[index].y + delta_bottom;
    }
}










PRECISION2 complex_multiply_CPU(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x * z2.x - z1.y * z2.y, z1.x * z2.y + z1.y * z2.x);
}

PRECISION2 flip_for_i_CPU(const PRECISION2 z)
{
    return MAKE_PRECISION2(-z.y, z.x);
}

PRECISION2 flip_for_neg_i_CPU(const PRECISION2 z)
{
    return MAKE_PRECISION2(z.y, -z.x);
}

PRECISION2 complex_conjugate_CPU(const PRECISION2 z1)
{
    return MAKE_PRECISION2(z1.x, -z1.y);
}

PRECISION2 complex_subtract_CPU(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x - z2.x, z1.y - z2.y);
}




