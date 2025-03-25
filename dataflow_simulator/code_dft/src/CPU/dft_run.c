#include "dft_run.h"

void dft_actor(int NUM_VISIBILITIES, int NUM_MAX_SOURCES, PRECISION3 *sources, PRECISION3 *vis_uvw_coords, int *num_sources, Config *config, PRECISION2 *visibilities)
{
	printf("UPDATE >>> Executing the Direct Fourier Transform algorithm...\n\n");

	printf("num_sources : \t\t%d\n", *num_sources);

	int vis_index;
	for(vis_index=0; vis_index < NUM_VISIBILITIES; vis_index++)
	{
    	const PRECISION two_PI = PI + PI;
    	const PRECISION3 vis = vis_uvw_coords[vis_index];
    	PRECISION3 src;
    	PRECISION2 theta_complex = MAKE_PRECISION2(0.0, 0.0);
    	PRECISION2 source_sum = MAKE_PRECISION2(0.0, 0.0);

    	// For all sources
    	for(int src_indx = 0; src_indx < *num_sources; ++src_indx)
    	{
    		src = sources[src_indx];
    		//Two formula below. uncomment if needed

    		// square root formula (most accurate method)
    		// 	PRECISION term = SQRT(1.0 - (src.x * src.x) - (src.y * src.y));
    		// 	PRECISION image_correction = term;
    		// 	PRECISION w_correction = term - 1.0;

    		// approximation formula - faster but less accurate
    		PRECISION term = 0.5 * ((src.x * src.x) + (src.y * src.y));
    		PRECISION w_correction = -term;
    		PRECISION image_correction = 1.0 - term;

    		PRECISION src_correction = src.z / image_correction;
    		PRECISION theta = (vis.x * src.x + vis.y * src.y + vis.z * w_correction) * two_PI;
//    		SINCOS(theta, &(theta_complex.y), &(theta_complex.x));
			theta_complex.y = SIN(theta);
			theta_complex.x = COS(theta);

    		source_sum.x += theta_complex.x * src_correction;
    		source_sum.y += -theta_complex.y * src_correction;
    	}

    	visibilities[vis_index] = MAKE_PRECISION2(source_sum.x, source_sum.y);

	}
     MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities);
}

