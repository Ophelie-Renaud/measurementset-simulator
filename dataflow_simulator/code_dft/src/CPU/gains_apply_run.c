
#include "gains_apply_run.h"

void reciprocal_transform_actor(int NUM_RECEIVERS,
                                PRECISION2* gains_in, Config *config, PRECISION2* gains_out)
{
    int i;
    PRECISION2 z;
    PRECISION real,imag;

    for (i=0; i<NUM_RECEIVERS; i++)
    {
        z=gains_in[i];
        real = z.x / (z.x * z.x + z.y * z.y);
        imag = z.y / (z.x * z.x + z.y * z.y);
        gains_out[i] = MAKE_PRECISION2(real, -imag);
    }
}

void apply_gains_actor(int CALIBRATION, int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES,
                       PRECISION2* measured_vis, PRECISION2* visibilities_in, PRECISION2* gains, int2* receiver_pairs, Config *config, PRECISION2* visibilities_out)
{
    printf("UPDATE >>> Applying gains... \n\n");
    printf("measured_vis\n");
    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, measured_vis);

    printf("gains\n");
    MD5_Update(sizeof(PRECISION2) * NUM_RECEIVERS, gains);

    if(!CALIBRATION)
    {
        //apply gain calibration to update gains between measured and predicted visibilities
        int vis_index;
        for (vis_index=0; vis_index<NUM_VISIBILITIES; vis_index++) {
            int baselineNumber = vis_index % NUM_BASELINES;

            //THIS ASSUMES THAT WE HAVE GAINS AS THE RECIPROCAL
            PRECISION2 gains_a_recip = gains[receiver_pairs[baselineNumber].x];
            PRECISION2 gains_b_recip_conj = complex_conjugate_apply(gains[receiver_pairs[baselineNumber].y]);


            PRECISION2 gains_product_recip = complex_multiply_apply(gains_a_recip, gains_b_recip_conj);
            PRECISION2 measured_with_gains = complex_multiply_apply(measured_vis[vis_index], gains_product_recip);
            visibilities_out[vis_index] = complex_subtract_apply(measured_with_gains, visibilities_in[vis_index]);
        }

    }
    else
    {
        int vis_index;
        for (vis_index=0; vis_index<NUM_VISIBILITIES ; vis_index++)
        {
            int baselineNumber =  vis_index % NUM_BASELINES;

            //THIS ASSUMES THAT WE HAVE GAINS AS THE RECIPROCAL
            PRECISION2 gains_a_recip = gains[receiver_pairs[baselineNumber].x];
            PRECISION2 gains_b_recip_conj = complex_conjugate_apply(gains[receiver_pairs[baselineNumber].y]);


            PRECISION2 gains_product_recip = complex_multiply_apply(gains_a_recip, gains_b_recip_conj);
            PRECISION2 measured_with_gains = complex_multiply_apply(measured_vis[vis_index], gains_product_recip);
            visibilities_out[vis_index] = measured_with_gains;
        }
    }


    printf("Sum d_visibilities\n");
    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities_out);
}

void subtract_from_measurements(int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES,
                    IN PRECISION2* measured_vis, IN PRECISION2* visibilities_in, IN PRECISION2* gains, IN int2* receiver_pairs, IN Config *config, OUT PRECISION2* visibilities_out){
    printf("UPDATE >>> Subtracting visibilities from measurements... \n\n");
    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, measured_vis);

    printf("gains\n");
    MD5_Update(sizeof(PRECISION2) * NUM_RECEIVERS, gains);

    int vis_index;
    for(vis_index=0; vis_index<NUM_VISIBILITIES ; vis_index++)
    {
        int baselineNumber =  vis_index % NUM_BASELINES;

        //THIS ASSUMES THAT WE HAVE GAINS AS THE RECIPROCAL
        PRECISION2 gains_a_recip = gains[receiver_pairs[baselineNumber].x];
        PRECISION2 gains_b_recip_conj = complex_conjugate_apply(gains[receiver_pairs[baselineNumber].y]);


        PRECISION2 gains_product_recip = complex_multiply_apply(gains_a_recip, gains_b_recip_conj);
        PRECISION2 measured_with_gains = complex_multiply_apply(measured_vis[vis_index], gains_product_recip);
        //printf("%f %f - %f %f\n", measured_with_gains.x, measured_with_gains.y, visibilities_in[vis_index].x, visibilities_in[vis_index].y);
        visibilities_out[vis_index] = complex_subtract_apply(measured_with_gains, visibilities_in[vis_index]);
    }

    printf("Sum d_visibilities\n");
    MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities_out);
}

void apply_gains_only(int NUM_RECEIVERS, int NUM_BASELINES, int NUM_VISIBILITIES, IN PRECISION2* measured_vis, IN PRECISION2* gains,
			IN int2* receiver_pairs, IN Config *config, OUT PRECISION2* visibilities_out){
	printf("UPDATE >>> Applying only gains on visibilities from measurements... \n\n");
	MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, measured_vis);
	MD5_Update(sizeof(PRECISION2) * NUM_RECEIVERS, gains);

	int vis_index;
	for(vis_index=0; vis_index<NUM_VISIBILITIES ; vis_index++)
	{
		int baselineNumber =  vis_index % NUM_BASELINES;

		PRECISION2 gains_a_recip = gains[receiver_pairs[baselineNumber].x];
		PRECISION2 gains_b_recip_conj = complex_conjugate_apply(gains[receiver_pairs[baselineNumber].y]);
		PRECISION2 gains_product_recip = complex_multiply_apply(gains_a_recip, gains_b_recip_conj);
		PRECISION2 measured_with_gains = complex_multiply_apply(measured_vis[vis_index], gains_product_recip);

		//ignore gains for now since we don't have the files and are just multiplying by 1 anyways
		visibilities_out[vis_index] = measured_vis[vis_index];//measured_with_gains;
	}

	MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities_out);
}

void subtract_from_measurements_finegrid(int NUM_VISIBILITIES, IN PRECISION2* gridded_measured_vis, IN PRECISION2* input_vis,
		int* num_finegrid_vis, OUT PRECISION2* delta_vis){
	printf("UPDATE >>> Subtracting finegrid visibilities from measurements... \n\n");
	MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, gridded_measured_vis);

	bool keepold = false;

	/*for(int i = 0; i < *num_finegrid_vis; ++i){
		if(input_vis[i].x > 0 || input_vis[i].y > 0){
			keepold = true;
		}
	}*/

	for(int i = 0; i < *num_finegrid_vis; ++i){
		if(keepold){
			delta_vis[i] = input_vis[i];
		}
		else{
			delta_vis[i] = complex_subtract_apply(gridded_measured_vis[i], input_vis[i]);
		}

	}

	MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, delta_vis);
}

PRECISION2 complex_multiply_apply(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x * z2.x - z1.y * z2.y, z1.x * z2.y + z1.y * z2.x);
}

PRECISION2 complex_conjugate_apply(const PRECISION2 z1)
{
    return MAKE_PRECISION2(z1.x, -z1.y);
}

PRECISION2 complex_subtract_apply(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x - z2.x, z1.y - z2.y);
}

