
#include "deconvolution_run.h"

unsigned int d_source_counter = 0;

void deconvolution_run(int GRID_SIZE, int CALIBRATION, int NUMBER_MINOR_CYCLES_CAL, int NUMBER_MINOR_CYCLES_IMG, int NUM_MAX_SOURCES, int NUM_MAJOR_CYCLES,
			PRECISION* dirty_image_in, PRECISION* psf_in, PRECISION3* sources_in, int* num_sources_in, Config* config, int* num_sources_out, PRECISION3* sources_out, PRECISION* residual_image)
{
	bool exit_early = false;
    //int col_index,row_index,index,idx,idy;
	int num_sources = 0;//*num_sources_in;
	bool exit_early_temp = false;

	PRECISION3* sources = (PRECISION3*) malloc(sizeof(PRECISION3) * NUM_MAX_SOURCES);
	memset(sources, 0, sizeof(PRECISION3) * NUM_MAX_SOURCES);
	if(num_sources > 0) // occurs only if has sources from previous major cycle
	{
		memcpy(sources, sources_in, sizeof(PRECISION3) * num_sources);
	}

    PRECISION3 max_locals[GRID_SIZE];
	memset(max_locals, 0, sizeof(PRECISION3) * GRID_SIZE);

	//memset(sources_out, 0, sizeof(PRECISION3) * NUM_MAX_SOURCES);

	// scale_dirty_image_by_psf
    for (int row_index =0; row_index< GRID_SIZE; row_index++)
    {
        for (int col_index =0; col_index< GRID_SIZE; col_index++)
        {
			residual_image[row_index * GRID_SIZE + col_index] = dirty_image_in[row_index * GRID_SIZE + col_index] / config->psf_max_value;
        }
    }
	printf("scale_dirty_image_by_psf MD5\t: ");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, residual_image);

	unsigned int cycle_number = 0;

	if(CALIBRATION)
		num_sources = 0;

	//convert existing sources to grid coords
	if(num_sources > 0)
	{
		printf("UPDATE >>> Performing grid conversion on previously found Source coordinates...\n\n");

		// image_to_grid_coords_conversion
        PRECISION cell_size = config->cell_size;
        for(int index=0; index < num_sources; index++)
        {
			sources[index].x = ROUND((sources[index].x / cell_size) + (GRID_SIZE / 2));
			sources[index].y = ROUND((sources[index].y / cell_size) + (GRID_SIZE / 2));
        }
		printf("image_to_grid_coords_conversion MD5\t: ");
		MD5_Update(sizeof(PRECISION3) * num_sources, sources);
	}

	int number_minor_cycles = (CALIBRATION) ? NUMBER_MINOR_CYCLES_CAL : NUMBER_MINOR_CYCLES_IMG;
	printf("UPDATE >>> Performing deconvolution, up to %d minor cycles...\n\n",number_minor_cycles);

	double weak_source_percent = (CALIBRATION) ? config->weak_source_percent_gc : config->weak_source_percent_img;
	while((cycle_number < number_minor_cycles) && (exit_early == false))
	{
		if(cycle_number % 10 == 0)
			printf("UPDATE >>> Performing minor cycle number: %u...\n\n", cycle_number);

		// Find local row maximum via reduction
        // find_max_source_row_reduction<<<reduction_blocks, reduction_threads>>> (d_image, max_locals, GRID_SIZE);
        for (int row_index=0; row_index < GRID_SIZE; row_index++)
        {
            PRECISION3 max = MAKE_PRECISION3(0.0, ABS(residual_image[row_index * GRID_SIZE]), residual_image[row_index * GRID_SIZE]);

            for(int col_index = 1; col_index < GRID_SIZE; col_index++)
            {
                PRECISION current = residual_image[row_index * GRID_SIZE + col_index];
                max.y += ABS(current);
                if(ABS(current) > ABS(max.z))
                {
                    // update m and intensity
                    max.x = (PRECISION) col_index;
                    max.z = current;
                }
            }
            max_locals[row_index] = max;
        }

        // Find final image maximum via column reduction (local maximums array)
        // find_max_source_col_reduction<<<1, 1>>> (d_sources, max_locals, cycle_number, GRID_SIZE, config->loop_gain, weak_source_percent, config->noise_factor);
        for(int col_index = 0; col_index < GRID_SIZE; col_index++)
        {
			if(col_index >= 1) // only single threaded
				continue;

            const double loop_gain=config->loop_gain;
            const double noise_factor=config->noise_factor;
            //obtain max from row and col and clear the y (row) coordinate.
            PRECISION3 max = max_locals[0];
            PRECISION running_avg = max_locals[0].y;
            max.y = 0.0;

            PRECISION3 current;

            for(int index = 1; index < GRID_SIZE; ++index)
            {
                current = max_locals[index];
                running_avg += current.y;
                current.y = index;

                if(ABS(current.z) > ABS(max.z))
                    max = current;
            }

            running_avg /= (GRID_SIZE * GRID_SIZE);
            max.z *= loop_gain;

            // determine whether we drop out and ignore this source

            bool extracting_noise = ABS(max.z) < noise_factor * running_avg * loop_gain;
            bool weak_source = ABS(max.z) < (ABS(sources[0].z) * weak_source_percent);
            exit_early_temp = extracting_noise || weak_source;


            if(!exit_early_temp)
            {
                // source was reasonable, so we keep it
				sources[num_sources] = max;
                ++num_sources;
            }
            else
                exit_early = true;
        }

        // substract_psf_from_image
        for(int idy=0; idy<GRID_SIZE; idy++)
        {
            for(int idx=0; idx<GRID_SIZE; idx++)
            {
                // Determine image coordinates relative to source location
                int2 image_coord = make_int2(sources[num_sources-1].x - GRID_SIZE/2 + idx, sources[num_sources-1].y - GRID_SIZE/2 + idy);

				// image coordinates fall out of bounds
				if(!(image_coord.x < 0 || image_coord.x >= GRID_SIZE || image_coord.y < 0 || image_coord.y >= GRID_SIZE))
                {
                    // Get required psf sample for subtraction
                    const PRECISION psf_weight = psf_in[idy * GRID_SIZE + idx];
                    // Subtract shifted psf sample from image
                    residual_image[image_coord.y * GRID_SIZE + image_coord.x] -= psf_weight  * sources[num_sources-1].z;
                }
            }
        }

        // compress_source
		/*PRECISION3 last_source = sources[num_sources - 1];
		for(int i = num_sources - 2; i >= 0; --i)
		{
		    if((int)last_source.x == (int)sources[i].x && (int)last_source.y == (int)sources[i].y)
		    {
				sources[i].z += last_source.z;
		        --num_sources;
		        break;
		    }
		}*/
//		printf("[%d] %lf, %lf, %lf\n", num_sources-1, sources_out[num_sources-1].x, sources_out[num_sources-1].y, sources_out[num_sources-1].z);

		if(exit_early)
		{
			printf(">>> UPDATE: Terminating minor cycles as now just cleaning noise, cycle number %u...\n\n", cycle_number);
		}

		cycle_number++;
	}

	// Determine how many compressed sources were found
//	cudaMemcpyFromSymbol(&num_sources, d_source_counter, sizeof(unsigned int), 0, cudaMemcpyDeviceToHost);
//	cudaDeviceSynchronize();
	if(num_sources > 0)
	{
		printf("UPDATE >>> Performing conversion on Source coordinates...\n\n");

		for(int index=0; index < num_sources; index ++)
		{
			sources[index].x = (sources[index].x - GRID_SIZE / 2) * config->cell_size;
			sources[index].y = (sources[index].y - GRID_SIZE / 2) * config->cell_size;
		}
	}

	// update outputs
    *num_sources_out = num_sources;
    memcpy(sources_out, sources, sizeof(PRECISION3) * num_sources);

    free(sources);

	printf("sources_out\t\t: ");
	MD5_Update(sizeof(PRECISION3) * num_sources, sources);
}

void hogbom_clean(int GRID_SIZE, int NUM_MINOR_CYCLES, int MAX_SOURCES, IN PRECISION* residual, IN PRECISION* psf, int2* partial_psf_halfdims,
		IN Config* config, IN PRECISION* current_model, OUT int* num_sources_out, OUT PRECISION3* sources_out, OUT PRECISION* output_model){

	memcpy(output_model, current_model, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	*num_sources_out = 0;
        /*
	PRECISION* temp_residual = (PRECISION*) malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	memcpy(temp_residual, residual, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	//maybe do the add in another node and only output the model/deconv of the current iteration here
	memcpy(output_model, current_model, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	*num_sources_out = 0;

	const double WEAK_SOURCE_PERC = config->weak_source_percent_img;
	const double GAIN = config->loop_gain;
	const double NOISE_FACTOR = config->noise_factor;

	for(int i = 0; i < NUM_MINOR_CYCLES; ++i){
		if(i % 50 == 0){
			printf("UPDATE >>> Performing minor cycle number: %u...\n\n", i);
		}

		//locate largest source
		PRECISION curr_max = 0.f;
		int2 max_pos = {.x = 0, .y = 0};
		PRECISION curr_average = 0.f;

		for(int y = 0; y < GRID_SIZE; ++y){
			for(int x = 0; x < GRID_SIZE; ++x){
				int curr_idx = y * GRID_SIZE + x;
				PRECISION curr_val = ABS(temp_residual[curr_idx]);
				curr_average = curr_idx == 0 ? curr_val : curr_average + (curr_val - curr_average) / curr_idx;

				if(curr_val > curr_max){
					curr_max = fabs(temp_residual[curr_idx]);
					max_pos.x = x; max_pos.y = y;
				}
			}
		}

		int max_idx = max_pos.y * GRID_SIZE + max_pos.x;

		//check to see if we should terminate
		bool extracting_noise = ABS(curr_max) < NOISE_FACTOR * curr_average * GAIN;
	    bool weak_source = *num_sources_out > 0 && ABS(curr_max) < (ABS(sources_out[0].z) * WEAK_SOURCE_PERC);

		if(extracting_noise || weak_source){
			break;
		}

		int2 psf_center = {.x = GRID_SIZE / 2, .y = GRID_SIZE / 2};

		float max_val = temp_residual[max_idx] / config->psf_max_value;

		//printf("%f\n", max_val);

		//add to source list. this is only used for the DFT degridder, also the y and x are swapped due to it expecting the row as the first and the column as the second
		sources_out[*num_sources_out].x = max_pos.x; sources_out[*num_sources_out].y = max_pos.y; sources_out[*num_sources_out].z = max_val * GAIN;
		*num_sources_out = *num_sources_out + 1;

		//add to model, which will be used in more recent degridders
		output_model[max_idx] += max_val * GAIN;

		//subtract partial PSF from current residual, offsets used to account for when the support of the psf falls outside of the image FIXXXX
		//int2 image_start_offset = {.x = MIN(max_pos.x, partial_psf_halfdims->x * 30), .y = MIN(max_pos.y, partial_psf_halfdims->y * 30)};
		//int2 image_end_offset = {.x = MIN(GRID_SIZE - max_pos.x, partial_psf_halfdims->x * 30 + 1), .y = MIN(GRID_SIZE - max_pos.y, partial_psf_halfdims->y * 30 + 1)};


		int2 image_start_offset = {.x = MIN(max_pos.x, GRID_SIZE/2), .y = MIN(max_pos.y, GRID_SIZE/2)};
		int2 image_end_offset = {.x = MIN(GRID_SIZE - max_pos.x, GRID_SIZE/2), .y = MIN(GRID_SIZE - max_pos.y, GRID_SIZE/2)};

		for(int y = max_pos.y - image_start_offset.y; y < max_pos.y + image_end_offset.y; ++y){
			for(int x = max_pos.x - image_start_offset.x; x < max_pos.x + image_end_offset.x; ++x){
				int image_idx = y * GRID_SIZE + x;
				int psf_idx = (y - max_pos.y + psf_center.y) * GRID_SIZE + (GRID_SIZE - (x - max_pos.x + psf_center.x));

				temp_residual[image_idx] -= GAIN * max_val * psf[psf_idx] * config->psf_max_value;
			}
		}

		/*for(int i = *num_sources_out - 2; i >= 0; --i)
		{
		    if(max_pos.x == (int)sources_out[i].x && max_pos.y == (int)sources_out[i].y)
		    {
				sources_out[i].z += sources_out[*num_sources_out - 1].z;
		        --*num_sources_out;
		        break;
		    }
		}*/
	/*}*/


	/*for(int i = 0; i < 10000; ++i){
		int x = rand() % GRID_SIZE;
		int y = rand() % GRID_SIZE;

		int idx = y * GRID_SIZE + x;

		output_model[idx] += config->psf_max_value * GAIN;
	}*/


	/*if(*num_sources_out > 0)
	{
		printf("UPDATE >>> Performing conversion on Source coordinates...\n\n");

		for(int index=0; index < *num_sources_out; index ++)
		{
			sources_out[index].x = (sources_out[index].x - GRID_SIZE / 2) * config->cell_size;
			sources_out[index].y = (sources_out[index].y - GRID_SIZE / 2) * config->cell_size;
		}
	}

    free(temp_residual);

	printf("sources_out\t\t: ");
	MD5_Update(sizeof(PRECISION3) * *num_sources_out, sources_out);

	 */
}
