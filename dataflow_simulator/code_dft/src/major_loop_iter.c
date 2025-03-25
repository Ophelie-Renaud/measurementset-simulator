
#include "major_loop_iter.h"


//naive implementation of 2d conv. Change later to FFT-based for more performance
//NB: This function assumes the same resolution for the image, filter, and result, and that the images are square. The part of the filter convolved is defined by the
//int2 half dimensions, and is centered on the center. Performs 0 padding on the edges
void conv2D(PRECISION* image, PRECISION* filter, PRECISION* result, int IMAGE_DIM, int2 filter_half_dims)
{
	memset(result, 0, sizeof(PRECISION) * IMAGE_DIM * IMAGE_DIM);

    for(int y = 0; y < IMAGE_DIM; ++y)
    {
        for(int x = 0; x < IMAGE_DIM; ++x)
        {
            int idx = y * IMAGE_DIM + x;

            if(image[idx] < 1e-5){
            	continue;
            }

            int center = IMAGE_DIM / 2 + 1;

            for(int filter_y = -filter_half_dims.y; filter_y < filter_half_dims.y; ++filter_y)
            {
                for(int filter_x = -filter_half_dims.x; filter_x < filter_half_dims.x; ++filter_x)
                {
                	int image_x = x + filter_x;
					int image_y = y + filter_y;

					if(image_x < 0 || image_x >= IMAGE_DIM || image_y < 0 || image_y >= IMAGE_DIM){
						continue;
					}

                	int corrected_filter_x = center + (2 * filter_half_dims.x - (filter_x + filter_half_dims.x) - filter_half_dims.x);
                	int corrected_filter_y = center + (2 * filter_half_dims.y - (filter_y + filter_half_dims.y) - filter_half_dims.y);

                	int filter_idx = corrected_filter_y * IMAGE_DIM + corrected_filter_x;
                	int image_idx = image_y * IMAGE_DIM + image_x;

                	result[image_idx] += image[idx] * filter[filter_idx];
                }
            }
        }
    }
}


void save_dirty_image_actor(__attribute__((unused)) int SAVE_DIRTY_IMAGE, int GRID_SIZE, PRECISION *dirty_image, int *cycle_count, Config *config)
{
	printf("Residual image to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, dirty_image);
	save_image_to_file(config, dirty_image, config->dirty_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle_count);
}

void save_estimated_gain_actor(__attribute__((unused)) int SAVE_ESTIMATED_GAIN, int NUM_RECEIVERS, PRECISION2 *estimated_gains, int *cycle_count, Config *config)
{
	
	Complex rotationZ = (Complex){
		.real = estimated_gains[0].x/SQRT(estimated_gains[0].x*estimated_gains[0].x + 
			estimated_gains[0].y * estimated_gains[0].y), 
		.imaginary = - estimated_gains[0].y/SQRT(estimated_gains[0].x*estimated_gains[0].x + 
			estimated_gains[0].y * estimated_gains[0].y)
	};

	printf("UPDATE >>> ROTATED GAINS ...\n\n");

	char buffer[256];
	snprintf(buffer,255,"%scycle_%d_%s",config->output_path,*cycle_count,config->output_gains_file);
	printf("UPDATE >>> Attempting to save image to %s... \n\n", buffer);

	FILE *f = fopen(buffer, "w");

	if(f == NULL)
	{	
		printf(">>> ERROR: Unable to save image to file %s, check file/folder structure exists...\n\n", buffer);
		return;
	}

	for(int i = 0; i < NUM_RECEIVERS; i++)
	{
		Complex rotatedGain = (Complex){
			.real = estimated_gains[i].x * rotationZ.real - estimated_gains[i].y * rotationZ.imaginary,
			.imaginary =  estimated_gains[i].x * rotationZ.imaginary + estimated_gains[i].y * rotationZ.real
		};

		#if SINGLE_PRECISION
            	fprintf(f, "%f %f\n", rotatedGain.real, rotatedGain.imaginary);
			#else
            	fprintf(f, "%lf %lf\n", rotatedGain.real, rotatedGain.imaginary);
		#endif 
	}
	fclose(f);
}

void save_predicted_visibilities_actor(__attribute__((unused)) int SAVE_PREDICTED_VISIBILITIES, int NUM_VISIBILITIES, PRECISION2 *predicted_visibilities, PRECISION3 *vis_uvw_coords, int *cycle_count, Config *config)
{
	char buffer[256];
	snprintf(buffer,255,"%scycle_%d_%s",config->output_path, *cycle_count, config->predicted_vis_output);
    FILE *f = fopen(buffer, "w");
    printf("UPDATE >>> TRYING to save here %s predicted visibilities...\n\n",config->predicted_vis_output);
    if(f == NULL)
    {
        printf("ERROR >>> Unable to save predicted visibilities to file, skipping...\n\n");
        return;
    }

    printf("UPDATE >>> Writing predicted visibilities to file...\n\n");

    PRECISION meters_to_wavelengths = config->frequency_hz / SPEED_OF_LIGHT;
    fprintf(f, "%d\n", NUM_VISIBILITIES);

    // Record individual visibilities
    for(int v = 0; v < NUM_VISIBILITIES; ++v)
    {
        Visibility current_uvw ;//= vis_uvw_coords[v];
        current_uvw.u = vis_uvw_coords[v].x;
        current_uvw.v = vis_uvw_coords[v].y;
        current_uvw.w = vis_uvw_coords[v].z;

        Complex current_vis ;//= predicted_visibilities[v];
        current_vis.real = predicted_visibilities[v].x;
        current_vis.imaginary = predicted_visibilities[v].y;


        current_uvw.u /= meters_to_wavelengths;
        current_uvw.v /= meters_to_wavelengths;
        current_uvw.w /= meters_to_wavelengths;

        if(config->right_ascension)
        {
            current_uvw.u *= -1.0;
            current_uvw.w *= -1.0;
        }

        // u, v, w, real, imag, weight (intensity)
        #if SINGLE_PRECISION
            fprintf(f, "%f %f %f %f %f %f\n", current_uvw.u, current_uvw.v,
                current_uvw.w, current_vis.real, current_vis.imaginary, 1.0);
        #else
            fprintf(f, "%lf %lf %lf %lf %lf %lf\n", current_uvw.u, current_uvw.v,
                current_uvw.w, current_vis.real, current_vis.imaginary, 1.0);
        #endif
    }

    fclose(f);
    printf("UPDATE >>> Predicted visibilities have been successfully written to file...\n\n");
}

void save_residual_image_actor(__attribute__((unused)) int SAVE_RESIDUAL_IMAGE, int GRID_SIZE, PRECISION *residual_image, int *cycle_count, Config *config)
{
	printf("Residual image to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, residual_image);
	save_image_to_file(config, residual_image, config->residual_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle_count);
}

void residual_image_sink(__attribute__((unused)) int SAVE_RESIDUAL_IMAGE, __attribute__((unused)) int GRID_SIZE, __attribute__((unused)) PRECISION *residual_image)
{
	return;
}

void token_sink(__attribute__((unused)) int *token_in)
{
	return;
}

//Saves output image to file  - used for testing
void save_image_to_file(Config *config, PRECISION *image, const char *file_name, int start_x, int start_y, int range_x, int range_y, int cycle)
{
    char buffer[256];
    snprintf(buffer,255,"%scycle_%d_%s",config->output_path,cycle,file_name);
    printf("UPDATE >>> Attempting to save image to %s... \n\n", buffer);

    MD5_Update(sizeof(PRECISION) * range_x * range_y, image);

    FILE *f = fopen(buffer, "w");

    if(f == NULL)
    {   
        printf(">>> ERROR: Unable to save image to file %s, check file/folder structure exists...\n\n", buffer);
        return;
    }

    for(int row = start_y; row < start_y + range_y; ++row)
    {
        for(int col = start_x; col < start_x + range_x; ++col)
        {
            PRECISION pixel = image[row * range_x + col];

            #if SINGLE_PRECISION
                fprintf(f, "%f, ", pixel);
            #else
                fprintf(f, "%lf, ", pixel);
            #endif  
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void save_extracted_sources_actor(int GRID_SIZE, __attribute__((unused)) int SAVE_EXTRACTED_SOURCES, __attribute__((unused)) int NUM_MAX_SOURCES,
		int *num_sources, PRECISION3 *sources, int *cycle_count, Config *config, PRECISION* image_in, PRECISION* clean_psf,
		int2* clean_psf_halfdims, PRECISION* image_out)
{
	memcpy(image_out, image_in, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	PRECISION cell_size = config->cell_size;

	for(int i = 0; i < *num_sources; ++i){
		int row = ROUND((sources[i].y / cell_size) + (GRID_SIZE / 2));
		int col = ROUND((sources[i].x / cell_size) + (GRID_SIZE / 2));

		int idx = col + row * GRID_SIZE;
		image_out[idx] += sources[i].z;
	}

	PRECISION* convolved_output = (PRECISION*)malloc(GRID_SIZE*GRID_SIZE*sizeof(PRECISION));

	conv2D(image_out, clean_psf, convolved_output, GRID_SIZE, *clean_psf_halfdims);

	save_image_to_file(config, convolved_output, config->model_sources_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle_count);

	free(convolved_output);
}

void save_extracted_sources(PRECISION3 *sources, int number_of_sources, const char *path, const char *output_file, int cycle)
{
	char buffer[256];
	snprintf(buffer,250,"%scycle_%d_%s",path,cycle,output_file);
	printf("UPDATE >>> Attempting to save sources to %s... \n\n", buffer);

	FILE *file = fopen(buffer, "w");

	if(file == NULL)
	{
		printf(">>> ERROR: Unable to save sources to file, moving on...\n\n");
		return;
	}

	fprintf(file, "%d\n", number_of_sources);
	for(int index = 0; index < number_of_sources; ++index)
	{
#if SINGLE_PRECISION
		fprintf(file, "%f, %f, %f\n", sources[index].x, sources[index].y, sources[index].z);
#else
		fprintf(file, "%.15f, %.15f, %.15f\n", sources[index].x, sources[index].y, sources[index].z);
#endif
	}

	fclose(file);
}

void save_output(int GRID_SIZE, IN PRECISION* residual, IN PRECISION* model, IN PRECISION* clean_psf, int2* clean_psf_halfdims, IN PRECISION* psf, IN Config *config, IN int* cycle){
	printf("Residual image to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, residual);
	save_image_to_file(config, residual, config->residual_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	printf("Dirty PSF to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, psf);
	save_image_to_file(config, psf, config->psf_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	printf("Clean PSF to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, clean_psf);
	save_image_to_file(config, clean_psf, config->clean_psf_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	printf("Model to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, model);
	save_image_to_file(config, model, config->model_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	PRECISION* final_image = (PRECISION*)malloc(GRID_SIZE * GRID_SIZE * sizeof(PRECISION));
	conv2D(model, clean_psf, final_image, GRID_SIZE, *clean_psf_halfdims);

	printf("Final image to file MD5:");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE, final_image);
	save_image_to_file(config, final_image, config->final_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	free(final_image);
}

void source_list_sink(__attribute__((unused))int MAX_SOURCES, __attribute__((unused)) IN PRECISION3* source_list, __attribute__((unused)) IN int* num_sources){
	return;
}

void delta_visibility_sink(__attribute__((unused))int NUM_VISIBILITIES, __attribute__((unused))IN PRECISION2* visibilities){
	return;
}

void psf_sink(__attribute__((unused))int GRID_SIZE, __attribute__((unused))IN PRECISION *psf)
{
	return;
}

void pass_through_image(int GRID_SIZE, IN PRECISION* image, OUT PRECISION* output_image){
	memcpy(output_image, image, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
}

void save_partial_psf(int GRID_SIZE, PRECISION* psf, int2* extents, int* cycle, Config* config){
	PRECISION* partial_psf = (PRECISION*)malloc(GRID_SIZE * GRID_SIZE * sizeof(PRECISION));
	memset(partial_psf, 0, GRID_SIZE * GRID_SIZE * sizeof(PRECISION));

	int center = GRID_SIZE / 2;

	for(int y = center - extents->y; y < center + extents->y; ++y){
		for(int x = center - extents->x; x < center + extents->x; ++x){
			int idx = y * GRID_SIZE + x;
			partial_psf[idx] = psf[idx];
		}
	}

	save_image_to_file(config, psf, config->clean_psf_image_output, 0, 0, GRID_SIZE, GRID_SIZE, *cycle);

	free(partial_psf);
}
