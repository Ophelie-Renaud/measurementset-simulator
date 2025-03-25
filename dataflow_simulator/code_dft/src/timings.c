#include "timings.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

void save_timings(int size, char* filename, clock_t* timings){
	FILE* file;
	file = fopen(filename, "w");

	for(int i = 0; i < size - 1; ++i){
		fprintf(file, "%ld, ", timings[i]);
	}
	fprintf(file, "%ld", timings[size - 1]);

	fclose(file);
}

void set_target(int norig, int ntarget, void* orig, void* target, int size){
	if(norig >= ntarget){
		memcpy(target, orig, ntarget * size);
		return;
	}

	int mult = ntarget / norig;
	int dif = ntarget % norig;

	for(int i = 0; i < mult; ++i){
		int offset = i * norig;
		memcpy(target + offset, orig, norig * size);
	}
	int offset = mult * norig;
	memcpy(target + offset, orig, dif * size);
}

void time_constant_setups(int NUM_SAMPLES){
	clock_t* config_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));
	clock_t* config_sequel_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));
	clock_t* dgkernel_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));
	clock_t* gkernel_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));
	clock_t* gains_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	Config config;
	Config config_out;
	double max_psf = 1.;

	int numkernels = 17;
	int numsamples = 108800;

	int num_receivers = 512;
	int num_baselines = num_receivers*(num_receivers-1)/2;

	int2* degridding_kernel_supports = (int2*)malloc(sizeof(int2) * numkernels);
	int2* gridding_kernel_supports = (int2*)malloc(sizeof(int2) * numkernels);
	PRECISION2* degridding_kernels = (PRECISION2*)malloc(sizeof(PRECISION2) * numsamples);
	PRECISION2* gridding_kernels = (PRECISION2*)malloc(sizeof(PRECISION2) * numsamples);
	PRECISION2* gains = (PRECISION2*)malloc(sizeof(PRECISION2) * num_receivers);
	int2* receiver_pairs = (int2*)malloc(sizeof(int2) * num_baselines);

	for(int i = 0; i < NUM_SAMPLES; ++i){

		start = clock();
		config_struct_set_up(2048, 17, &config);
		end = clock();
		config_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;

		start = clock();
		config_struct_set_up_sequel(&config, &max_psf, &config_out);
		end = clock();
		config_sequel_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;

		start = clock();
		degridding_kernel_host_set_up(numkernels, numsamples, &config_out, degridding_kernel_supports, degridding_kernels);
		end = clock();
		dgkernel_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;

		start = clock();
		kernel_host_set_up(numkernels, numsamples, &config_out, gridding_kernel_supports, gridding_kernels);
		end = clock();
		gkernel_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;

		start = clock();
		gains_host_set_up(num_receivers, num_baselines, &config, gains, receiver_pairs);
		end = clock();
		gains_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "config_timings", config_timings);
	save_timings(NUM_SAMPLES, "config_sequel_timings", config_sequel_timings);
	save_timings(NUM_SAMPLES, "dgkernel_timings", dgkernel_timings);
	save_timings(NUM_SAMPLES, "gkernel_timings", gkernel_timings);

	free(degridding_kernel_supports);
	free(gridding_kernel_supports);
	free(degridding_kernels);
	free(gridding_kernels);
	free(config_timings);
	free(config_sequel_timings);
	free(dgkernel_timings);
	free(gkernel_timings);
	free(gains_timings);
	free(gains);
	free(receiver_pairs);
}

void time_gridsize_setups(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* prolate_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION* prolate = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE / 2);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		correction_set_up(GRID_SIZE, prolate);
		end = clock();
		prolate_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "prolate_setup_timings", prolate_timings);
	free(prolate_timings);
	free(prolate);
}

void time_visibility_setups(int NUM_SAMPLES, int NUM_VISIBILITIES){
	clock_t* visibility_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	Config config;
	config_struct_set_up(2048, 17, &config);

	PRECISION2* visibilities;
	PRECISION3* vis_coords;

	for(int j = 0; j < 5; ++j){
		char input_vis_filename[100];

		sprintf(input_vis_filename, "../data/input/GLEAM_small_visibilities_corrupted_%d.csv", j+1);
		config.visibility_source_file = input_vis_filename;
		printf("allocating\n");
		visibilities = (PRECISION2*)malloc(NUM_VISIBILITIES * (j+1) * sizeof(PRECISION2));
		vis_coords = (PRECISION3*)malloc(NUM_VISIBILITIES * (j+1) * sizeof(PRECISION3));

		for(int i = 0; i < NUM_SAMPLES; ++i){
			start = clock();
			visibility_host_set_up(NUM_VISIBILITIES * (j + 1), &config, vis_coords, visibilities);
			end = clock();
			visibility_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
		}

		char output_filename_vis[200];

		sprintf(output_filename_vis, "visibility_setup_timings_%d", j+1);

		save_timings(NUM_SAMPLES, output_filename_vis, visibility_timings);
		printf("Freeing\n");
		free(visibilities);
		free(vis_coords);
	}

	free(visibility_timings);
}

void time_save_output(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* saveoutput_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION* residual = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	PRECISION* model = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	memset(model, 0, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	for(int i = 0; i < 100; ++i){
		int idx = (rand() % GRID_SIZE) * GRID_SIZE + (rand() % GRID_SIZE);
		model[idx] = 1.;
	}

	PRECISION* clean_psf = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* psf = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	int2 clean_psf_halfdims = {.x = 100, .y = 100};

	Config config;
	config_struct_set_up(2048, 17, &config);

	int cycle = 100;

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		save_output(GRID_SIZE, residual, model, clean_psf, &clean_psf_halfdims, psf, &config, &cycle);
		end = clock();
		saveoutput_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "save_output_timings", saveoutput_timings);
	free(saveoutput_timings);
	free(residual);
	free(model);
	free(clean_psf);
	free(psf);
}

void time_dft(int NUM_SAMPLES, int NUM_MINOR_CYCLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES){
	PRECISION2* orig_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_ACTUAL_VISIBILITIES);
	PRECISION2* target_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION3* orig_vis_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_ACTUAL_VISIBILITIES);
	PRECISION3* target_vis_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	clock_t* dft_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	int grid_size = 2458;

	Config config;
	config_struct_set_up(grid_size, 17, &config);

	visibility_host_set_up(NUM_ACTUAL_VISIBILITIES, &config, orig_vis_coords, orig_vis);

	set_target(NUM_ACTUAL_VISIBILITIES, NUM_VISIBILITIES, (void*)orig_vis_coords, (void*)target_vis_coords, sizeof(PRECISION3));
	set_target(NUM_ACTUAL_VISIBILITIES, NUM_VISIBILITIES, (void*)orig_vis, (void*)target_vis, sizeof(PRECISION2));

	PRECISION3* sources = (PRECISION3*)malloc(NUM_MINOR_CYCLES * sizeof(PRECISION3));
	for(int i = 0; i < NUM_MINOR_CYCLES; ++i){
		sources[i].z = 1.f;
		sources[i].x = rand() % grid_size;
		sources[i].y = rand() % grid_size;
	}
	int nsources = NUM_MINOR_CYCLES;

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		dft_actor(NUM_VISIBILITIES, nsources, sources, target_vis_coords, &nsources, &config, target_vis);
		end = clock();
		dft_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "dft_timings", dft_timings);

	free(sources);
	free(orig_vis);
	free(target_vis);
	free(orig_vis_coords);
	free(target_vis_coords);
	free(dft_timings);
}

void time_gains_application(int NUM_SAMPLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES){
	PRECISION2* orig_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_ACTUAL_VISIBILITIES);
	PRECISION2* target_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION3* orig_vis_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_ACTUAL_VISIBILITIES);
	clock_t* gains_apply_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));
	clock_t* recip_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	Config config;
	config_struct_set_up(2458, 17, &config);

	int num_receivers = 512;
	int num_baselines = num_receivers*(num_receivers-1)/2;
	PRECISION2* gains = (PRECISION2*)malloc(sizeof(PRECISION2) * num_receivers);
	PRECISION2* gains_out = (PRECISION2*)malloc(sizeof(PRECISION2) * num_receivers);
	int2* receiver_pairs = (int2*)malloc(sizeof(int2) * num_baselines);

	visibility_host_set_up(NUM_ACTUAL_VISIBILITIES, &config, orig_vis_coords, orig_vis);
	gains_host_set_up(num_receivers, num_baselines, &config, gains, receiver_pairs);

	set_target(NUM_ACTUAL_VISIBILITIES, NUM_VISIBILITIES, (void*)orig_vis, (void*)target_vis, sizeof(PRECISION2));

	PRECISION2* sub_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	memcpy(sub_vis, target_vis, sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION2* out_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		apply_gains_actor(0, num_receivers, num_baselines, NUM_VISIBILITIES, target_vis, sub_vis, gains, receiver_pairs, &config, out_vis);
		end = clock();
		gains_apply_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;

		start = clock();
		reciprocal_transform_actor(num_receivers, gains, &config, gains_out);
		end = clock();
		recip_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "gains_apply_timings", gains_apply_timings);
	save_timings(NUM_SAMPLES, "gains_reciprocal_transform_timings", recip_timings);

	free(gains_out);
	free(orig_vis);
	free(sub_vis);
	free(out_vis);
	free(target_vis);
	free(orig_vis_coords);
	free(gains);
	free(receiver_pairs);
	free(gains_apply_timings);
	free(recip_timings);
}

void time_add_visibilities(int NUM_SAMPLES, int NUM_VISIBILITIES){
	clock_t* addvis_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION2* v1 = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION2* v2 = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION2* out = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		add_visibilities(NUM_VISIBILITIES, v1, v2, out);
		end = clock();
		addvis_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "addvis_timings", addvis_timings);
	free(addvis_timings);
	free(v1);
	free(v2);
	free(out);
}

void time_prolate(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* prolate_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION* dirty_image_in = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* dirty_image_out = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* prolate = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE / 2);

	Config config;
	config_struct_set_up(GRID_SIZE, 17, &config);

	correction_set_up(GRID_SIZE, prolate);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		execute_convolution_correction_actor(GRID_SIZE, dirty_image_in, prolate, &config, dirty_image_out);
		end = clock();
		prolate_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "prolate_timings", prolate_timings);
	free(prolate_timings);
	free(dirty_image_in);
	free(dirty_image_out);
	free(prolate);
}

void time_finegrid(int NUM_SAMPLES, int NUM_VISIBILITIES, int NUM_ACTUAL_VISIBILITIES){
	clock_t* finegrid_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	int grid_size = 2458;
	int oversampling_factor = 16;
	int num_output_visibilities;

	PRECISION2* orig_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_ACTUAL_VISIBILITIES);
	PRECISION2* target_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION3* orig_vis_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_ACTUAL_VISIBILITIES);
	PRECISION3* target_vis_coords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);
	PRECISION2* out_vis = (PRECISION2*)malloc(sizeof(PRECISION2) * NUM_VISIBILITIES);
	PRECISION3* out_viscoords = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_VISIBILITIES);

	Config config;
	config_struct_set_up(2458, 17, &config);

	visibility_host_set_up(NUM_ACTUAL_VISIBILITIES, &config, orig_vis_coords, orig_vis);

	set_target(NUM_ACTUAL_VISIBILITIES, NUM_VISIBILITIES, (void*)orig_vis_coords, (void*)target_vis_coords, sizeof(PRECISION3));
	set_target(NUM_ACTUAL_VISIBILITIES, NUM_VISIBILITIES, (void*)orig_vis, (void*)target_vis, sizeof(PRECISION2));

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		correct_to_finegrid(NUM_VISIBILITIES, grid_size, oversampling_factor, 1, target_vis_coords, target_vis, &config, out_vis, out_viscoords, &num_output_visibilities);
		end = clock();
		finegrid_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "correct_to_finegrid_timings", finegrid_timings);

	free(out_viscoords);
	free(out_vis);
	free(orig_vis);
	free(target_vis);
	free(orig_vis_coords);
	free(target_vis_coords);
	free(finegrid_timings);
}

void time_subtract_ispace(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* subimage_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION* i1 = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* i2 = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* out = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		subtract_image_space(GRID_SIZE, i1, i2, out);
		end = clock();
		subimage_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "subtraction_imagespace_timings", subimage_timings);
	free(subimage_timings);
	free(i1);
	free(i2);
	free(out);
}

void time_fftshift(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* fftshift_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION2* in = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
	PRECISION2* out = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

	Config config;
	config_struct_set_up(2458, 17, &config);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		fft_shift_complex_to_complex_actor(GRID_SIZE, in, &config, out);
		end = clock();
		fftshift_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "fftshift_timings", fftshift_timings);
	free(fftshift_timings);
	free(in);
	free(out);
}

void time_fft(int NUM_SAMPLES, int GRID_SIZE){
	clock_t* fft_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION2* in = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);
	for(int i = 0; i < GRID_SIZE * GRID_SIZE; ++i){
		in[i].x = (float)rand()/(float)RAND_MAX;
		in[i].y = (float)rand()/(float)RAND_MAX;
	}

	PRECISION2* out = (PRECISION2*)malloc(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

	Config config;
	config_struct_set_up(2458, 17, &config);

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		CUFFT_EXECUTE_FORWARD_C2C_actor(GRID_SIZE, in, out);
		end = clock();
		fft_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "fft_timings", fft_timings);
	free(fft_timings);
	free(in);
	free(out);
}

void time_hogbom(int NUM_SAMPLES, int GRID_SIZE, int NUM_MINOR_CYCLES){
	clock_t* clean_timings = (clock_t*)malloc(NUM_SAMPLES * sizeof(clock_t));

	clock_t start, end;
	clock_t CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

	PRECISION* in = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	for(int i = 0; i < GRID_SIZE * GRID_SIZE; ++i){
		in[i] = (float)rand()/(float)RAND_MAX;
	}

	PRECISION* out = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION* psf = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	for(int i = 0; i < GRID_SIZE * GRID_SIZE; ++i){
		psf[i] = (float)rand()/(float)RAND_MAX;
	}
	int2 psf_halfdims = {.x = 50, .y = 50};
	PRECISION* current_model = (PRECISION*)malloc(sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	memset(current_model, 0, sizeof(PRECISION) * GRID_SIZE * GRID_SIZE);
	PRECISION3* sources_out = (PRECISION3*)malloc(sizeof(PRECISION3) * NUM_MINOR_CYCLES);

	int num_sources_out;

	Config config;
	config_struct_set_up(GRID_SIZE, 17, &config);
	config.weak_source_percent_gc = 0;
	config.weak_source_percent_img = 0;
	config.psf_max_value = 1.f;

	for(int i = 0; i < NUM_SAMPLES; ++i){
		start = clock();
		hogbom_clean(GRID_SIZE, NUM_MINOR_CYCLES, NUM_MINOR_CYCLES, in, psf, &psf_halfdims,	&config, current_model, &num_sources_out, sources_out, out);
		end = clock();
		clean_timings[i] = ((double) (end - start)) / CLOCKS_PER_MS + 0.5;
	}

	save_timings(NUM_SAMPLES, "clean_timings", clean_timings);
	free(psf);
	free(current_model);
	free(sources_out);
	free(clean_timings);
	free(in);
	free(out);
}
