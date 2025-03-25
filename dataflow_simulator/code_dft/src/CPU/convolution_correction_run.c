
#include "convolution_correction_run.h"

void execute_convolution_correction_actor(int GRID_SIZE, PRECISION* dirty_image_in, PRECISION* prolate, Config *config, PRECISION* dirty_image_out)
{
    int grid_square = GRID_SIZE * GRID_SIZE;
    int row_index, col_index;

    printf("prolate\n");
    MD5_Update(sizeof(PRECISION) * GRID_SIZE / 2, prolate);

    printf("dirty_image_in\n");
    MD5_Update(sizeof(PRECISION) * grid_square, dirty_image_in);

    printf("grid_size: %d\n", GRID_SIZE);

    for(row_index=0; row_index < GRID_SIZE; row_index++)
	{
		for (col_index = 0; col_index < GRID_SIZE; col_index++)
		{
			const int image_index = row_index * GRID_SIZE + col_index;
			const int half_image_size = GRID_SIZE / 2;
			PRECISION taper = prolate[abs(col_index - half_image_size)] * prolate[abs(row_index - half_image_size)];

			if (abs(col_index - half_image_size) == 1229 || abs(row_index - half_image_size) == 1229)
				taper = 0.0;

			dirty_image_out[image_index] = (ABS(taper) > (1E-10)) ? dirty_image_in[image_index] / taper : 0.0;

		}
	}

    printf("convolution_correction MD5\t: ");
    MD5_Update(sizeof(PRECISION) * grid_square, dirty_image_out);
}

void execute_inv_convolution_correction_actor(int GRID_SIZE, PRECISION* dirty_image_in, PRECISION* prolate, Config *config, PRECISION* dirty_image_out)
{
	int grid_square = GRID_SIZE * GRID_SIZE;
	int row_index, col_index;

	printf("prolate\n");
	MD5_Update(sizeof(PRECISION) * GRID_SIZE / 2, prolate);

	printf("dirty_image_in\n");
	MD5_Update(sizeof(PRECISION) * grid_square, dirty_image_in);

	printf("grid_size: %d\n", GRID_SIZE);

	for(row_index=0; row_index < GRID_SIZE; row_index++)
	{
		for (col_index = 0; col_index < GRID_SIZE; col_index++)
		{
			const int image_index = row_index * GRID_SIZE + col_index;
			const int half_image_size = GRID_SIZE / 2;
			PRECISION taper = prolate[abs(col_index - half_image_size)] * prolate[abs(row_index - half_image_size)];

			if (abs(col_index - half_image_size) == 1229 || abs(row_index - half_image_size) == 1229)
				taper = 0.0;

			dirty_image_out[image_index] = dirty_image_in[image_index];//(ABS(taper) > (1E-10)) ? dirty_image_in[image_index] / taper : dirty_image_in[image_index];
		}
	}

	printf("convolution_correction MD5\t: ");
	MD5_Update(sizeof(PRECISION) * grid_square, dirty_image_out);
}
