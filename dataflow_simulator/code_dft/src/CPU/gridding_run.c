
#include "gridding_run.h"


void gridding_actor(int GRID_SIZE, int NUM_VISIBILITIES, int NUM_KERNELS, int TOTAL_KERNEL_SAMPLES,
        PRECISION2* kernels, int2* kernel_supports, PRECISION3* vis_uvw_coords, PRECISION2* visibilities, Config *config, PRECISION2* uv_grid)
{
    memset(uv_grid, 0, sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE);

    printf("UPDATE >>> Gridding on CPU, for %d visibilities...\n\n", NUM_VISIBILITIES);

   // printf("visibilities MD5\t: ");
    //MD5_Update(sizeof(PRECISION2) * NUM_VISIBILITIES, visibilities);

    // kernels MD5
    //printf("kernels MD5\t\t\t: ");
    //MD5_Update(sizeof(PRECISION2) * TOTAL_KERNEL_SAMPLES, kernels);

    // vis_uvw_coords MD5
    //printf("vis_uvw_coords MD5\t: ");
   // MD5_Update(sizeof(PRECISION3) * NUM_VISIBILITIES, vis_uvw_coords);

    // kernel_supports MD5
   // printf("kernel_supports MD5\t: ");
   // MD5_Update(sizeof(int2) * NUM_KERNELS, kernel_supports);

    // num_visibilities
    printf("num_visibilities\t: %d\n", NUM_VISIBILITIES);

    // oversampling
    printf("oversampling\t\t: %d\n", config->oversampling);

    // grid_size
    printf("grid_size\t\t\t: %d\n", GRID_SIZE);

    // uv_scale
    printf("uv_scale\t\t\t: %lf\n", config->uv_scale);

    // w_scale
    printf("w_scale\t\t\t\t: %lf\n", config->w_scale);

    gridding_CPU(uv_grid, (PRECISION2*) kernels, kernel_supports,
        (PRECISION3*) vis_uvw_coords, (PRECISION2*) visibilities, NUM_VISIBILITIES, 16,
        GRID_SIZE, config->uv_scale, config->w_scale);

    printf("UPDATE >>> Gridding complete...\n\n");


    printf("d_uv_grid MD5 \t\t: ");
    MD5_Update(sizeof(PRECISION2) * GRID_SIZE * GRID_SIZE, uv_grid);

}


void gridding_CPU(PRECISION2 *grid, const PRECISION2 *kernel, const int2 *supports,
                  const PRECISION3 *vis_uvw, const PRECISION2 *vis, const int num_vis, const int oversampling,
                  const int grid_size, const double uv_scale, const double w_scale)
{
    // const unsigned int vis_index = blockIdx.x * blockDim.x + threadIdx.x;

    // if(vis_index >= num_vis)
    // 	return;
    for (unsigned int vis_index = 0; vis_index < num_vis; vis_index++)
    {

        // Represents index of w-projection kernel in supports array
        const int plane_index = (int) ROUND(SQRT(ABS(vis_uvw[vis_index].z * w_scale)));

        // Scale visibility uvw into grid coordinate space
        const PRECISION2 grid_coord = MAKE_PRECISION2(
                vis_uvw[vis_index].x * uv_scale,
                vis_uvw[vis_index].y * uv_scale
        );
        const int half_grid_size = grid_size / 2;
        const int half_support = supports[plane_index].x;

        PRECISION conjugate = (vis_uvw[vis_index].z < 0.0) ? -1.0 : 1.0;

        const PRECISION2 snapped_grid_coord = MAKE_PRECISION2(
                ROUND(grid_coord.x * oversampling) / oversampling,
                ROUND(grid_coord.y * oversampling) / oversampling
        );

        const PRECISION2 min_grid_point = MAKE_PRECISION2(
                CEIL(snapped_grid_coord.x - half_support),
                CEIL(snapped_grid_coord.y - half_support)
        );

        const PRECISION2 max_grid_point = MAKE_PRECISION2(
                FLOOR(snapped_grid_coord.x + half_support),
                FLOOR(snapped_grid_coord.y + half_support)
        );
        // PRECISION2 grid_point = MAKE_PRECISION2(0.0, 0.0);
        PRECISION2 convolved = MAKE_PRECISION2(0.0, 0.0);
        PRECISION2 kernel_sample = MAKE_PRECISION2(0.0, 0.0);
        int2 kernel_uv_index = make_int2(0, 0);

        int grid_index = 0;
        int kernel_index = 0;
        int w_kernel_offset = supports[plane_index].y;

        // printf("%lf \t %lf\n", max_grid_point.x - min_grid_point.x, max_grid_point.y - min_grid_point.y);

        for(int grid_v = min_grid_point.y; grid_v <= max_grid_point.y; ++grid_v)
        {
        	if(grid_v < -half_grid_size || grid_v >= half_grid_size){
				continue;
			}

            kernel_uv_index.y = abs((int)ROUND((grid_v - snapped_grid_coord.y) * oversampling));

            for(int grid_u = min_grid_point.x; grid_u <= max_grid_point.x; ++grid_u)
            {
            	if(grid_u < -half_grid_size || grid_u >= half_grid_size){
					continue;
				}

                kernel_uv_index.x = abs((int)ROUND((grid_u - snapped_grid_coord.x) * oversampling));

                kernel_index = w_kernel_offset + kernel_uv_index.y * (half_support + 1)
                                                 * oversampling + kernel_uv_index.x;
                kernel_sample = MAKE_PRECISION2(kernel[kernel_index].x, kernel[kernel_index].y  * conjugate);

                grid_index = (grid_v + half_grid_size) * grid_size + (grid_u + half_grid_size);

                convolved = complex_mult_CPU(vis[vis_index], kernel_sample);

                grid[grid_index].x += convolved.x;
                grid[grid_index].y += convolved.y;
                // atomicAdd(&(grid[grid_index].x), convolved.x);
                // atomicAdd(&(grid[grid_index].y), convolved.y);
            }
        }
    }

    printf("FINISHED GRIDDING\n");
}

void add_visibilities(int NUM_VISIBILITIES, IN PRECISION2* v1, IN PRECISION2* v2, OUT PRECISION2* output){
	printf("Adding visibilities...\n\n");

	for(int i = 0; i < NUM_VISIBILITIES; ++i){
		output[i].x = v1[i].x + v2[i].x;
		output[i].y = v1[i].y + v2[i].y;
	}
}

PRECISION2 complex_mult_CPU(const PRECISION2 z1, const PRECISION2 z2)
{
    return MAKE_PRECISION2(z1.x * z2.x - z1.y * z2.y, z1.y * z2.x + z1.x * z2.y);
}
