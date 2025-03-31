# Script usage

You can directly use the script instead of downloading the full project:

1. Generate a custom "true sky" model: `python fits_to_csv.py <grid_size> <num_sources> <output_true_sky_fits_path> <output_true_sky_csv_path>`
2. Generate a custom PSF: `python psf.py <grid_size> <output_psf_csv_path>`  (ToDo: faire des PSF uniform weighting, natural weighting )
3. Generate custom gridding/degridding kernels: `python kernel.py <oversampling_factor> <output_k_csv_path>` (ToDo: faire)
4. Display visibility:  `python plot_vis.py <visibility_csv_path>`
5. Display pipelines generated image: `python plot_im <num_major_cycle> <im_csv_folder>`

