# Dataflow Simulator

This project was initially conceived sequentially, following dataflow semantics, but was not designed to exibit parallelism, although this is inherent to the model. Here, we're looking for functionality first.

## Usage

- Follow the steps in the notebook to produce visibility via a degridder dataflow and reconstruct the image via the DFT pipeline :arrow_right: `jupyter notebook`.

## Script usage

You can directly use the script instead of downloading the full project:

1. Generate a custom "true sky" model :arrow_right: `python fits_to_csv.py <grid_size> <num_sources> <output_true_sky_fits_path> <output_true_sky_csv_path>`
   - **INPUT**: parameters `grid_size` *(e.g.: 512)*, `num_sources` *(e.g.: 3)*.
   - **OUTPUT**: Display figure $I(l,m,Intensity)$ and generate `your_true_sky_name.csv` and `your_true_sky_name.fits`.

2. Generate a custom PSF :arrow_right: `python psf.py <grid_size> <output_psf_csv_path>`  (ToDo: faire des PSF uniform weighting, natural weighting )
   - **INPUT**: parameters `grid_size` *(e.g.: 512)*.
   - **OUTPUT**: Display figure $I(l,m,Intensity)$ and generate `your_psf_name.csv` .

3. Generate custom gridding/degridding kernels :arrow_right: `python kernel.py <grid_size> <num_kernels> <oversampling_factor> <kernel_support> <beta> <kernel_width> <baseline_max> <frequency_hz> <output_k_csv_path>` .
   - **INPUT**: parameters `grid_size` *(e.g.: 512)*, `num_kernels` *(e.g.: 17)*, `oversampling_factor` *(e.g.: 16)*, `kernel_support` *(e.g.: 6)*, `beta` *(e.g.: 2)*, `kernel_width` *(e.g.: 4)*, `baseline_max` *(e.g.: 1000)*, `frequency_hz` *(e.g.: 1.4E9)*.
   - **OUTPUT**: Display figure $k(u,v,real, imag,Intensity)$ and generate `your_kernel_name.csv`.

4. Display visibility :arrow_right:  `python plot_vis.py <visibility_csv_path>`
   - **INPUT**: path `your_visibilities.csv`.
   - **OUTPUT**: Display figure $V(u,v,intensity)$.

5. Display pipelines generated image :arrow_right: `python plot_im <num_major_cycle> <im_csv_folder>`
   - **INPUT**: path `your_reconstructed_image.csv`.
   - **OUTPUT**: Display figure $I(l,m,Intensity)$.

