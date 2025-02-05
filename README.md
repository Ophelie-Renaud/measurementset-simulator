<div align="center">
  <table width="100%" border="0">
    <tr>
      <td style="text-align: center; vertical-align: middle; width: 33%;">
        <strong>HACKATHON : 3-5/02/2025</strong>
      </td>
      <td style="text-align: center; vertical-align: middle; width: 33%;">
        <img src="https://avalon.ens-lyon.fr/wp-content/uploads/2024/03/Eclat_ecusson.png" height="150">
      </td>
      <td style="text-align: center; vertical-align: middle; width: 33%;">
        <strong>TEAM simulation</strong>
      </td>
    </tr>
  </table>
</div>



# Generating distributed MeasurementSet

![](https://raw.githubusercontent.com/Ophelie-Renaud/vis-generator/refs/heads/main/proj.png)

This project aim to create visibilities from output image which correspond to the reverse steps of standard radio-astronomy imaging pipeline. All the process are contained in the notebook except the `ska_sdp_datamodels`  and  `ska-sdp-func-python` install them before benefiting from this project.

```bash
git clone https://gitlab.com/ska-telescope/sdp/ska-sdp-datamodels.git
cd ska_sdp_datamodels
pip install --extra-index-url https://artefact.skao.int/repository/pypi-internal/simple ska-telmodel
pip install .
cd ..

git clone https://gitlab.com/ska-telescope/sdp/ska-sdp-func-python.git
cd ska-sdp-func-python
pip install .
cd..

pip install python-casacore

pip install notebook
```
Then run the python notebook
```bash
python notebook
```

## Repository structure

```plaintext
├── SGRA_full_gt.fits        # Input image considered as the true sky  
├── spectral_fits/           # Folder containing the image as a spectral cube  
│   ├── *.fits  
├── instrumented_vis_fits/   # Folder containing visibilities with instrument noise  
│   ├── *.fits  
├── spectral_ms/             # Folder containing the distributed MeasurementSet  
│   ├── *.ms  
├── fits_to_vis.ipynb        # Notebook to generate all files from the input image  
```

## Contact  

For questions or feedback, please contact:  
- [Ophélie Renaud](mailto:ophelie.renaud@ens-paris-saclay.fr)

*This project is part of the ECLAT labcom hackathon.*

[![Demo](https://img.shields.io/badge/Live-Demo-blue)](https://ophelie-renaud.github.io/vis-generator/wast.html)

