# vis-generator

This project aim to create visibilies from output image which correspond to the reverse steps of standard radio-astronomy imaging pipeline. All the pocess are contained in the notebook except the `ska_sdp_datamodels` install it before benefiting from the project.
```bash
git clone https://gitlab.com/ska-telescope/sdp/ska-sdp-datamodels.git
cd ska_sdp_datamodels
pip install --extra-index-url https://artefact.skao.int/repository/pypi-internal/simple ska-telmodel
pip install .
```
Then run the python notebook
```bash
python notebook
```

## Contact  

For questions or feedback, please contact:  
- [Ophélie Renaud](mailto:ophelie.renaud@ens-paris-saclay.fr)

*This project is part of the ECLAT labcom hackathon.*

