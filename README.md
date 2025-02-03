# vis-generator

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



## Contact  

For questions or feedback, please contact:  
- [Ophélie Renaud](mailto:ophelie.renaud@ens-paris-saclay.fr)

*This project is part of the ECLAT labcom hackathon.*

[![Demo](https://img.shields.io/badge/Live-Demo-blue)](https://ophelie-renaud.github.io/vis-generator/wast.html)

