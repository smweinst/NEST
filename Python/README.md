# NEST
 The Python implementation of NEST 

## Getting Started:
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
### Prerequisites
```
Python
scikit-learn
nibabel
nilearn
matplotlib
numpy
pandas
```
Install requirement dependents
```
pip3 install scipy scikit-learn pandas matplotlib nibabel nilearn
```
Then install the NEST via pip

```
pip install nest-sw
```

Example can be found at [here](https://github.com/smweinst/NEST/blob/main/Python/example.ipynb).

## Details about `statFun` options
`statFun` determines how location-specific test statistics ("T(v)") are computed. We've provided several built-in options (e.g., `statFun='lm'` or `statFun='gam_deltaRsq'`), and there is also an option to specify a custom one (`statFun='custom'` -- see example [here](https://github.com/smweinst/NEST/blob/main/Python/example.ipynb).

**Note**: for `statFun='gam_deltaRsq'`, the result may be different than the corresponding implementation in R due to differences between R and python implementations of GAM, so please use this version with caution (or make adjustments as needed by using a custom function).

## Citation
If you use our method in any project or publication, please use the following citation [Network Enrichment Significance Testing in Brain-Phenotype Association Studies](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.26714).

```
@article{weinstein2024network,
  title={Network enrichment significance testing in brain--phenotype association studies},
  author={Weinstein, Sarah M and Vandekar, Simon N and Li, Bin and Alexander-Bloch, Aaron F and Raznahan, Armin and Li, Mingyao and Gur, Raquel E and Gur, Ruben C and Roalf, David R and Park, Min Tae M and others},
  journal={Human Brain Mapping},
  volume={45},
  number={8},
  pages={e26714},
  year={2024},
  publisher={Wiley Online Library}
}
```
