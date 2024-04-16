# NEST
 The Python implementation of NEST 

# Illustration of NEST

 <p align="center">
  <img width="100%" height="auto" src="https://github.com/smweinst/NEST/blob/main/Python/NEST_illustration.jpeg">
</p>



## Getting Started:
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
### Prerequisites
```
Python
sci-learn
nibabel
nilearn.plotting as plotting
matplotlib
hcp_utils
numpy
pandas
```
Install requirement dependents
```
pip3 install scipy sklearn pandas matplotlib nibabel nilearn hcp_utils
```
Then install the NEST via pip

```
pip install nest-sw
```

Run the NEST method. An example like

```
from nest-sw import nest

args = {
    'X': data,
    'y': phenotype,
    'Z': covariate,
    'type': 'coef',
    'FL': False,
    'getNull': True,
    'n_perm': 5
}
pval,ES_obs,ES_null,running_sum = NEST(statFun='lm',args=args,net_maps=net_6,n_cores=1,seed=123)
```


## Contact
You can reach out to us regarding your questions , suggestions, and possible collaboration:

Prof. Sarah M. Weinstein: sarah.m.weinstein@temple.edu

## Citation (need update)
If you use our model in any project or publication, please cite our paper [Network Enrichment Significance Testing in Brain-Phenotype Association Studies](https://www.biorxiv.org/content/10.1101/2023.11.10.566593v1.abstract).

```
@article{weinstein2023network,
  title={Network Enrichment Significance Testing in Brain-Phenotype Association Studies},
  author={Weinstein, Sarah M and Vandekar, Simon N and Alexander-Bloch, Aaron F and Raznahan, Armin and Li, Mingyao and Gur, Raquel E and Gur, Ruben C and Roalf, David R and Park, Min Tae M and Chakravarty, Mallar and others},
  journal={bioRxiv},
  pages={2023--11},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```
