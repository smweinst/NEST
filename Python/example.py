# Step 1: import the NEST package after the installation.
from nest-sw import nest

# Step 2: Preparing the data and define the following dictionary of arguments passed to the nest method. 

# The args can be defined as
args = {
    'X': data,  # Hypothetical predictor variables, e.g., the 2-D brain-scan imaging data shaped as (sample_size, feature_size) in a numpy array format.
    'y': phenotype, # Hypothetical response variable, e.g., phenotypic data for samples corresponding to images, shaped as (sample_size, 1).
    'Z': covariate,  # Optionally specifying covariates, set to NULL here if you want to explicitly state no covariates.
    'FL': False,  # Optional, set to True to use Freedman-Lane permutation.
    'getNull': True,  # Optional, defult is True to explicitly stating to get a null distribution.
    'n_perm': 5   # Optional, default is 999. Override to specify the number of permutations. 
}

# Step 3: Loading the labels of networks corresponding to each location in X, shaped as a binary array (feature_size,) for a specific parcellation. e.g., the default network from yeo Yeo2011_7networks

# network_labels = load_some_data()

# Step 4: Calling the NEST method and get a list of return values:

output = nest(statFun='lm',args=args,net_maps=network_labels)

print('The p-value from tests of enrichment is', output[0])

# output[0]: p-value
# output[1]: observed enrichment score
# output[2]: enrichment score for null distribution