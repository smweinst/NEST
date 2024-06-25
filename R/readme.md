Get started to use NEST in R
================

#### Step 1: Install the package

``` r
if (!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("smweinst/NEST")
library(NEST)
```

#### Step 2: Prepare the data

Import necessary libraries for NEST:

parallel for efficient parallel computing across multiple processor
cores or machines to speed up computations.

``` r
library(parallel)
library(NEST)
```

This code loads and preprocesses the brain imaging data and phenotype
(demographic) information:

- Brain imaging data (cifti_data) is loaded from the function
  get_CT_data().
- Network labels (network_labels) are loaded from yeo_7 network and
  filtered to remove unspecified indices (-1).
- A list (net) of a vector (net_7) is created for the specific default
  network, indicating vertex membership.

``` r
cifti_data <- load_dataset()
network_labels <- load_yeo_7()
idx <- which(network_labels != -1)
network_label <- network_label[idx]
labeled_data <- cifti_data[, idx]
net_7 <- as.integer(network_label == 7)

net <- list(
    Default <- net_7
)
```

- Phenotype data (age, sex) are loaded and combined into a single list
  (pheno).

``` r
age <- load_pheno("age")
sex <- load_pheno("sex")

pheno <- list(age = age, sex = sex)
```

After loading the data, we can print the shapes of the processed data to
verify their dimensions, ensuring they are correct for subsequent
analysis. In this example, the sample size is 440 and the number of
vertices is 58606.

``` r
print(dim(labeled_data))  # should print 440 58606
print(dim(pheno$age))  # should print 440
print(dim(pheno$sex))  # should print 440
print(dim(net$net_7))  # should print 58606
```

#### Step 3: Define the following dictionary of arguments passed to the nest method. The args can be defined as follows, assuming vertex-wise linear models will be fit to estimate local brain-phenotype associations (i.e., specifying statFun=“lm” in step 4.).

- X: N x P matrix of P imaging features (e.g., vertices) for N
  participants (e.g., N=440, P=58606).
- y: N-dimensional vector of the phenotype of interest (i.e., testing
  enrichment of X-y associations).
- Z: Optional. Specify one or more covariates (matrix with N rows and q
  columns for q covariates). Default is nNULL (no covariates to be
  included).
- FL: Optional. Default is FALSE, set to TRUE to use the Freedman-Lane
  procedure to account for dependence between covariates in permutation.
- n_perm: Optional. Default is 999, with the smallest possible p-value
  of 1/1000.

``` r
args <- list(
    X = labeled_data,
    y = pheno$age,
    Z = pheno$sex,
    type = "coef",
    n.perm = 1
)
```

#### Step 4: Apply NEST to test enrichment of brain-phenotype associations in specified networks.

We can set a variable what_to_return to specify what should be returned
from the NEST function. (Default is to return p-values only) E.g.,
specify “everything” to make nest return p-value, enrichment score and
statistics for null distribution.

``` r
what_to_return <- c("everything")

out <- nest(statFun = "lm", # Use linear regression to get vertex-level test statistics.
    args = args, # Arguments specified above (specific to statFun="lm").
    net.maps = net, # List of binary indicating locations inside (1) or outside (1) network(s) of interest.
    one.sided = TRUE, # Determines whether the enrichment score calculation should consider only the positive alignment (True) or both directions (False).
    n.cores = 1, # Specify the number of CPU cores to be employed for parallel processing tasks within the function
    seed = NULL, # Random seed for reproducible permutation. Default is None. 
    what.to.return = what_to_return # Specify what should be returned.
)
```
