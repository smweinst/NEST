import numpy as np

def check_args(args, required_args, optional_args):
    '''
    args: a dictionary of arguments provided to the function.
    required_args: a list of strings representing the names of required arguments.
    optional_args: a dictionary where keys are the names of optional arguments and values are their default values.
    '''
    # Create a copy of args to avoid modifying the original dictionary
    args = args.copy()

    # Check for missing required arguments
    missing_args = [arg for arg in required_args if arg not in args]
    if missing_args:
        print("\n1 or more missing or misspecified args:")
        for missing_arg in missing_args:
            print(f"\n{missing_arg}\n")
        return False
    else:
        # Apply default values for optional arguments if they are not in args
        missing_args = [arg for arg in optional_args if arg not in args]
        if missing_args:
            print("Using the following default settings: \n")
        # Additional logic for specific default settings
        if 'Z' in optional_args and 'Z' not in args:
            args['Z'] = np.ones(args['y'].shape)
            print("\nargs['Z'] = 1 ---> defaulting to no covariates, just intercept\n")

        if 'type' in optional_args and 'type' not in args:
            args['type'] = 'coef'
            print("\nargs['type'] = 'coef' ---> for statFun = 'lm', defaulting to lm coefficients as T(v)\n")

        if 'FL' in optional_args and 'FL' not in args:
            args['FL'] = False
            print("\nargs['FL'] = False ---> defaulting not to use Freedman-Lane permutation\n")

        if 'getNull' in optional_args and 'getNull' not in args:
            args['getNull'] = True
            print("\nargs['getNull'] = True ---> defaulting to get a null distribution\n")

        if 'n_perm' in optional_args and 'n_perm' not in args:
            args['n_perm'] = 999
            print("\nargs['n_perm'] = 999 ---> defaulting to 999 permutations\n")

        return args