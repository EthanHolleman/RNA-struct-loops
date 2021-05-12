import numpy as np


def struct_score(DSi, SSi):
    '''Calculate structure score based on Li et al 2012
    method.

    Args:
        DSi (float): Normalized dsRNA-seq (RNaseI resistant) coverage.
        SSi (float): Arcsinh transformed ssRNA-seq (RNaseVI resistant) coverage.
    '''
    return np.log2(DSi + np.sqrt(1 + np.pow(DSi, 2))) - np.log2(SSi + np.sqrt(1 + np.pow(SSi, 2)))

