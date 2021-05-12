import numpy as np


def struct_score(DSi, SSi):
    '''Calculate structure score based on Li et al 2012
    method.

    Args:
        DSi (float): Normalized dsRNA-seq (RNaseI resistant) coverage.
        SSi (float): Arcsinh transformed ssRNA-seq (RNaseVI resistant) coverage.
    '''
    return np.log2(DSi + np.sqrt(1 + np.pow(DSi, 2))) - np.log2(SSi + np.sqrt(1 + np.pow(SSi, 2)))


def struct_score_overlap(line):
    # column 3 and 7 should be scores
    # we are going to assume first is DSi and second is SSi
    return struct_score(line[3], line[7])


def struct_score_single(line, DSi=0, SSi=0):
    # calculate struct score for a region that has no overlapping
    # complement RNA-seq
    if DSi:
        DSi = DSi[3]
    elif SSi:
        SSi = SSi[3]
    else:
        raise Exception('Either DSi or SSi must be a list!')
    return struct_score(DSi=DSi, SSi=SSi)




