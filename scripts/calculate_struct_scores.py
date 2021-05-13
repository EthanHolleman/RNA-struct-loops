import numpy as np
import csv
import argparse

# assumes scores have already beed normalized / scaled
# problem with this kind of thinking is that overlaps
# will not be complete so need to make sure that the
# first file requires complete overlap?

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('overlap', 
                        help='Path to intersescted dsRNA-seq ssRNA-seq bedgraph with the intersection in that respective order.'
                        )
    parser.add_argument('dsRNA_seq',
                        help='Path to bedgraph containing dsRNA_seq regions with no overlap in the ssRNA-seq bedgraph.'
                        )
    parser.add_argument('ssRNA_seq',
                        help='Path to bedgraph containing ssRNA_seq regions with no overlap in the dsRNA_seq bedgraph.'
                        )
    parser.add_argument('output', help='Path to write output bedgraph to')

    return parser.parse_args()


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
    a_start, a_end = line[1], line[2]
    b_start, b_end = line[5], line[6]

    left_gap = a_start - b_start

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


def read_bedgraph(filepath):
    with open(filepath) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            yield row






def main():
    args = get_args()
    with open(args.output, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        # do intersection
        for row in read_bedgraph(args.overlap):
            ss = struct_score_overlap(row)
            writer.writerow(, s)  # region needs to come from one of the files









