import csv
import argparse
import numpy as np
import sys

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bedgraph', help='Path to bedgraph to scale')
    parser.add_argument('output_bedgraph', help='Path to write output bedgraph')
    parser.add_argument('--min_max_norm', default=False, type=bool, help='Apply min max normalization to scores')
    parser.add_argument('--arcsinh', default=False, type=bool, help='Apply asinh transcormation to scores')

    return parser.parse_args()

def read_bedgraph(filepath):
    lines = []
    with open(filepath) as handle:
        reader = csv.reader(handle, sep='\t')
        for row in reader:
            lines.append(row)
    return lines


def reassign_scores(lines, scores):
    for i, line in enumerate(lines):
        line[-1] = scores[i]
    return lines


def min_max_norm(lines):
    scores = [line[-1] for line in lines]
    norm_scores = list(np.linalg.norm(scores))
    return reassign_scores(lines, norm_scores)


def arcsinh_transform(lines):
    scores = [line[-1] for line in lines]
    norm_scores = list(np.arcsinh(scores))
    return reassign_scores(lines, norm_scores)


def write_bedgraph_from_lines(lines, output_path):
    with open(output_path, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerows(lines)
    return output_path


def main():
    args = get_args()
    if not args.arcsinh and not args.min_max_norm:
        # need to set at least one of these
        print('Need to set either arcsinh or min max norm!')
        sys.exit(1)
    elif args.arcsinh and args.min_max_norm:
        print('Need to set either arcsinh or min max norm not both!')
        sys.exit(1)
    else:
        if args.arcsinh:
            trans_func = arcsinh_transform
        else:
            trans_func = min_max_norm
        lines = read_bedgraph(args.input_bedgraph)
        trans_lines = trans_func(lines)

        assert len(lines) == len(trans_lines)
        assert len(lines[0]) == len(trans_lines[0])

        write_bedgraph_from_lines(trans_lines)


if __name__ == '__main__':
    main()




