import pandas as pd

SAMPLES = pd.read_csv(
    'samples/SraRunTable.csv').set_index('Run', drop=False)
STRUCT_SCORE_TREATMENTS = ['RNaseOne', 'Rnase V1']
STRUCT_SCORE_SAMPLES = SAMPLES.loc[SAMPLES['treatment'].isin(STRUCT_SCORE_TREATMENTS)]
dsRNA_SAMPLES = SAMPLES.loc[SAMPLES['treatment'] == STRUCT_SCORE_TREATMENTS[0]]
ssRNA_SAMPLES = SAMPLES.loc[SAMPLES['treatment'] == STRUCT_SCORE_TREATMENTS[1]]


ADAPTERS = pd.read_table('samples/adapters.tsv')


wildcard_constraints:
   sample = '(' + '|'.join(SAMPLES['Run'].tolist()) + ')',
   sample_dsRNA = '(' + '|'.join(SAMPLES['Run'].tolist()) + ')',
   sample_ssRNA = '(' + '|'.join(SAMPLES['Run'].tolist()) + ')',


include: 'rules/download.smk'
include: 'rules/coverage.smk'
include: 'rules/trim_adapters.smk'
include: 'rules/trim_mapped.smk'
include: 'rules/map.smk'
include: 'rules/struct_score.smk'

rule all:
    input:
        #'output/bowtie2/map_and_convert_struct_reads.done'
        'output/struct_score/final_score_files.done'



