import pandas as pd

SAMPLES = pd.read_csv(
    'samples/SraRunTable.csv').set_index('Sample Name', drop=False)
STRUCT_SCORE_TREATMENTS = ['RNaseOne', 'Rnase V1']
STRUCT_SCORE_SAMPLES = SAMPLES.loc[SAMPLES['treatment'].isin(STRUCT_SCORE_TREATMENTS)]
ADAPTERS = pd.read_table('samples/adapters.tsv')


include: 'rules/download.smk'
include: 'rules/coverage.smk'
include: 'rules/trim_adapters.smk'
include: 'rules/trim_mapped.smk'
include: 'rules/map.smk'

rule all:
    input:
        'output/bowtie2/map_and_convert_struct_reads.done'



