import pandas as pd

SAMPLES = pd.read_csv(
    'samples/SraRunTable.csv').set_index('Sample_Name', drop=False)

rule all:
    input:


