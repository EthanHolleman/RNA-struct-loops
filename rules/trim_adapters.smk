
ALL_ADAPTER_PARAMS = ' '.join([f'-a {seq}' for seq in ADAPTERS['Sequence'].tolist()])

rule trim_all_adapters_single:
    conda:
        '../envs/cutadapt.yml'
    input:
        'output/sample/fastq/single/{sample}.fastq'
    output:
        'output/trimmed/single/{sample}.trim.fastq'
    params:
        all_adapters=lambda wildcards: ALL_ADAPTER_PARAMS
    shell:'''
    cutadapt {params.all_adapters} -o {output} {input}
    '''


rule trim_all_adapters_paired:
    conda:
        '../envs/cutadapt.yml'
    input:
        'output/sample/fastq/paired/{sample}_{mate}.fastq'
    output:
        'output/trimmed/paired/{sample}_{mate}.trim.fastq'
    params:
        all_adapters=lambda wildcards: ALL_ADAPTER_PARAMS
    shell:'''
    cutadapt {params.all_adapters} -o {output} {input}
    '''
