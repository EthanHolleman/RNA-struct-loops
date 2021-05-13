
ALL_ADAPTER_PARAMS = ' '.join([f'-a {seq}' for seq in ADAPTERS['Sequence'].tolist()])


def paired_end_adapter_params(*args, **kwargs):
    fwd_adapters = ADAPTERS.loc[ADAPTERS['direction'] == 'fwd']['Sequence'].tolist()
    rev_adapters = ADAPTERS.loc[ADAPTERS['direction'] == 'rev']['Sequence'].tolist()
    all_adapters = [f'-a {seq}' for seq in fwd_adapters] + [f'-A {seq}' for seq in rev_adapters]
    return ' '.join(all_adapters)


rule trim_all_adapters_single:
    conda:
        '../envs/cutadapt.yml'
    input:
        fq = 'output/sample/fastq/single/{sample}.fastq.gz'
    output:
        'output/trimmed/single/{sample}.trim.fastq.gz'
    params:
        all_adapters=lambda wildcards: ALL_ADAPTER_PARAMS
    shell:'''
    cutadapt {params.all_adapters} -o {output} {input.fq}
    '''


rule trim_all_adapters_paired:
    conda:
        '../envs/cutadapt.yml'
    input:
        mate_1='output/sample/fastq/paired/{sample}_1.fastq.gz',
        mate_2='output/sample/fastq/paired/{sample}_2.fastq.gz',
        qc='output/qc/{sample}'
    output:
        mate_1='output/trimmed/paired/{sample}_1.trim.fastq.gz',
        mate_2='output/trimmed/paired/{sample}_2.trim.fastq.gz'
    params:
        all_adapters=lambda wildcards: paired_end_adapter_params()
    threads: 8
    shell:'''
    mkdir -p output/trimmed/paired/
    cutadapt {params.all_adapters} -j {threads}  \
    -o {output.mate_1} -p {output.mate_2} {input.mate_1} {input.mate_2} 
    '''
