

rule dump_fastq_single:
    output:
        fastq='output/sample/fastq/single/{sample}.fastq',
        marker='output/sample/fastq/{sample}.SINGLE'
    params:
        out_dir='output/sample/fastq/single',
        sra_run=lambda wildcards: SAMPLES.loc[wildcards.sample]['Run']
    shell:'''
    fastq-dump {params.sra_run} --outdir {params.out_dir}
    touch {output.marker}
    '''

rule dump_fastq_paired:
    output:
        one='output/sample/fastq/paired/{sample}_1.fastq',
        two='output/sample/fastq/paired/{sample}_2.fastq',
        marker='output/sample/fastq/{sample}.PAIRED'
    params:
        out_dir='output/sample/fastq/paired',
        sra_run=lambda wildcards: SAMPLES.loc[wildcards.sample]['Run']
    shell:'''
    mkdir -p {params.out_dir}
    fastq-dump {params.sra_run} --split-files --outdir {params.out_dir}
    touch {output.marker}
    '''

    




