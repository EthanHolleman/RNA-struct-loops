rule run_fastqc_paired:
    conda:
        '../envs/fastqc.yml'
    input:
        mate_1='output/sample/fastq/paired/{sample}_1.fastq.gz',
        mate_2='output/sample/fastq/paired/{sample}_2.fastq.gz'
    output:
        directory('output/qc/{sample}')
    threads: 4
    shell:'''
    mkdir -p {output}
    fastqc -t {threads} -o {output} {input}
    '''

