
rule remove_duplicates_single:
    conda:
        '../envs/bbtools.yml'
    input:
        'output/trimmed/single/{sample}.trim.fastq'
    output:
        'output/clumped/single/{sample}.trim.clumped.fastq'
    shell:'''
    clumpify.sh in={input} out={output} groups=16
    '''


rule remove_duplicates_paired:
    conda:
        '../envs/bbtools.yml'
    input:
        'output/trimmed/paired/{sample}_{mate}.trim.fastq'
    output:
        'output/clumped/paired/{sample}_{mate}.trim.clumped.fastq'
    shell:'''
    clumpify.sh in={input} out={output} groups=16
    '''

