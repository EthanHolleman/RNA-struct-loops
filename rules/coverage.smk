
rule calculate_coverage_fwd:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/trim_mapped/{sample}.trim.nomito.norRna.bam'
    output:
        'output/coverage/{sample}.fwd.bedgraph'
    shell:'''
    bedtools genomecov -i -bg strand "+" {input}
    '''


rule calculate_coverage_rev:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/trim_mapped/{sample}.trim.nomito.norRna.bam'
    output:
        'output/coverage/{sample}.rev.bedgraph'
    shell:'''
    bedtools genomecov -i -bg strand "-" {input}
    '''


rule remove_low_coverage_regions:
    input:
        'output/coverage/{sample}.{strand}.bedgraph'
    output:
        'output/coverage/{sample}.{strand}.trim.bedgraph'
    shell:'''
    awk '$4 > 9' {input} > {output}
    '''





