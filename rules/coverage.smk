

rule sort_trimmed_bam:
    conda:
        '../envs/samtools.yml'
    input:
        'output/trim_mapped/{sample}.{layout}.trim.nomito.norRna.bam'
    output:
        temp('output/trim_mapped/{sample}.{layout}.trim.nomito.norRna.sorted.bam')
    shell:'''
    samtools sort {input} {output}
    '''

rule calculate_coverage_fwd:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/trim_mapped/{sample}.{layout}.trim.nomito.norRna.sorted.bam'
    output:
        temp('output/coverage/{sample}.{layout}.fwd.bedgraph')
    shell:'''
    bedtools genomecov -i -bg strand "+" {input}
    '''


rule calculate_coverage_rev:
    conda:
        '../envs/bedtools.yml'
    input:
        'output/trim_mapped/{sample}.{layout}.trim.nomito.norRna.sorted.bam'
    output:
        temp('output/coverage/{sample}.{layout}.rev.bedgraph')
    shell:'''
    bedtools genomecov -i -bg strand "-" {input}
    '''


rule remove_low_coverage_regions:
    input:
        'output/coverage/{sample}.{layout}.{strand}.bedgraph'
    output:
        'output/coverage/{sample}.{layout}.{strand}.trim.bedgraph'
    shell:'''
    awk '$4 > 8' {input} > {output}
    '''





