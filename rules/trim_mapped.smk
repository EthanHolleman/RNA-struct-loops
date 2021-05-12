rule trim_mapped_reads:
    conda:
        '../envs/samtools.yml'
    input:
        'output/bowtie2/{sample}.bam'
    output:
        'output/trim_mapped/{sample}.trim.bam'
    shell:'''
    samtools view -q 2 -bhu {input} > {output}
    '''


rule trim_mitochonrial_reads:
    conda:
        '../envs/samtools.yml'
    input:
        'output/trim_mapped/{sample}.trim.bam'
    output:
        'output/trim_mapped/{sample}.trim.nomito.bam'
    shell:'''
    samtools idxstats {input} | cut -f 1 | grep -v MT | xargs samtools view -b {input} > {output}
    '''


rule trim_rRNA_reads:
    # trimming based on this post
    # https://www.biostars.org/p/159959/
    conda:
        '../envs/bbtools.yml'
    input:
        reads='output/trim_mapped/{sample}.trim.nomito.bam',
        rRNAs='data/ribokmers.fa.gz'
    output:
        'output/trim_mapped/{sample}.trim.nomito.norRna.bam'
    shell:'''
    bbduk.sh in={input.reads} out={output} k=31 ref={input.rRNAs}
    '''