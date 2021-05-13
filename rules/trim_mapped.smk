rule trim_mapped_reads:
    conda:
        '../envs/samtools.yml'
    input:
        'output/bowtie2/{sample}.{layout}.bam'
    output:
        temp('output/trim_mapped/{sample}.{layout}.trim.bam')
    shell:'''
    samtools view -q 2 -bhu {input} > {output}
    '''


rule index_trimmed_bam:
    conda:
        '../envs/samtools.yml'
    input:
        'output/bowtie2/{sample}.{layout}.bam'
    output:
        'output/bowtie2/{sample}.{layout}.bai'
    threads: 4
    shell:'''
    samtools index -b --threads {threads} {input} {output}
    '''


rule trim_mitochonrial_reads:
    conda:
        '../envs/samtools.yml'
    input:
        bam='output/trim_mapped/{sample}.{layout}.trim.bam',
        bai='output/bowtie2/{sample}.{layout}.bai'
    output:
        temp('output/trim_mapped/{sample}.{layout}.trim.nomito.bam')
    shell:'''
    samtools idxstats {input.bam} | cut -f 1 | grep -v MT | xargs samtools view -b {input.bam} > {output}
    '''


rule trim_rRNA_reads:
    # trimming based on this post
    # https://www.biostars.org/p/159959/
    conda:
        '../envs/bbtools.yml'
    input:
        reads='output/trim_mapped/{sample}.{layout}.trim.nomito.bam',
        rRNAs='data/ribokmers.fa.gz'
    output:
        temp('output/trim_mapped/{sample}.{layout}.trim.nomito.norRna.bam')
    shell:'''
    bbduk.sh in={input.reads} out={output} k=31 ref={input.rRNAs}
    '''