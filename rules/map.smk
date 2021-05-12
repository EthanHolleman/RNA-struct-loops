PAIRED_STRUCT_MAPS = expand(
    'output/bowtie2/{sample}.{layout}.bam',
    zip, sample=STRUCT_SCORE_SAMPLES['Run'].tolist(),
    layout=STRUCT_SCORE_SAMPLES['LibraryLayout']
)

rule map_and_convert_struct_reads:
    input:
        PAIRED_STRUCT_MAPS
    output:
        'output/bowtie2/map_and_convert_struct_reads.done'
    shell:'''
    touch {output}
    '''

rule download_hg19_bt_index:
    output:
        index_dir=directory('output/bowtie2/hg19_index'),
        downloaded_zip='output/bowtie2/hg19_index/hg19.zip',
    shell:'''
    mkdir --parents {output.index_dir}
    curl https://genome-idx.s3.amazonaws.com/bt/hg19.zip -o {output.downloaded_zip}
    unzip {output.downloaded_zip} -d {output.index_dir}
    '''



rule map_reads_paired:
    conda: 
        '../envs/bowtie2.yml'
    input:
        mate_1='output/trimmed/paired/{sample}_1.trim.fastq',
        mate_2='output/trimmed/paired/{sample}_2.trim.fastq',
        bt_index='output/bowtie2/hg19_index',
    output:
        'output/bowtie2/{sample}.PAIRED.sam'
    threads: 12
    shell:'''
    mkdir -p output/bowtie2
    bowtie2 -q -x {input.bt_index}/hg19 -1 {input.mate_1} -2 {input.mate_2} \
    -p {threads} -S {output}
    '''


rule sam_to_bam:
    conda:
        '../envs/samtools.yml'
    input:
        'output/bowtie2/{sample}.{layout}.sam'
    output:
        'output/bowtie2/{sample}.{layout}.bam'
    shell:'''
    samtools view -bhSu -o {output} {input}
    '''







