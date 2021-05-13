

# rule collect_samples:
#     input:
#         'output/coverage/{sample}.{strand}.trim.bedgraph'
#     output:
#         'output/struct_score/{treatment}/{sample}.{strand}.{treatment}.bedgraph'
#     params:
#         out_dir=lambda wildcards: f'output/struct_score/{wildcards.treatment}'
#     shell:'''
#     mkdir -p {params.out_dir}
#     ln -sf {input} {output}
#     '''

# Currently using PAIRED where a wildcard for library layout should be to get
# the pipeline to run. Fix in the future. This works for now because all
# the sample we need for structure scores are paired end reads.

rule normalize_dsRNA_scores:
    input:
        'output/coverage/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        'output/struct_score/norm/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph'
    shell:'''
    python scripts/normalize.py {input} {output} --min_max_norm
    '''


rule arcsinh_ssRNA_scores:
    input:
        'output/coverage/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        'output/struct_score/asinh/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    shell:'''
    python scripts/normalize.py {input} {output} --arcsinh
    '''


rule dsRNA_intersect_ssRNA:
    conda:
        '../envs/bedtools.yml'
    input:
        dsRNA='output/struct_score/norm/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph',
        ssRNA='output/struct_score/asinh/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        'output/struct_score/intersect_dsRNA/{sample_dsRNA}.intersect.{sample_ssRNA}.{strand}.bed'
    shell:'''
    bedtools intersect -a {input.dsRNA} -b {input.ssRNA} > {output}
    '''


rule ssRNA_intersect_dsRNA:
    conda:
        '../envs/bedtools.yml'
    input:
        dsRNA='output/struct_score/norm/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph',
        ssRNA='output/struct_score/asinh/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        'output/struct_score/intersect_ssRNA/{sample_ssRNA}.intersect.{sample_dsRNA}.{strand}.bed'
    shell:'''
    bedtools intersect -b {input.dsRNA} -a {input.ssRNA} > {output}
    '''


rule concat_dsRNA_ssRNA_intersections:
    conda:
        '../envs/bedtools.yml'
    input:
        dsRNA='output/struct_score/intersect_dsRNA/{sample_dsRNA}.intersect.{sample_ssRNA}.{strand}.bed',
        ssRNA='output/struct_score/intersect_ssRNA/{sample_ssRNA}.intersect.{sample_dsRNA}.{strand}.bed'
    output:
        'output/struct_score/concat_intersect/{sample_dsRNA}.concat.{sample_ssRNA}.{strand}.tsv'
    shell:'''
    bedtools map -a {input.dsRNA} -b {input.ssRNA} -c 4 -o concat > {output}
    '''


rule non_overlaping_regions_dsRNA:
    # creating a combined score file first score should be the dsRNA always
    # and second is the ssRNA always since these regions did not overlap
    # we know the score for ssRNA to be 0 so can fill that in to have
    # the same format as the intersected and concat bed file
    conda:
        '../envs/bedtools.yml'
    input:
        dsRNA='output/struct_score/norm/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph',
        ssRNA='output/struct_score/asinh/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        pre_awk='output/struct_score/dsRNA_no_overlap/{sample_dsRNA}.{sample_ssRNA}.{strand}.preawk.bedgraph',
        final='output/struct_score/dsRNA_no_overlap/{sample_dsRNA}.nolap.{sample_ssRNA}.{strand}.tsv'
    shell:'''
    bedtools intersect -v -a {input.dsRNA} -b {input.ssRNA} > {output.pre_awk}
    awk '{{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "0"}}' > {output.final}
    '''


rule non_overlaping_regions_ssRNA:
    conda:
        '../envs/bedtools.yml'
    input:
        dsRNA='output/struct_score/norm/{sample_dsRNA}.PAIRED.{strand}.trim.bedgraph',
        ssRNA='output/struct_score/asinh/{sample_ssRNA}.PAIRED.{strand}.trim.bedgraph'
    output:
        pre_awk='output/struct_score/ssRNA_no_overlap/{sample_ssRNA}.{sample_dsRNA}.{strand}.preawk.bedgraph',
        final='output/struct_score/ssRNA_no_overlap/{sample_ssRNA}.nolap.{sample_dsRNA}.{strand}.tsv'
    shell:'''
    bedtools intersect -v -a {input.dsRNA} -b {input.ssRNA} > {output.pre_awk}
    awk '{{print $1 "\t" $2 "\t" $3 "\t" "0" "\t" $4}}' > {output.final}
    '''


rule concatenate_all_score_tsv_files:
    input:
        intersection='output/struct_score/concat_intersect/{sample_dsRNA}.concat.{sample_ssRNA}.{strand}.tsv',
        dsRNA_non_overlap='output/struct_score/dsRNA_no_overlap/{sample_dsRNA}.nolap.{sample_ssRNA}.{strand}.tsv',
        ssRNA_non_overlap='output/struct_score/ssRNA_no_overlap/{sample_ssRNA}.nolap.{sample_dsRNA}.{strand}.tsv'
    output:
        'output/struct_score/complete_bed/{sample_dsRNA}.{sample_ssRNA}.{strand}.complete.tsv'
    shell:'''
    cat {input.intersection} {input.dsRNA_non_overlap} \
    {input.ssRNA_non_overlap} > {output}
    '''


rule make_all_final_score_files:
    input:
        expand(
            expand(
                'output/struct_score/complete_bed/{sample_dsRNA}.{sample_ssRNA}.{strand}.complete.tsv',
                zip, sample_dsRNA=dsRNA_SAMPLES.index.tolist(), sample_ssRNA=ssRNA_SAMPLES.index.tolist(), allow_missing=True
            ),
        strand=['fwd', 'rev']
        )
    output:
        'output/struct_score/final_score_files.done'
    shell:'''
    touch {output}
    '''


