

rule run_whatshap:
    input:
        truth = ReferenceGenome.path(file_ending="") + "/variants.vcf.gz",
        dipcall = HifiasmResults.path() + "/dipcall.dip.cleaned.vcf.gz"
    output:
        HifiasmResults.path() + "/whatshap.tsv"
    conda:
        "../envs/whatshap.yml"
    shell:
        """
        whatshap compare --names truth,sample --tsv-pairwise {output} {input} 
        """




rule make_hic_heatmap_for_scaffolds:
    input:
        scaffolds = ScaffoldingResults.path() + "_scaffolds_final.fa",
        hic_mapped_to_scaffolds = ScaffoldingResults.path() + "_scaffolds_final.bam",
    output:
        ScaffoldingResults.path() + "_scaffolds_final_heatmap.png"
    shell:
        """
        bnp_assembly heatmap {input} {output}
        """


rule run_edison:
    input:
        assembly = ScaffoldingResults.path() + "_scaffolds_final.fa",
        true_reference = ReferenceGenome.path(file_ending="") + "/haplotype0.fa"
    output:
        ScaffoldingResults.path() + ".edison.txt"
    conda: "../envs/edison.yml"
    shell:
        """
        python edison/edit_distance.py -a {input.assembly} -r {input.true_reference} > {output}
        """

