
rule run_yahs:
    input:
        contigs=HifiasmResults.path() + "/hifiasm.hic.p_ctg.fa",
        contigs_index=HifiasmResults.path() + "/hifiasm.hic.p_ctg.fa.fai",
        hic_to_contig_mappings=HifiasmResults.path() + "/p_ctg.bam",
    output:
        ScaffoldingResults.path(scaffolder="yahs") + "_scaffolds_final.fa"
    conda:
        "../envs/yahs.yml"
    params:
        out_prefix = lambda wildcards, input, output: output[0].replace("_scaffolds_final.fa", "")
    shell:
        """
        yahs -o {params.out_prefix} {input.contigs} {input.hic_to_contig_mappings} 
        """


