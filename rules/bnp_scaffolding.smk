



rule run_bnp_scaffolding:
    input:
        contigs=HifiasmResults.path() + "/hifiasm.hic.p_ctg.fa",
        contigs_index=HifiasmResults.path() + "/hifiasm.hic.p_ctg.fa.fai",
        hic_to_contig_mappings=HifiasmResults.path() + "/hic.sorted_by_read_name.bam",
    output:
        ScaffoldingResults.path(scaffolder="bnp_scaffolding") + "_scaffolds_final.fa"
    shell:
        """
        bnp_assembly {input.contigs} {input.hic_to_contig_mappings} > {output}
        """








