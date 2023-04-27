import os


rule run_gfase:
    input:
        unitig_graph=HifiasmResults.path() + "/hifiasm.hic.p_utg.gfa",
        sorted_hic_reads=HifiasmResults.path() + "/hic.sorted_by_read_name.bam"
    output:
        PhasingResults.path(phaser="gfase") + "/phased.fa"
    params:
        out_dir=lambda wildcards, input, output: os.path.sep.join(output[0].split(os.path.sep)[:-1])
    shell:
        """
        rm -rf {params.out_dir} && 
        # remove old index file if has been created before
        rm -f {input.unitig_graph}i && 
        phase_contacts_with_monte_carlo -t 4 -i {input.sorted_hic_reads} -g {input.unitig_graph} -o {params.out_dir} --use_homology
        """