import os


rule run_quast:
    input:
        assembly=ScaffoldingResults.path() + "_scaffolds_final.fa",
        true_reference=ReferenceGenome.path(file_ending="") + "/haplotype0.fa"
    output:
        ScaffoldingResults.path() + "_quast_report/report.tsv"
    params:
        quast_dir=lambda wildcards, input, output: os.path.sep.join(output[0].split(os.path.sep)[:-2])
    conda:
        "envs/quast.yml"
    shell:
        """
        quast {input.assembly} -r {input.true_reference} -o {params.quast_dir}
        """
