import os


# A scaffolder giving correct results
rule truth_scaffolder:
    input:
        true_reference= ReferenceGenome.path(file_ending="") + "/haplotype0.fa"
    output:
        ScaffoldingResults.path(scaffolder="true_scaffolder") + "_scaffolds_final.fa"
    shell:
        """
        cp {input} {output}
        """


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
    threads:
        10000000  # hack: cannot be run in parallel because of temporary files
    params:
        edison_tmp_agp_file = lambda wildcards, input, output: input.assembly.split(os.path.sep)[-1].replace(".fa", "_assembly.agp"),
        edison_pdf_alignment = lambda wildcards, input, output: input.assembly.split(os.path.sep)[-1].replace(".fa", "_alignment.pdf"),
    shell:
        """
        rm -f {output} &&
        rm -f {params.edison_tmp_agp_file} &&
        echo {params.edison_tmp_agp_file} && 
        python edison/edit_distance.py -a {input.assembly} -r {input.true_reference} > {output} && 
        cat {output} && gio open {params.edison_pdf_alignment}
        """



# heatmap, edison
rule full_evaluation:
    input:
        ScaffoldingResults.path() + "_scaffolds_final_heatmap.png",
        ScaffoldingResults.path() + ".edison.txt",
        ScaffoldingResults.path() + "_quast_report/report.pdf",
    output:
        touch(ScaffoldingResults.path() + ".evaluation.txt")




rule accuracy:
    input:
        edison_results = ScaffoldingResults.path() + ".edison.txt"
    output:
        touch(ScaffolderAccuracy.path())
    run:
        with open(input[0]) as f:
            line = [l for l in f if "Accuracy:" in l][0]
            accuracy = float(line.split(": ")[1].replace("%", ""))

        with open(output[0], "w") as f:
            f.write(str(accuracy) + "\n")
