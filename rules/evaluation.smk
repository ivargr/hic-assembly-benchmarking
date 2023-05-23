import os


# A scaffolder giving correct results
rule truth_scaffolder:
    input:
        true_reference= ReferenceGenome.path(file_ending="") + "/haplotype0.fa"
    output:
        ScaffoldingResults.path(scaffolder="true_scaffolder") + "/scaffolds.fa"
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
        scaffolds = ScaffoldingResults.path() + "/scaffolds.fa",
        agp = ScaffoldingResults.path() + "/scaffolds.agp",
        hic_mapped_to_scaffolds = ScaffoldingResults.path() + "/scaffolds.sorted_by_read_name.bam",
        _ = ScaffoldingResults.path() + "/scaffolds.fa.fai",
    output:
        ScaffoldingResults.path() + "/heatmap.png"
    shell:
        """
        bnp_assembly heatmap {input.scaffolds} {input.hic_mapped_to_scaffolds} {input.agp} {output}
        """


rule run_edison:
    input:
        assembly = ScaffoldingResults.path() + "/scaffolds.fa",
        true_reference = ReferenceGenome.path(file_ending="") + "/haplotype0.fa"
    output:
        txt_report = ScaffoldingResults.path() + "/edison.txt",
        alignment_viz = ScaffoldingResults.path() + "/alignment.pdf",
        agp = ScaffoldingResults.path() + "/assembly.agp",
    conda: "../envs/edison.yml"
    threads:
        10000000  # hack: cannot be run in parallel because of temporary files
    params:
        edison_agp_file = lambda wildcards, input, output: input.assembly.split(os.path.sep)[-1].replace(".fa", "_assembly.agp"),
        edison_pdf_alignment = lambda wildcards, input, output: input.assembly.split(os.path.sep)[-1].replace(".fa", "_alignment.pdf"),
    shell:
        """
        rm -f {output} &&
        rm -f {params.edison_agp_file} &&
        echo {params.edison_agp_file} && 
        python edison/edit_distance.py -a {input.assembly} -r {input.true_reference} > {output.txt_report} && 
        mv {params.edison_pdf_alignment} {output.alignment_viz} &&
        mv {params.edison_agp_file} {output.agp} &&
        cat {output.txt_report}
        """
        #&& gio open {params.edison_pdf_alignment}



# heatmap, edison
rule full_evaluation:
    input:
        ScaffoldingResults.path() + "/heatmap.png",
        ScaffoldingResults.path() + "/edison.txt",
        ScaffoldingResults.path() + "/quast_report/report.pdf",
    output:
        touch(ScaffoldingResults.path() + "/evaluation.txt")




rule accuracy:
    input:
        edison_results = ScaffoldingResults.path() + "/edison.txt"
    output:
        touch(ScaffolderAccuracy.path())
    run:
        with open(input[0]) as f:
            line = [l for l in f if "Accuracy:" in l][0]
            accuracy = float(line.split(": ")[1].replace("%", ""))

        with open(output[0], "w") as f:
            f.write(str(accuracy) + "\n")
