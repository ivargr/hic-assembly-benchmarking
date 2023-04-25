

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
