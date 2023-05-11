

rule simulate_hic_for_haplotype:
    input:
        reference = ReferenceGenome.path(file_ending="") + "/haplotype{haplotype}.fa",
    output:
        reads1 = HiCReadsHaplotype.path() + "/1.fq.gz",
        reads2= HiCReadsHaplotype.path() + "/2.fq.gz",
    conda:
        "../envs/sim3c.yml"
    params:
        abundance_profile = lambda wildcards, input, output: "/".join(output.reads1.split("/")[:-1]) + "/profile.tsv",
        tmp_output = lambda wildcards, input, output: output.reads1.replace(".fq.gz", ".tmp")
    shell:
        """
        rm -f {params.abundance_profile} && 
        sim3C --seed 123 --dist uniform -n {wildcards.n_reads} -l 150 -e NlaIII --insert-mean 1000 -m hic {input.reference} {params.tmp_output} && 
        seqtk seq -1 {params.tmp_output} | gzip -c > {output.reads1} && 
        seqtk seq -2 {params.tmp_output} | gzip -c > {output.reads2} 
        """



rule merge_hic_haplotype_reads:
    input:
        haplotype0 = HiCReadsHaplotype.path(haplotype=0) + "/{pair}.fq.gz",
        haplotype1 = HiCReadsHaplotype.path(haplotype=1) + "/{pair}.fq.gz",
    output:
        HiCReads.path() + "/reads{pair}.fq.gz",
    shell:
        "zcat {input} | gzip -c > {output}"
