

rule map_hic:
    input:
        reads1=HiCReads.path() + "/reads1.fq.gz",
        reads2=HiCReads.path() + "/reads2.fq.gz",
        primary_assembly=HifiasmResults.path() + "/hifiasm.hic.p_ctg.fa"
    output:
        HifiasmResults.path() + "/hic.bam"
    conda:
        "../envs/hic_mapping.yml"
    params:
        out_dir=lambda wildcards, input, output: output[0].replace("hic.bam", "")
    shell:
        """
        arima_hic_mapping_pipeline/01_mapping_arima.sh {input} {params.out_dir}
        """



rule sort_hic_mapped_reads_by_name:
    input:
        HifiasmResults.path() + "/hic.bam"
    output:
        HifiasmResults.path() + "/hic.sorted_by_read_name.bam"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools sort -n {input} -o {output}
        """


