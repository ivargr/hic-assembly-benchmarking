

rule map_hic:
    input:
        reads1=HiCReads.path() + "/reads1.fq.gz",
        reads2=HiCReads.path() + "/reads2.fq.gz",
        primary_assembly=HifiasmResults.path() + "/hifiasm.hic.{graph}.fa"
    output:
        HifiasmResults.path() + "/{graph}.bam"
    conda:
        "../envs/hic_mapping.yml"
    params:
        out_dir=lambda wildcards, input, output: output[0].replace(wildcards.graph + ".bam", "")
    shell:
        """
	bwa mem -t {config[n_threads]} -5SPM {input.primary_assembly} \
	{input.reads1} {input.reads2} \
	|samtools view -buS - |samtools sort -n -O bam - \
	|samtools fixmate -mr - -|samtools sort -O bam - |samtools markdup -rsS - {output}
        """



rule sort_hic_mapped_reads_by_name:
    input:
        HifiasmResults.path() + "/{name,\w+}.bam"
    output:
        HifiasmResults.path() + "/{name,\w+}.sorted_by_read_name.bam"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools sort -n {input} -o {output}
        """


