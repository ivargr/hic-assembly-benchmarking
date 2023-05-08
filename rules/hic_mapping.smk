rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("{genome}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.2/bio/bwa/index"


rule map_hic:
    input:
        reads1=HiCReads.path() + "/reads1.fq.gz",
        reads2=HiCReads.path() + "/reads2.fq.gz",
        primary_assembly=HifiasmResultsWithExtraSplits.path() + "/{assembly}.fa",
        #bwa_index = multiext(HifiasmResultsWithExtraSplits.path() + "/{assembly}", ".fa.amb",".fa.ann",".fa.bwt",".fa.pac",".fa.sa")
    output:
        HifiasmResultsWithExtraSplits.path() + "/{assembly}.bam"
    conda:
        #"../envs/hic_mapping.yml"
        "../envs/picard.yml"
    params:
        out_dir=lambda wildcards, input, output: output[0].replace(wildcards.assembly + ".bam", "")
    shell:
        """
    arima_hic_mapping_pipeline/01_mapping_arima.sh {input.reads1} {input.reads2} {input.primary_assembly} {params.out_dir} {wildcards.assembly}
	#bwa mem -t {config[n_threads]} -5SPM {input.primary_assembly} \
	#{input.reads1} {input.reads2} \
	#| samtools view -buS - | samtools sort -n -O bam - \
	#| samtools fixmate -mr - -| samtools sort -O bam - | samtools markdup -rsS - {output}
        #"""



rule sort_hic_mapped_reads_by_name:
    input:
        "{file}.bam"
    output:
        "{file}.sorted_by_read_name.bam"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools sort -n {input} -o {output}
        """


