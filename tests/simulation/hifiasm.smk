


rule run_hifiasm_with_hic_reads:
    input:
        hic1 = HiCReads.path() + "/reads1.fq.gz",
        hic2 = HiCReads.path() + "/reads2.fq.gz",
        hifi = HifiReads.path() + ".fq"

    output:
        multiext(HifiasmResults.path() + "/out.asm", ".hic.hap1.p_ctg.gfa", ".hic.hap2.p_ctg.gfa")
    conda:
        "envs/hifiasm.yml"
    shell:
        "hifiasm -o {output} -t6 --h1 {input.hic1} --h2 {input.hic2} {input.hifi}"