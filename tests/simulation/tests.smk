


rule test_hifiasm:
    input:
        HifiasmResults.from_flat_params(
            genome_build="sacCer3",
            individual="simulated",
            depth=1,
            dataset_size="small",
            n_reads=100
        ).file_path() + "/whatshap.tsv"
    output:
        touch("test_hifiasm")


rule test_yahs:
    input:
        ScaffoldingResults.from_flat_params(scaffolder="yahs").file_path() + "_scaffolds_final.fa"
    output:
        touch("test_yahs")
