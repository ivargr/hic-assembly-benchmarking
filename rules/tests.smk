


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


rule test_bnp_scaffolding:
    input:
        ScaffoldingResults.from_flat_params(scaffolder="bnp_scaffolding", depth=2, n_reads=50000).file_path() + "_scaffolds_final.fa"
    output:
        touch("test_bnp_scaffolding")


rule test_quast:
    input:
        ScaffoldingResults.from_flat_params(scaffolder="yahs", depth=2, n_reads=50000).file_path() + "_quast_report/report.tsv",
        ScaffoldingResults.from_flat_params(scaffolder="bnp_scaffolding", depth=2, n_reads=50000).file_path() + "_quast_report/report.tsv"
    output:
        touch("test_quast")


rule test_pbsim:
    input:
        SingleHaplotypeAndChromosomeHifiReads.from_flat_params().file_path() + "_0001.sam",
    output:
        touch("test_pbsim")


