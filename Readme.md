

## Snakemake pipeline for testing HiC phasing and scaffolding on simulated data

Overview:

* This pipeline can simulate a diploid individual with snps, indels and optinally structural variation on a reference genome of choice. Specify the individual in `config/config.yaml`. 
* Hifi + HiC reads can be simulated from this genome
* Hifiasm is run using the simulated reads
* Yahs and other tools can be run to scaffold the assembly graph created by Hifiasm
* The resulting assembly is evaluated by Quast


### Installation

1. Clone this repo. 
2. You unfortunately need to manually install pbsim3. This can be done by running `./install_pbsim.sh`. GFAse also needs to be installed manually (but this tool is not used by the main pipeline).
2. Install Python requirements: `pip install -r requirements.txt`
3. Test that the pipeline works: `snakemake --use-conda -F test_hifiasm` (should finish without errors in ~2 minutes)


### Run scaffolding with yahs and quality assesment with quast

See the rule `test_quast` in `rules/test.smk`. 

```snakemake
snakemake --use-conda -F test_quast
```



### Make a test plot
```snakemake
snakemake --use-conda plots/test.png
```

### Missing functionality / todo
* Genome simulation is hardcoded in config and the same for all simulations. Should be possible to configure for a given run
* Heterozygosity is fixed, should be tunable. Every variant is now default heterozygous. 