from typing import Literal
from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"
configfile: "config/plots.yaml"

workflow.use_conda = True

from snakehelp import parameters, result


@parameters
class GenomeBuild:
    genome_build: str = "sacCer3"

@parameters
class Individual:
    genome_build: GenomeBuild
    individual: str = "simulated"


@parameters
class SingleHaplotypeAndChromosomeHifiReads:
    individual: Individual
    depth: float = 1.0
    chromosome: str = "chrIV"
    haplotype: Literal[0, 1] = 0


@parameters
class HifiReads:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"]
    depth: float = 10


@parameters
class HiCReads:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"] = "small"
    hic: Literal["hic"] = "hic"
    n_reads: int = 100


@parameters
class HiCReadsHaplotype:
    reads: HiCReads
    haplotype: Literal[0, 1] = 0


@parameters
class HifiasmResults:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"] = "small"
    depth: float = 10
    n_reads: int = 500


@parameters
class HifiasmResultsWithExtraSplits:
    hifiasm_results: HifiasmResults
    source: Literal["assembled_from_hifi", "not_assembled"] = "assembled_from_hifi"
    extra_splits: int = 20


@parameters
class ScaffoldingResults:
    assembly_graph: HifiasmResultsWithExtraSplits
    scaffolder: Literal["yahs", "salsa2", "custom", "bnp_scaffolding", "true_scaffolder"] = "bnp_scaffolding"


@result
class ScaffolderAccuracy:
    scaffolding_results: ScaffoldingResults



@parameters
class PhasingResults:
    assembly_graph: HifiasmResultsWithExtraSplits
    phaser: Literal["gfase"]



#include: github("bioinf-benchmarking/mapping-benchmarking", "rules/reference_genome.smk", branch="master")
include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/reference_genome.smk"
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/read_simulation.smk", branch="master")
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/mason.smk", branch="master")
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/plotting.smk", branch="master")
#include: "/home/ivargry/dev/sync/mapping-benchmarking/rules/plotting.smk"
include: "rules/hifi_simulation.smk"
include: "rules/hic_simulation.smk"
include: "rules/hic_mapping.smk"
include: "rules/hifiasm.smk"
include: "rules/yahs.smk"
include: "rules/gfase.smk"
include: "rules/bnp_scaffolding.smk"
include: "rules/quast.smk"
include: "rules/evaluation.smk"
include: "rules/tests.smk"
include: "rules/real_data.smk"
include: "rules/salsa2.smk"


