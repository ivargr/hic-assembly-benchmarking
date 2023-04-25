from typing import Literal
from snakemake.utils import min_version
min_version("6.0")

configfile: "config/config.yaml"

from snakehelp import parameters


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
    depth: int = 1
    chromosome: str = "chrIV"
    haplotype: Literal[0, 1] = 0


@parameters
class HifiReads:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"]
    depth: int = 5


@parameters
class HiCReads:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"] = "small"
    hic: Literal["hic"] = "hic"
    n_reads: int = 100


@parameters
class HifiasmResults:
    individual: Individual
    dataset_size: Literal["small", "medium", "big"] = "small"
    depth: int = 5
    n_reads: int = 500


@parameters
class ScaffoldingResults:
    assembly_graph: HifiasmResults
    scaffolder: Literal["yahs", "custom"]



include: github("bioinf-benchmarking/mapping-benchmarking", "rules/reference_genome.smk", branch="master")
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/read_simulation.smk", branch="master")
include: github("bioinf-benchmarking/mapping-benchmarking", "rules/mason.smk", branch="master")
include: "rules/hifi_simulation.smk"
include: "rules/hic_simulation.smk"
include: "rules/hic_mapping.smk"
include: "rules/hifiasm.smk"
include: "rules/yahs.smk"
include: "rules/quast.smk"
include: "rules/evaluation.smk"
include: "rules/tests.smk"




