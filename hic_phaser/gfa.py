from dataclasses import dataclass
from typing import List
import logging

@dataclass
class GfaSegment:
    id: str
    sequence: str

    @classmethod
    def from_gfa_line(cls, line):
        l = line.split()
        assert l[0] == "S"
        sequence = l[2]
        return cls(l[1], sequence)


@dataclass
class GfaSegments:
    segments: List[GfaSegment]

    @classmethod
    def from_gfa_file(cls, file_name):
        with open(file_name) as f:
            segments = [GfaSegment.from_gfa_line(l)
                        for l in f if l.startswith("S")]

        return cls(segments)


def index(gfa_file_name, k=31):
    data = bnp_open(gfa_file_name, chunk_size=10000000)

    # first find all unique kmers
    hasher = TwoBitHash(k)
    kmers = np.concatenate(
        [hasher.get_kmer_hashes(chunk.sequence).ravel() for chunk in data]
    )
    logging.info("Counting kmers")
    unique_kmers = np.unique(kmers)
    logging.info("%d unique kmers" % len(unique_kmers))
    kmer_counter = Counter(unique_kmers)
    kmer_counter.count(kmers)

    kmers_with_frequency_1 = unique_kmers[kmer_counter[unique_kmers] == 1]
    logging.info("%d kmers have frequency 1" % len(kmers_with_frequency_1))

    # make a lookup from kmer to node segment
    kmer_index = HashTable(kmers_with_frequency_1, 0)
    data = bnp_open(gfa_file_name, chunk_size=10000000)
    for chunk in data:
        kmers = hasher.get_kmer_hashes(chunk.sequence)
        for sequence_id, sequence_kmers in enumerate(kmers):
            kmer_index[sequence_kmers[np.where(kmer_counter[sequence_kmers] == 1)[0]]] = sequence_id

    logging.info("Made hashmap from unique kmers to sequence ids")

    if out_file_name is not None:
        to_file(kmer_index, out_file_name)
        logging.info("Wrote kmer index to %s" % out_file_name)
    else:
        return kmer_index

