from dataclasses import dataclass
from .gfa import GfaSegments
from npstructures import HashTable

class GfaKmerIndex:
    """
    Simple lookup from kmers to gfa segment IDs
    """
    index: HashTable  # keys are kmers, values are numeric segment ids

    @classmethod
    def from_gfa_segments(cls, gfa_segments: GfaSegments):
        # for every segment, find kmers
        pass