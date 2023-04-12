import pytest
from npstructures import HashTable, HashSet, RaggedArray
from hic_phaser.gfa_kmer_index import GfaKmerIndex
import numpy as np


@pytest.fixture
def gfa_kmer_index():
    index = HashTable([
        123, 5, 20, 10, 1
    ],
    [
        1, 2, 3, 4, 5
    ])
    in_index = HashSet([123, 5, 20, 10, 1])

    return GfaKmerIndex(index, in_index)


def test_gfa_kmer_index(gfa_kmer_index):
    kmers = RaggedArray(
        [
            [200, 4, 10, 1000, 1],
            [5, 1],
            [2],
            [123, 10, 1, 1, 10],
            [25, 10]
        ]
    )

    hits = gfa_kmer_index.get_as_ragged_array(kmers)
    assert hits == RaggedArray([
        [4, 5],
        [2, 5],
        [],
        [1, 4, 5, 5, 4],
        [4]
    ])
