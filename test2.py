from bionumpy.kmers import TwoBitHash
from bionumpy.files import bnp_open
k = 31
hasher = TwoBitHash(k)
chunks = bnp_open("hic_chr8_1.fastq.gz")
for chunk in chunks:
    print(chunk)
    print(hasher.get_kmer_hashes(chunk.sequence))




