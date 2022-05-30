import logging
import tqdm
import numpy as np

from hic_phaser.gfa_kmer_index import GfaKmerIndex
logging.basicConfig(level=logging.INFO)
import typer
from bionumpy.files import bnp_open
from bionumpy.kmers import TwoBitHash
from shared_memory_wrapper import from_file, to_file
from bionumpy.kmers import fast_hash

app = typer.Typer()


def main():
    app()


@app.command()
def index_gfa(gfa: str, out_file_name: str, k: int=31):
    index = GfaKmerIndex.from_gfa(gfa, k)
    to_file(index, out_file_name)
    logging.info("Wrote kmer index to %s" % out_file_name)


@app.command()
def map_reads(kmer_index: str, hic1: str, hic2: str, out_file_name, k: int=31, limit_to_n_reads: int=-1):
    kmer_index = from_file(kmer_index)
    hasher = TwoBitHash(k)
    read_chunks = [bnp_open(hic1), bnp_open(hic2)]
    # using Two-Bit-hash wich is a tiny bit faster but may have some weird bugs now
    #kmer_chunks = [(hasher.get_kmer_hashes(chunk.sequence) for chunk in chunks)
    #               for chunks in read_chunks]

    # using "fast-hash" which should be safe to use
    kmer_chunks = [(fast_hash(chunk.sequence, k) for chunk in chunks)
                   for chunks in read_chunks]


    mapped_node_ids = [[], []]

    for readset in (0, 1):
        n_reads_processed = 0
        for chunk_idx, chunk in enumerate(kmer_chunks[readset]):
            logging.info("Processing kmer chunk %d" % chunk_idx)
            kmer_hits = np.unique(kmer_index.get_as_ragged_array(chunk), axis=-1)  # only store unique kmer hits for each read
            mapped_node_ids[readset].append(kmer_hits)
            n_reads_processed += len(kmer_hits)

            if limit_to_n_reads > 0 and n_reads_processed > limit_to_n_reads:
                logging.info("Stopping at %d reads" % n_reads_processed)
                break

    mapped_node_ids = [np.concatenate(m)[0:limit_to_n_reads] for m in mapped_node_ids]
    to_file(mapped_node_ids, out_file_name)


@app.command()
def get_link_counts(mapped_reads: str, out_file_name: str, max_node_id: int=100):
    mappings = from_file(mapped_reads)
    print(mappings)
    counts = np.zeros((max_node_id, max_node_id))

    for mapped1, mapped2 in tqdm.tqdm(zip(*mappings), total=len(mappings[0])):
        for node1 in mapped1:
            for node2 in mapped2:
                counts[node1, node2] += 1

    logging.info("Total links: %d" % np.sum(counts))
    to_file(counts, out_file_name)




