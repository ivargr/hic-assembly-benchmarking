import bionumpy as bnp
import numpy as np
import logging
np.random.seed(1)

logging.basicConfig(level=logging.INFO)

extra_splits = int(snakemake.wildcards.extra_splits)


contigs = bnp.open(snakemake.input[0]).read()

new_contig_names = []
new_contig_sequences = []
new_contig_id = 0


total_size = np.sum(contigs.sequence.shape[1])

weights = contigs.sequence.shape[1]/total_size
logging.info(f"Using weights {weights}")
splits_at_contig = np.bincount(np.sort(
    np.random.choice(np.arange(len(contigs)), extra_splits, p=weights)), minlength=len(contigs))

logging.info(f"Will split contigs n times: {splits_at_contig}")


def random_spaced_locations(start, stop, n, min_space=1000):
    assert stop > min_space
    min_space = min(min_space, stop-start-1)
    candidates = np.arange(start, stop, min_space)
    assert len(candidates) >= n
    np.random.shuffle(candidates)
    return candidates[0:n]


min_contig_size = 15000

for contig_id, n_splits in enumerate(splits_at_contig):
    if n_splits == 0:
        new_contig_names.append(f"contig{new_contig_id}")
        new_contig_sequences.append(contigs.sequence[contig_id])
        new_contig_id += 1
        continue

    # make n_splits new contigs (between the splits)
    old_contig_sequence = contigs.sequence[contig_id]
    split_positions = np.sort(random_spaced_locations(min_contig_size, len(old_contig_sequence)-min_contig_size,
                                                      n_splits, min_space=min_contig_size))
    split_positions = np.insert(split_positions, 0, 0)
    split_positions = np.append(split_positions, len(old_contig_sequence))
    print(split_positions)
    logging.info(f"Splitting contig {contig_id} between {split_positions}")
    for start, end in zip(split_positions[0:-1], split_positions[1:]):
        logging.info(f"New contig at old contig {contig_id} between {start} and {end}")
        new_contig_names.append(f"contig{new_contig_id}")
        new_contig_sequences.append(contigs.sequence[contig_id][start:end])
        new_contig_id += 1

new_fasta = bnp.datatypes.SequenceEntry.from_entry_tuples(
    zip(new_contig_names, new_contig_sequences)
)
logging.info(f"Ended up with {len(new_fasta)} contigs")

with bnp.open(snakemake.output[0], "w") as f:
    f.write(new_fasta)
