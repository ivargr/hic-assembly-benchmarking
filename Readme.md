
# HicPhaser

## Installation

1. Install BioNumpy by cloning the latest version of BioNumpy
2. Install hic_phaser:

```bash
git clone https://github.com/ivargr/hicphaser.git
cd hicphaser
python3 -m pip install -e .
```


## How to use
Typical work-flow:

### Step 1: Index kmers in a GFA
This step creates a lookup (in practice a [HashTable](https://github.com/knutdrand/npstructures)) from
unique kmers to the ID of the node that this kmer is found in. Nodes represent segments in the GFA, and are given numeric
ids from 0 to the number of segments (based on the order).

The GFA (for now, because of lack of full GFA-support in bionumpy) needs to contain only segments, so start by making a new GFA
with only segments:

```bash
grep "^S" my.gfa > only_segments.gfa
```

Then create the index:
```bash
hicphaser index-gfa --k 31 only_segments.gfa index
```

An index.npz file should be created after reunning the above


### Step 2: Map HiC-read kmers to this index
For every read, kmers are extracted, and reads with kmers are assigned a node ID based on where the kmers match
(a read may have kmers matching multiple nodes, this is ignored for now and only one node is chosen for simplicity).

```bash
hicphaser map-reads --limit-to-n-reads 500000 index hic_1.fastq.gz hic_2.fastq.gz mapping
```

Set `--limit-to-n-reads` to -1 to map all reads (a number other than -1 is useful only for testing/debugging to make things faster).

A file `mapping.npz` is created. This file contains two numpy arrays, one for each HiC input read set. In each array, the position represents the read ID and the value at that position is the node that the read "mapped" to.


### Step 3: Make a link-matrix
Given the two arrays from the previous step, we can count how many times two nodes in the graph are "linked" by a read, i.e. how often a pair of reads maps to the two nodes.

```bash
hicphaser get-link-counts mapping link_counts
```

We now have a matrix representing the link counts, which can be analyzed further:

```python
from shared_memory_wrapper import from_file
counts = from_file("link_counts")
print(counts)
```

