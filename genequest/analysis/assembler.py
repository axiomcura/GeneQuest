#------------------------------
# aligment.py
# Erik Serrano
# erik.serrano@cuanschutz.edu
#
# Module containing functions for geneome aseembly by
# by using De bruijn graphs.
# paper explaining development of de brujin graphs: https://doi.org/10.1038/nbt.2023
#------------------------------
from collections import defaultdict
import itertools

def kmerize(sequence, k=3):
    """Breaks the sequences into fragments of size k

    Parameters
    ----------
    sequence : str
        DNA sequence
    k : int, optional
        framgment size, by default 3

    Returns
    -------
    dict
        containing kmer and n_instances as key value pairs
    """

    kmers = defaultdict(lambda: None)

    for i in range(len(sequence)):
        kmer_fragment = sequence[i:i+k]

        # checking if the size of the kmer is equal to 3 if not skip
        if len(kmer_fragment) != 3:
            continue

        # keep track how many times the fragment appears
        if kmer_fragment in kmers:
            kmers[kmer_fragment] += 1
        else:
            kmers[kmer_fragment]

    return kmers


def generate_edges_from_kmers(kmers):
    kmer_list = list(kmers.keys())



    # creating I J combinatino
    kmer_pairs = itertools.product(kmer_list, 2)

    # iterate through pairs and find overlaping kmers
    edges = {}
    for kmer1, kmer2 in kmer_pairs:
        # getting prefixes and sufixes
        if kmer1[1:] == kmer2[:-1]:
            edges.add((kmer1[1:], kmer2[:-1]))
        if kmer1[:-1] == kmer2[1:]:
            edges.add((kmer2[:-1], kmer2[:-1]))

    return edges


def run_de_brujin(sequences, k):
    """Builds a de brujin graph to solve assemble the sequence
    from reads.

    Parameters
    ----------
    sequences : str
        DNA sequence
    k : int
        DNA fragment size
    """

    pass
