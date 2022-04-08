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

def kmerize(entries, k=3):
    """Breaks the sequences into fragments of size k

    Parameters
    ----------
    sequence : list
        DNA sequences reads
    k : int, optional
        framgment size, by default 3

    Returns
    -------
    dict
        containing kmer and n_instances as key value pairs
    """

    kmers = defaultdict(lambda: 1)

    for entry in entries:
        for i in range(len(entry.seq)-k+1):
            kmer_fragment = entry.seq[i:i+k]

            # keep track how many times the fragment appears
            if kmer_fragment in kmers:
                kmers[kmer_fragment] += 1
            else:
                kmers[kmer_fragment]

    return kmers


def generate_nodes(kmers):
    """generate nodes based on given edges

    Parameters
    ----------
    edges : list
        list of edges
    """
    nodes = set()
    for kmer in kmers:
        node = kmer[:len(kmer)-1]
        nodes.add(node)

    return list(nodes)



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



def run_de_brujin(sequences, k=3):
    """Builds a de brujin graph to solve assemble the sequence
    from reads.

    Parameters
    ----------
    sequences : str
        DNA sequence
    k : int
        DNA fragment size
    """
    # group sequences
    grouped_sequences = sequences.group_by_scaffold()

    # # # breaking the sequences into fragments
    kmerized_seqs = []
    for s_id, entries in grouped_sequences.items():
        kmers = kmerize(entries, k)
        result = (s_id, kmers)
        kmerized_seqs.append(result)

    grouped_nodes = []
    for s_id, edges in kmerized_seqs:
        edge_list = list(edges.keys())
        nodes = generate_nodes(edge_list)
        result = (s_id, nodes)
        grouped_nodes.append(result)

    

    return grouped_nodes
