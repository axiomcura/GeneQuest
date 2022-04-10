# ------------------------------
# assembler.py
# Erik Serrano
# erik.serrano@cuanschutz.edu
#
# Module containing functions for genome assembly by
# by using De bruijn graphs.
# paper explaining development of de bruijn graphs: https://doi.org/10.1038/nbt.2023
# ------------------------------
from copy import deepcopy
from collections import defaultdict


class Node:
    """Class Node to represent a vertex in the de bruijn graph"""

    # TODO: add reads id, start and ending
    def __init__(self, label):
        self.label = label
        self.indegree = 0
        self.outdegree = 0


class Edge:
    # TODO: add reads id, start and ending
    def __init__(self, label):
        self.label = label


def create_bruijn_graph(reads, k=3) -> tuple:
    """Generate a

    Parameters
    ----------
    reads : _type_
        _description_
    k : int, optional
        _description_, by default 3
    """

    edges = dict()
    nodes = dict()
    for read in reads:
        for i in range(len(read.seq) - k + 1):
            edge1 = read.seq[i : i + k]
            edge2 = read.seq[i + 1 : i + k + 1]
            if edge1 in edges.keys():
                nodes[edge1].outdegree += 1
                edges[edge1] += [Edge(edge2)]
            else:
                nodes[edge1] = Node(edge1)
                nodes[edge1].outdegree += 1
                edges[edge1] = [Edge(edge2)]

            if edge2 in edges.keys():
                nodes[edge2].indegree += 1
            else:
                nodes[edge2] = Node(edge2)
                nodes[edge2].indegree += 1
                edges[edge2] = []
            i += 1

    return (nodes, edges)


def generate_contig(starting_edge, edges):
    """
    Assembles genomes through eulerian walks
    """
    # identifying all possible start nodes that have
    # Conducting eulerian walks
    contig = starting_edge
    current = starting_edge
    while len(edges[current]) > 0:
        next = edges[current][0]
        del edges[current][0]
        contig += next.label[-1]
        current = next.label

    return contig


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
    # group sequences
    grouped_sequences = sequences.group_by_scaffold()

    generated_contigs = defaultdict(lambda: None)
    for scafold, reads in grouped_sequences.items():
        nodes, edges = create_bruijn_graph(reads[:1000], k)

        # finding starting nodes
        starting_edges = []
        for kmer in nodes.keys():
            if nodes[kmer].indegree == 0:
                starting_edges.append(kmer)

        labeled_contigs = defaultdict(None)
        for idx, starting_edge in enumerate(starting_edges):
            d_edges = deepcopy(edges)
            contig = generate_contig(starting_edge, d_edges)
            labeled_contigs[f"conting{idx+1}"] = contig

        generated_contigs[scafold] = labeled_contigs

    return generated_contigs
