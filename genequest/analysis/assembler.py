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
from genequest.io_handler.parser import FastaReader
from genequest.common.errors import KmerSizeError, NoNodesFoundError
from genequest.io_handler.gene_io import save_contigs


class Node:
    """Class Node to represent a vertex in the de bruijn graph"""

    def __init__(self, label):
        self.label = label
        self.indegree = 0
        self.outdegree = 0


class Edge:
    def __init__(self, label):
        self.label = label


def create_bruijn_graph(reads, k=3) -> tuple:
    """Generate a

    Parameters
    ----------
    reads : _type_
        _description_
    k : int, optional
        sequence fragment length, by default 3
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


def run_de_bruijn(sequences: FastaReader, k: str, save=False) -> defaultdict:
    """Builds a de brujin graph to solve assemble the sequence
    from reads.

    Parameters
    ----------
    sequences : str
        DNA sequence
    k : int
        DNA fragment size
    save : bool, optional
        saves contigs data into a pickle file. Default = False. If set to true
        a contigs_data/ directory will be created contining all the pickle files.
        To open the files, import the genequest.io_handler.gene_io and import the
        load_contigs() function load in the contigs data. This is to prevent re-
        running the assembly algorithm

    Raises
    ------
    KmerSizeError
        Raised if the kmer size is larger than the length of the shortest sequence read

    NoNodesFoundError
        raised when no starting nodes are found for assembly
    """

    # kmer length checking (cannot be larger than the smallest reads)
    min_length = min([seq.end_pos for seq in sequences])
    if not k < min_length:
        raise KmerSizeError(
            f"k size is larger than the smallest read must be smaller than {min_length}"
        )

    # group sequences by scaffold
    grouped_sequences = sequences.group_by_scaffold()

    generated_contigs = defaultdict(None)
    for scaffold, reads in grouped_sequences.items():

        # generating graph only using 10000 reads per scaffold
        nodes, edges = create_bruijn_graph(reads[:10000], k)

        # finding starting nodes
        starting_edges = []
        for kmer in nodes.keys():
            if nodes[kmer].indegree == 0:
                starting_edges.append(kmer)

        if len(starting_edges) == 0:
            raise NoNodesFoundError(
                "No starting nodes were found, try increasing k size"
            )

        # assembly
        labeled_contigs = defaultdict(None)
        for idx, starting_edge in enumerate(starting_edges):
            d_edges = deepcopy(edges)
            contig = generate_contig(starting_edge, d_edges)
            labeled_contigs[f"conting{idx+1}"] = contig

        generated_contigs[scaffold] = labeled_contigs

    if save is True:
        save_contigs(generated_contigs)

    return generated_contigs
