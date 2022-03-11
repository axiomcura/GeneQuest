import sys
import argparse

# genequest imports

sys.path.append("../GeneQuest")
import genequest
from genequest.io_handler.parser import FastaReader
# from genequest.io.parser import FastaReader

if __name__ == '__main__':

    # CLI arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    # required arguments
    required.add_argument("-r" "--reads", type=str, required=True,
                        help="FASTA file containg reads")
    required.add_argument("-q", "--gene_query", type=str, required=True,
                        help="FASTA file gene of interest")

    # optional arguments
    # NOTE: these are not used yet. Used for alignment only
    optional.add_argument("-k", "--kmer_size", type=int, required=False, default=3,
                          help="size of the sequence fragments")
    optional.add_argument("-ms", "--match_score", type=int, required=False, default=-5,
                          help="Score applied for matching nucleotide")
    optional.add_argument("-gs", "--gap_score", type=int, required=False, default=0,
                          help="Penalty score applied for gaps")
    optional.add_argument("-mis", "--mistmatch_score", type=int, required=False, default=-5,
                          help="Penalty score applied for nucleotide mismatch")

    args = parser.parse_args()

    # parsing fasta file containg reads and gene query
    # -- printing sequence information
    reads = FastaReader(args.reads)
    gene_query = FastaReader(args.gene_query)

    print("Loaded Sequence information:")
    print("query sequence length: {}".format(len(gene_query)))
    print("total number of reads: {}".format(len(reads)))

    # assemble the genome
    # -- group the geomes based on scaffold_id
    # -- apply the run_de_brujin() function

    # next step is the alignment
    # -- iterate all generated contigs into the run_local_alignment(query, contig)
        # returns alignment score a positional alignment data

    # write out the results 


