import argparse

# genequest imports
from genequest.io.parser import FastaReader
from ge




if __name__ == '__main__':

    # CLI arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("Required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument("-r" "--reads", type=str, required=True,
                        help="FASTA file containg reads")
    required.add_argument("-q" "--gene_query", type=str, required=True,
                        help="FASTA file gene of interest")

    # optional arguments
    # NOTE: these are not used yet. Used for alignment only
    optional.add_argument("-s" "--match_score", type=int, required=False, default=-5,
                          help="Score applied for matching nucleotide")
    optional.add_argument("-p" "--gap_penalty", type=int, required=False, default=0,
                          help="Penalty score applied for gaps")
    optional.add_argument("-m" "--mistmatch_penalty", type=int, required=False, default=-5,
                          help="Penalty score applied for nucleotide mismatch")

    args = parser.parse_args()

    # parsing fasta file containg reads and gene query
    reads = FastaReader(args.reads)
    gene_query = FastaReader(args.gene_query)

    # next is to execute the aligment
