# README

# GeneQuest

`GeneQuest` is a command line (CLI) tool that searches for genes in unassembled genome. `GeneSearch` utilizes unassembled genomes and gene targets in the`FASTA` format and produces  two output `FASTA` files:

- `ALLESES.fasta` - File that contains the largest constructed contig along with the aligned gene query.
- `ALLESES.aln` - Tab delimited file that contains the alignment of sequence

## Requirements

`GeneQuests` dependencies are written in the `requirements.txt` file: 

```
memory_profiler==0.58.0
numpy==1.22.2
pandas==1.4.1
packaging==21.3
```

## Installation

Next the installation assumes that you have `python3` and `pip` 

First, the source code must be downloaded through GitHub but using `git clone`:

```
git clone https://github.com/axiomcura/GeneSearch.git

```

Once downloaded go into the `GeneSearch` folder and use `pip`  package manager to install all dependencies by typing:

```
pip install -e .
```

To check if the installation was completed, execute the `run_tests.py` script to tests `GeneQuest`'s dependencies and functionality .

```
python ./tests/run_tests.py
```

*NOTE*: If an error is caught, a detailed message should appear indicating what the error is and which testing step if failed in. 

Otherwise, a complete message should is displayed if all tests we successfully completed. 

## Usage

### Documentation

One can access the `GeneQuest` documentation by using the `-h` or `--help` flag to display it on the terminal. 

```
python genequest.py --help

usage: genequest.py [-h] -r--reads R__READS -q GENE_QUERY [-k KMER_SIZE] [-ms MATCH_SCORE] [-gs GAP_SCORE] [-mis MISTMATCH_SCORE]

options:
  -h, --help            show this help message and exit

Required arguments:
  -r--reads R__READS    FASTA file containg reads
  -q GENE_QUERY, --gene_query GENE_QUERY
                        FASTA file gene of interest

Optional arguments:
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        size of the sequence fragments
  -ms MATCH_SCORE, --match_score MATCH_SCORE
                        Score applied for matching nucleotide
  -gs GAP_SCORE, --gap_score GAP_SCORE
                        Penalty score applied for gaps
  -mis MISTMATCH_SCORE, --mistmatch_score MISTMATCH_SCORE
                        Penalty score applied for nucleotide mismatch
```

### Basic Use

Below is an example using two FASTA file where one contains the unassembled reads `READS.fasta` and gene of interest sequence `QUERY.fasta` 

First we uses the sequences as inputs for `GeneQuest` to use:

```
python genequest.py -r READS.fasta -q QUERY.fasta
```

Then messages were start to print out indicating which stages the analysis is currently in:

```
MESSAGE: parsing FASTA file 
-- number of reads:  32342

MESSAGE: Running genome assembly 
MESSAGE: Running local alimgnet with query gene
-- query gene and contig2_<Scaffold_id>
-- Algnment score: 12321

MESSAGE: Saving Files 
 -- saved: absolute/path/to/file/ALLELES.fasta
 -- saved: absolute/path/to/file/ALLELES.aln

Complete!
```

## Additional Notes

- A `ValueError` is thrown out when then `-k` parameter has a value higher than the length of the reads
- `FormatError` is raised if indicating that that it is an unsupported file time or
- Users can change the penalty scores of the local alignment by changing the `--gap_score`, `--match_score` and `--mismatch_score`. However, this requires prior knowledge on your sequences dictating how strict or lenient you want the scoring to be. If the gene happens to be vary conserved as more stricter scoring parameters are required. If the genome is distantly related from the organism where the gene sequence was obtained, then a more lenient score is required.