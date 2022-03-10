from memory_profiler import profile
from genequest.io.parser import FastaReader

@profile
def usecase_1():
    """ Test and measure memory usage"""
    reader = FastaReader("./test_data/test_fasta_seq.fasta")
    entries = reader.entries
    for entry in reader:
        data = entry


    del reader

if __name__ == "__main__":
    r = usecase_1()


