import unittest
import logging
import random
import numpy as np

# genequest imports
from genequest.io.parser import FastaEntry, FastaReader
from genequest.analysis.aligment import generate_scoring_matrix
#====================
# data generator functions
#====================
def generate_random_seq(n_nucleotides):
    dna_nuc = ["A", "T", "C", "G"]
    random_dna_seq = "".join(random.choices(dna_nuc, k=n_nucleotides))

    return random_dna_seq


# TODO: add readable logger
class TestFunctions(unittest.TestCase):
    logger = logging.getLogger(__name__)
    logging.basicConfig(format = '%(asctime)s %(module)s %(levelname)s: %(message)s',
                        datefmt = '%m/%d/%Y %I:%M:%S %p', level = logging.INFO)


    def test_FastaEntry(self):
        """ testing creating of FastaEntry Datatype"""

        # creating FastaEntry type
        entry = FastaEntry(
            header_id = "2S43D:04730:00438",
            scaffold_id = "2S43D",
            seq = "AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA",
            )

        test_data_type = str(entry)
        test_header_id = entry.header_id
        test_scaffold_id = entry.scaffold_id
        test_seq = entry.seq

        data_type = "FastaEntry(header_id='2S43D:04730:00438', scaffold_id='2S43D', seq='AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA')"
        header_id = "2S43D:04730:00438"
        scaffold_id = "2S43D"
        seq = "AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA"

        self.assertEqual(data_type, test_data_type)
        self.assertEqual(header_id, test_header_id)
        self.assertEqual(scaffold_id, test_scaffold_id)
        self.assertEqual(seq, test_seq)

    def test_FastaReader(self):
        """ Tests reading fasta file """
        reader = FastaReader("./test_data/test_fasta_seq.fasta")
        test_repr = str(repr(reader))
        repr_str = "FastaReader(filename=./test_data/test_fasta_seq.fasta, entries=9)"
        self.assertEqual(test_repr, repr_str)

    def test_FastaReader_iter(self):
        """ Tests FastaReader iteration"""
        entries_obj = [
                    FastaEntry(header_id='2S43D:03629:08794', scaffold_id='>2S43D', seq='TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT'),
                    FastaEntry(header_id='2S43D:08938:01257', scaffold_id='>2S43D', seq='GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCCAGGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTGGCTCCGACAGAAGCTACGGCATCCTATCGAGTGCACTGGGCTCGTGGTGGGAAGCAGGACA'),
                    FastaEntry(header_id='2S43D:05292:10188', scaffold_id='>2S43D', seq='GGGTGGTCTCCTTTACTTGTAACTTGTCCTAAGTCGTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCAGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTTGGCTTCCACAGAAGCTACGGCATCCTGATG'),
                    FastaEntry(header_id='2S43D:03619:08385', scaffold_id='>2S43D', seq='CAACAGGGTTTTGGAAATTTGCCCATTTGCATGGCGAAGACCACCTCTCTCTCTCTCATCGACCT'),
                    FastaEntry(header_id='2S43D:08782:12110', scaffold_id='>2S43D', seq='CCCCCCTCCTTTATTTTGTTGATTATTGAGTTTGGCATTCTGTTCTTGTGGCTCTCTTCTTTTGTTTCGTTTGAGGAATACTTCTTGGCTTTTTCTACTGGGCGTGAGTTTTCTTGGTCCTTGATTATTGGGTT'),
                    FastaEntry(header_id='2S43D:09644:04759', scaffold_id='>2S43D', seq='TAGGGCTGGAGGCTGGGGTAGTGTAACACATCCTACACGTGGCAGGCAGAGACAGGATGAACCTGATGACTTGGAGGCCAGCTTGATTTATGTAGCGAGTTTAGGTCATCCAAAGCTATACAGTGAGAACCTGTCTGAAAAAAACCAACAACCGAAATGAAAGAAAGAAAGAAAAGAA'),
                    FastaEntry(header_id='IDIDID:04730:00438', scaffold_id='>IDIDID', seq='AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA'),
                    FastaEntry(header_id='IDIDID:06986:00601', scaffold_id='>IDIDID', seq='GCGGGGGAGAGGATGGAGGCTTCTGAGTGAAAACGAGGAAGGGACTAAATTTCAAATGTAAATAAAGACAATATCTAATAAAAATAAAATAAAATTAATGGGGGGGGA'),
                    FastaEntry(header_id='IDIDID:02506:06718', scaffold_id='>IDIDID', seq='GGCTGACATGTATCTATGTTTAAATTAAGGTGCCCTGTCCTCCAATGTCTGCATTGCACTCAGAAGGGAGCCAAGTGCTGCTTGTAAAATGGAATCACTA')
                    ]

        entries = FastaReader("./test_data/test_fasta_seq.fasta")
        test_entries = entries.entries

        # merging entries together
        self.assertEqual(len(entries_obj), len(entries_obj))

        # iterating
        merged_entries = zip(entries_obj, test_entries)
        for entry, test_entry in merged_entries:
            self.assertEqual(str(entry), str(test_entry))


    def test_FastaReader_indexing(self):
        """ Tests FastaReader indexing"""
        entries_obj = [
                    FastaEntry(header_id='2S43D:03629:08794', scaffold_id='>2S43D', seq='TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT'),
                    FastaEntry(header_id='2S43D:08938:01257', scaffold_id='>2S43D', seq='GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCCAGGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTGGCTCCGACAGAAGCTACGGCATCCTATCGAGTGCACTGGGCTCGTGGTGGGAAGCAGGACA'),
                    FastaEntry(header_id='2S43D:05292:10188', scaffold_id='>2S43D', seq='GGGTGGTCTCCTTTACTTGTAACTTGTCCTAAGTCGTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCAGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTTGGCTTCCACAGAAGCTACGGCATCCTGATG'),
                    FastaEntry(header_id='2S43D:03619:08385', scaffold_id='>2S43D', seq='CAACAGGGTTTTGGAAATTTGCCCATTTGCATGGCGAAGACCACCTCTCTCTCTCTCATCGACCT'),
                    FastaEntry(header_id='2S43D:08782:12110', scaffold_id='>2S43D', seq='CCCCCCTCCTTTATTTTGTTGATTATTGAGTTTGGCATTCTGTTCTTGTGGCTCTCTTCTTTTGTTTCGTTTGAGGAATACTTCTTGGCTTTTTCTACTGGGCGTGAGTTTTCTTGGTCCTTGATTATTGGGTT'),
                    FastaEntry(header_id='2S43D:09644:04759', scaffold_id='>2S43D', seq='TAGGGCTGGAGGCTGGGGTAGTGTAACACATCCTACACGTGGCAGGCAGAGACAGGATGAACCTGATGACTTGGAGGCCAGCTTGATTTATGTAGCGAGTTTAGGTCATCCAAAGCTATACAGTGAGAACCTGTCTGAAAAAAACCAACAACCGAAATGAAAGAAAGAAAGAAAAGAA'),
                    FastaEntry(header_id='IDIDID:04730:00438', scaffold_id='>IDIDID', seq='AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA'),
                    FastaEntry(header_id='IDIDID:06986:00601', scaffold_id='>IDIDID', seq='GCGGGGGAGAGGATGGAGGCTTCTGAGTGAAAACGAGGAAGGGACTAAATTTCAAATGTAAATAAAGACAATATCTAATAAAAATAAAATAAAATTAATGGGGGGGGA'),
                    FastaEntry(header_id='IDIDID:02506:06718', scaffold_id='>IDIDID', seq='GGCTGACATGTATCTATGTTTAAATTAAGGTGCCCTGTCCTCCAATGTCTGCATTGCACTCAGAAGGGAGCCAAGTGCTGCTTGTAAAATGGAATCACTA')
                    ]

        entries = FastaReader("./test_data/test_fasta_seq.fasta")
        test_entries = entries.entries

        self.assertEqual(str(entries_obj[0]), str(test_entries[0]))
        self.assertEqual(str(entries_obj[1]), str(test_entries[1]))
        self.assertEqual(str(entries_obj[2]), str(test_entries[2]))


    def test_zero_matrix_builder(self):
        """ Builds a matrix of zero-filled matrix """
        test_contig = ""
        test_query = ""
        test_shape = (len(test_query) + 1, len(test_contig) + 1)

        # build the matrix (query size, contig size)

