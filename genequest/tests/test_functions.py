from cmath import exp
import unittest
import logging
import random
import numpy as np

# genequest imports
from genequest.io_handler.parser import FastaEntry, FastaReader
from genequest.analysis.alignment import generate_scoring_matrix, score_alignment, trace_back

# ====================
# data generator functions
# ====================
def generate_random_seq(n_nucleotides):
    dna_nuc = ["A", "T", "C", "G"]
    random_dna_seq = "".join(random.choices(dna_nuc, k=n_nucleotides))

    return random_dna_seq


# TODO: add readable logger
class ParserFunctions(unittest.TestCase):

    # creating a stdout logger
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="%(asctime)s %(module)s %(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.INFO,
    )

    def test_FastaEntry(self):
        """testing creating of FastaEntry Datatype"""

        # creating FastaEntry type
        entry = FastaEntry(
            header_id="2S43D:04730:00438",
            scaffold_id="2S43D",
            seq="AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA",
        )

        test_data_type = str(entry)
        test_header_id = entry.header_id
        test_scaffold_id = entry.scaffold_id
        test_seq = entry.seq

        data_type = "FastaEntry(header_id='2S43D:04730:00438', scaffold_id='2S43D', seq='AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA', beg_pos='0', end_pos='52')"
        header_id = "2S43D:04730:00438"
        scaffold_id = "2S43D"
        seq = "AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA"

        try:
            self.assertEqual(data_type, test_data_type)
            self.assertEqual(header_id, test_header_id)
            self.assertEqual(scaffold_id, test_scaffold_id)
            self.assertEqual(seq, test_seq)
        except:
            self.logger.error("FastEntry Test: FAILED")
            self.fail("Please read Exception above")

        self.logger.info("FastaEntry Test: PASSED")

    def test_FastaReader(self):
        """Tests reading fasta file"""
        reader = FastaReader("./test_data/test_fasta_seq.fasta")
        test_repr = str(repr(reader))
        repr_str = "FastaReader(filename=./test_data/test_fasta_seq.fasta, entries=9)"

        try:
            self.assertEqual(test_repr, repr_str)
        except:
            self.logger.error("FastaReader Test: FAILED")
            self.fail("FastaEntry elements in FastaReader are not the same")

        self.logger.info("FastaReader Test: PASSED")

    def test_FastaReader_indexing(self):
        """Tests FastaReader indexing"""
        entries_obj = [
            FastaEntry(
                header_id="2S43D:03629:08794",
                scaffold_id="2S43D",
                seq="TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT",
            ),
            FastaEntry(
                header_id="2S43D:08938:01257",
                scaffold_id="2S43D",
                seq="GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCCAGGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTGGCTCCGACAGAAGCTACGGCATCCTATCGAGTGCACTGGGCTCGTGGTGGGAAGCAGGACA",
            ),
            FastaEntry(
                header_id="2S43D:05292:10188",
                scaffold_id="2S43D",
                seq="GGGTGGTCTCCTTTACTTGTAACTTGTCCTAAGTCGTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCAGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTTGGCTTCCACAGAAGCTACGGCATCCTGATG",
            ),
            FastaEntry(
                header_id="2S43D:03619:08385",
                scaffold_id="2S43D",
                seq="CAACAGGGTTTTGGAAATTTGCCCATTTGCATGGCGAAGACCACCTCTCTCTCTCTCATCGACCT",
            ),
            FastaEntry(
                header_id="2S43D:08782:12110",
                scaffold_id="2S43D",
                seq="CCCCCCTCCTTTATTTTGTTGATTATTGAGTTTGGCATTCTGTTCTTGTGGCTCTCTTCTTTTGTTTCGTTTGAGGAATACTTCTTGGCTTTTTCTACTGGGCGTGAGTTTTCTTGGTCCTTGATTATTGGGTT",
            ),
            FastaEntry(
                header_id="2S43D:09644:04759",
                scaffold_id="2S43D",
                seq="TAGGGCTGGAGGCTGGGGTAGTGTAACACATCCTACACGTGGCAGGCAGAGACAGGATGAACCTGATGACTTGGAGGCCAGCTTGATTTATGTAGCGAGTTTAGGTCATCCAAAGCTATACAGTGAGAACCTGTCTGAAAAAAACCAACAACCGAAATGAAAGAAAGAAAGAAAAGAA",
            ),
            FastaEntry(
                header_id="IDIDID:04730:00438",
                scaffold_id="IDIDID",
                seq="AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA",
            ),
            FastaEntry(
                header_id="IDIDID:06986:00601",
                scaffold_id="IDIDID",
                seq="GCGGGGGAGAGGATGGAGGCTTCTGAGTGAAAACGAGGAAGGGACTAAATTTCAAATGTAAATAAAGACAATATCTAATAAAAATAAAATAAAATTAATGGGGGGGGA",
            ),
            FastaEntry(
                header_id="IDIDID:02506:06718",
                scaffold_id="IDIDID",
                seq="GGCTGACATGTATCTATGTTTAAATTAAGGTGCCCTGTCCTCCAATGTCTGCATTGCACTCAGAAGGGAGCCAAGTGCTGCTTGTAAAATGGAATCACTA",
            ),
        ]

        entries = FastaReader("./test_data/test_fasta_seq.fasta")
        test_entries = entries.entries

        try:
            self.assertEqual(str(entries_obj[0]), str(test_entries[0]))
            self.assertEqual(str(entries_obj[1]), str(test_entries[1]))
            self.assertEqual(str(entries_obj[2]), str(test_entries[2]))
        except:
            self.logger.error("FastaReader Iteration Test: FAILED")
            self.fail("Failed to index FastaReader object")
        self.logger.info("FastaReader Iteration Indexing Test: PASSED")

    def test_FastaReader_iter(self):
        """Tests FastaReader iteration"""
        entries_obj = [
            FastaEntry(
                header_id="2S43D:03629:08794",
                scaffold_id="2S43D",
                seq="TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTGTTTT",
            ),
            FastaEntry(
                header_id="2S43D:08938:01257",
                scaffold_id="2S43D",
                seq="GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAGTCGTTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCCAGGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTGGCTCCGACAGAAGCTACGGCATCCTATCGAGTGCACTGGGCTCGTGGTGGGAAGCAGGACA",
            ),
            FastaEntry(
                header_id="2S43D:05292:10188",
                scaffold_id="2S43D",
                seq="GGGTGGTCTCCTTTACTTGTAACTTGTCCTAAGTCGTTTCTTTAGCCCATGGTGTTGGTGGGGTTCACAGAAACACCCAGAGTTCACCTGAGCCTTTAACCAATCCCAGCCAGGAGCCAGAGCCCAGGCACAGGTGCAGGACCACGGCAGGCCCAGTATTTGGCTTCCACAGAAGCTACGGCATCCTGATG",
            ),
            FastaEntry(
                header_id="2S43D:03619:08385",
                scaffold_id="2S43D",
                seq="CAACAGGGTTTTGGAAATTTGCCCATTTGCATGGCGAAGACCACCTCTCTCTCTCTCATCGACCT",
            ),
            FastaEntry(
                header_id="2S43D:08782:12110",
                scaffold_id="2S43D",
                seq="CCCCCCTCCTTTATTTTGTTGATTATTGAGTTTGGCATTCTGTTCTTGTGGCTCTCTTCTTTTGTTTCGTTTGAGGAATACTTCTTGGCTTTTTCTACTGGGCGTGAGTTTTCTTGGTCCTTGATTATTGGGTT",
            ),
            FastaEntry(
                header_id="2S43D:09644:04759",
                scaffold_id="2S43D",
                seq="TAGGGCTGGAGGCTGGGGTAGTGTAACACATCCTACACGTGGCAGGCAGAGACAGGATGAACCTGATGACTTGGAGGCCAGCTTGATTTATGTAGCGAGTTTAGGTCATCCAAAGCTATACAGTGAGAACCTGTCTGAAAAAAACCAACAACCGAAATGAAAGAAAGAAAGAAAAGAA",
            ),
            FastaEntry(
                header_id="IDIDID:04730:00438",
                scaffold_id="IDIDID",
                seq="AAGATTTTTGTTGTAGATAGTGATAAACCAGCTACCCCATCCTAGTCTTAAA",
            ),
            FastaEntry(
                header_id="IDIDID:06986:00601",
                scaffold_id="IDIDID",
                seq="GCGGGGGAGAGGATGGAGGCTTCTGAGTGAAAACGAGGAAGGGACTAAATTTCAAATGTAAATAAAGACAATATCTAATAAAAATAAAATAAAATTAATGGGGGGGGA",
            ),
            FastaEntry(
                header_id="IDIDID:02506:06718",
                scaffold_id="IDIDID",
                seq="GGCTGACATGTATCTATGTTTAAATTAAGGTGCCCTGTCCTCCAATGTCTGCATTGCACTCAGAAGGGAGCCAAGTGCTGCTTGTAAAATGGAATCACTA",
            ),
        ]

        entries = FastaReader("./test_data/test_fasta_seq.fasta")
        test_entries = entries.entries

        # testing length property
        try:
            self.assertEqual(len(entries_obj), len(entries_obj))
        except:
            self.logger.error("FastaReader Iteration Test: FAILED")
            self.fail("Cannot identify the length FastaReader object")
        self.logger.info("FastaReader Iteration Length attribute Test: PASSED")

        # checking iteration properties in FastaReader Object
        try:
            merged_entries = zip(entries_obj, test_entries)
            for entry, test_entry in merged_entries:
                self.assertEqual(str(entry), str(test_entry))
        except:
            self.logger.error("FastaReader Iteration Test: FAILED")
            self.fail("Iteration could to be conducted with FastaReader object")
        self.logger.info("FastaReader Iteration Test: PASSED")

    class AlignmentTest(unittest.TestCase):
        """testes all functions and cases when assembling"""

        # creating a stdout logger
        logger = logging.getLogger(__name__)
        logging.basicConfig(
            format="%(asctime)s %(module)s %(levelname)s: %(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
            level=logging.INFO,
        )

    def test_zero_matrix_builder(self):
        """Builds a matrix of zero-filled matrix"""
        test_contig = "CCCCTACATGTTGTTATAGACAATCAGTGGAAACCCAGTGCCAGACGATG"
        test_query = "CAATCAGTGGAAACCCAGTG"

        score_matrix = generate_scoring_matrix(len(test_query), len(test_contig))
        expected_size = (20, 50)

        # testing if matrices are equal
        try:
            self.assertEqual(expected_size, score_matrix.shape)
        except:
            self.logger.error("Zero Matrix Builder Test: FAILED")
            self.fail("Generated matrix does not have the expected shape")

        try:
            self.assertEqual(0, score_matrix.max())
        except:
            self.logger.error("Zero Matrix Builder Test: FAILED")
            self.fail("Non-zero value found in scoring matrix")

        self.logger.info("Matrix Builder: PASSED")

    # TODO: create test
    def test_alignment_score(self):
        """Produce a message if the alignment score"""

        # contig and gene as inputs
        contig = "ACGTA"
        query = "ACG"

        expected_alignment_score = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 10.0, 5.0, 0.0, 0.0, 10.0],
                [0.0, 5.0, 20.0, 15.0, 10.0, 5.0],
                [0.0, 0.0, 15.0, 30.0, 25.0, 20.0],
            ]
        )

        # construct scoring matrix
        test_alignment_score = score_alignment(contig, query)

        # checking if all values in the matrix are the same.
        # -- if one value is different, then value_check == False, else == True
        value_check = np.all(expected_alignment_score == test_alignment_score)

        try:
            self.assertEqual(True, value_check)
        except:
            self.logger.error("Local Alignment Test: FAILED")
            self.fail("Alignment score did not match with the expected scores")

        self.logger.info("Local Alignment Test: PASSED")

    def test_traceback(self):
        """Tests for tracing back along the local alignment scoring
        matrix to find the best alignment
        """

        # expected values
        contig = "ACGTA"
        query = "ACG"

        expected_alignment_score = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 10.0, 5.0, 0.0, 0.0, 10.0],
                [0.0, 5.0, 20.0, 15.0, 10.0, 5.0],
                [0.0, 0.0, 15.0, 30.0, 25.0, 20.0],
            ]
        )

        expected_position, expected_scores = (
            [(3, 3), (2, 2), (1, 1), (0, 1)],
            [30.0, 20.0, 10.0, 0.0],
        )
        
        expected_results = None        

        # construct scoring matrix
        test_alignment_matrix = score_alignment(contig, query)
        test_positions, test_scores = trace_back(query, contig, test_alignment_matrix)

        # checking for equal alignment score matrices
        equal_check = np.all(expected_alignment_score == test_alignment_matrix)
        try:
            self.assertEqual(True, equal_check)
        except:
            self.logger.error("Traceback Test: FAILED")
            self.fail("Alignment score did not match with the expected scores")
        self.logger.info("Traceback Test -- scoring matrix: PASSED")


        # check if traceback positions is the same as the expected
        try:
            self.assertEqual(expected_position, test_positions)
        except:
            self.logger.error("Traceback Test - positions: FAILED")
            self.fail("Trace back algorithm failed to obtain accurate positions")
        self.logger.info("Traceback Test - positions: PASSED")

        # testing if the obtains scores are equal to the expected
        try:
            self.assertEqual(expected_scores, test_scores)
        except:
            self.logger.error("Traceback Test - scores: FAILED")
            self.fail("Obtained scores are not equal")
        self.logger.info("Traceback Test - scores: PASSED")

        # check if the same sequence contig position is obtained
