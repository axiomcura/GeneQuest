class FormatError(Exception):
    """Raised when files containing incorrect formats"""


class format_failure_msg(Exception):
    """allows failure messages from unittesting to become more readable"""


class KmerSizeError(Exception):
    """Raised if kmer size affects runtime"""


class NoNodesFoundError(Exception):
    """Raised if no nodes were found for assembly"""
