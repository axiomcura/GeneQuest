from collections import defaultdict
from genequest.common.errors import FormatError

class FastaEntry:
    """ Class object that contains the FASTA entry information"""
    __slot__ = ["header_id", "scafold_id", "seq", "beg_pos", "end_pos"]

    def __init__(self, header_id, scaffold_id, seq):
        self.header_id = header_id
        self.scaffold_id = scaffold_id
        self.seq = seq
        self.beg_pos = 0
        self.end_pos = len(seq)

    def __str__(self):
        """String representation of the data type"""
        str_rep = f"FastaEntry(header_id='{self.header_id}', scaffold_id='{self.scaffold_id}', seq='{self.seq}', beg_pos='{self.beg_pos}', end_pos='{self.end_pos}')"
        return str_rep

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.seq)


class FastaReader():
    """ Class that is used for parsing FASTA files. Contains function to conduct
    simple edits.

    parameters
    ----------
    filename:  str
        path that leads to FASTA file


    methods
    -------
    reverse(fasta_entry):method return reversed sequences along with its position
    """

    def __init__(self, filename):
        self.__index = 0
        self.filename = filename
        self.entries = None
        self.n_entries = None
        self.ids = None
        self.scaffold_ids = None

        # populating attributes
        self.__parse_fasta()


    # TODO: make it only support FastaEntry datatypes
    @staticmethod
    def reverse(entry):
        """ Returns the reverse position of the data"""
        if not isinstance(entry, FastaEntry):
            TypeError("entry must be a FastaEntry type, you provided {}".format(type(entry)))
        pass

    # TODO: grouping function
    def group_by_scaffold(self):
        """ Groups all entries based on a scaffold"""
        grouped_entries = defaultdict(list)
        for entry in self.entries:
            grouped_entries[entry.scaffold_id].append(entry)
        
        return grouped_entries

    # TODO: useful function to search
    def search(self, header_id):
        """Obtains sequence read with given header_id

        Parameters
        ----------
        header_id : str
            FASTA header id

        Returns
        -------
        str
            sequence associated with given fasta header_id

        Raises
        ------
        ValueError
            raised if the provided header_id does not exist
        """
        pass


    # -----------------
    # Private functions
    # -----------------
    def __parse_fasta(self):
        """Parses Fasta files and returns a list of FastaEntry objects

        Parameters
        ----------
        filename : str
            path to FASTA file

        Returns
        -------
        list
            list of FastaEntry objects

        Raises
        ------
        FormatError
            raised if input file is not a FASTA file
        """
        with open(self.filename, "r") as f:
            data = f.readlines()
            if not data[0].startswith(">"):
                raise FormatError("Invalid FASTA file")

            raw_entries = [tuple(data[i:i+2]) for i in range(0, len(data), 2)]

            entries = []
            scaffold_set = set()
            for header, seq in raw_entries:

                # removing unwanted formating
                header_id = header.strip().replace("\n", "").replace(">", "")
                scaffold_id = header.split(':')[0].replace(">", "")
                seq = seq.strip().replace("\n", "")

                # Collecting data and converting it into FastaEntry type
                entry = FastaEntry(header_id, scaffold_id, seq)
                entries.append(entry)

                scaffold_set.add(scaffold_id)


        self.entries = entries
        self.n_entries = len(entries)
        self.scaffold_ids = list(scaffold_set)
        self.n_scaffolds = len(list(scaffold_set))


    #---------------------
    # class attributes
    #---------------------
    # TODO: data type transformation functions
    # Data Type transformation support
    # NOTE: to list
    def to_list(self):
       """ converts FASTA entries into python list"""
       pass

    def to_pandas(self):
        """Converts FASTA entries into pandas DataFrame"""
        pass

    def to_dict():
        """ Converts FASTA entries into python dictionary"""
        pass

    # Allow python functionallity support
    # -- indexing and iterating support
    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.entries[self.__index]
        except IndexError:
            self.__index = 0
            raise StopIteration

        self.__index += 1
        return result

    def __getitem__(self, val):
        return self.entries[val]

    # -- printing support
    def __repr__(self):
        return f"FastaReader(filename={self.filename}, entries={self.n_entries})"

    def __str__(self):
        return f"FastaReader: Filename: '{self.filename}' has {self.n_entries} entries"





