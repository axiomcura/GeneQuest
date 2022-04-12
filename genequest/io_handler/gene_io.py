import os
import glob
import pickle
from pathlib import Path

# genequest imports
import genequest.common.utils as utils


def save_contigs(contigs_data: dict, outfile="contigs_data"):
    """Serializes contigs data into a pickle file

    Parameters
    ----------
    contigs_data : dict
        dictionary containing contig data
    outfile : str, optional
        name of the compressed file, by default "contigs_data"
    """
    path_check = os.path.join(os.getcwd(), f"{outfile}.pickle")
    if os.path.exists(path_check) is False:
        e_msg = "File exists, use other outname"
        FileExistsError(e_msg)

    # saving contigs data
    path = Path("./contig_data")
    path.mkdir(exist_ok=True)

    uid = utils.generate_unique_id()
    outname = f"{path.stem}/{outfile}_{uid}.pickle"
    with open(outname, "wb") as outfile:
        pickle.dump(contigs_data, outfile)


def load_contigs(contigs_data_path=None) -> dict:
    """Loads in contig data from pickle files

    Parameters
    ----------
    contigs_data_path : str, NoneType, optional
        Path to pickle file containing contig data. if None, the latest
        contig file will be open. , by default None

    Returns
    -------
    dict
        contig data

    Raises
    ------
    FileNotFoundError
        unable to find specified path to contig file
    RuntimeError
        if more than latest files were found
    """

    # if contigs_data_path is None, by default is searches for the latest
    if contigs_data_path is None:

        # searching for latest file
        search_path = "contig_data/*.pickle"
        pickle_files = glob.glob(search_path)
        latest_file = max(pickle_files, key=os.path.getctime)
        if len(latest_file) == 1:
            raise RuntimeError("More than 2 latest files found")

        # opening latest file and returning dictionary of contig data
        print(f"loading contig data: {latest_file}")
        with open(latest_file, "rb") as infile:
            contigs_data = pickle.load(infile)
            return contigs_data

    # loading pickle file and returning a dictionary
    path_check = os.path.exists(contigs_data_path)
    if path_check is False:
        e_msg = "unable to find pickle file"
        raise FileNotFoundError(e_msg)

    with open(contigs_data_path, "rb") as infile:
        contigs_data = pickle.load(infile)
        return contigs_data
