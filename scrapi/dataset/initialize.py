import os
import argparse

from ..utils import fileio
from .Gene_Expression_Dataset import Gene_Expression_Dataset


def get_arguments():

    argparser = argparse.ArgumentParser(
        description="Initialize a SCRAP workspace")
    argparser.add_argument("--workspace_path", "-w",
                           help="Path to the workspace", default=None)
    argparser.add_argument("--h5_file_path", "-f",
                           help="Path to the H5 file", default=None)
    argparser.add_argument("--csv_file_path", "-c",
                           help="Path to a CSV file", default=None)
    argparser.add_argument("--mtx_directory_path", "-m",
                           help="Path to an mtx directory", default=None)
    argparser.add_argument("--sdt_file_path", "-s",
                           help="Path to an SDT file", default=None)

    args = argparser.parse_args()

    return args


def initialize():

    args = get_arguments()

    if args.workspace_path is None:
        workspace_path = os.path.join(os.getcwd(), "workspace")
    else:
        workspace_path = args.workspace_path

    if os.path.isfile(workspace_path):
        raise ValueError("Cannot use file as a workspace")

    if not os.path.isdir(workspace_path):
        os.makedirs(workspace_path, exist_ok=True)

    if args.h5_file_path is None and \
            args.csv_file_path is None and \
            args.mtx_directory_path is None and \
            args.sdt_file_path is None:
        raise ValueError("Must specify path to file with -f, -c, -m, or -s")

    raw_transcript_counts_file_path = os.path.join(
        workspace_path, "source_barcode_transcript_counts.sdt")

    if args.sdt_file_path is not None:
        raw_transcript_counts_file_path = args.sdt_file_path
    elif args.h5_file_path is not None:
        fileio.convert_h5_to_sdt(
            args.h5_file_path,
            raw_transcript_counts_file_path)
    elif args.csv_file_path is not None:
        raise NotImplementedError("Haven't implemented csv to sdt")
    else:
        sdt = fileio.load_mtx(
            os.path.join(args.mtx_directory_path, "matrix.mtx"),
            os.path.join(args.mtx_directory_path, "features.tsv"),
            os.path.join(args.mtx_directory_path, "barcodes.tsv")
        )
        sdt.save(raw_transcript_counts_file_path)

    seed_matrices_file_paths = [raw_transcript_counts_file_path]

    Gene_Expression_Dataset.initialize_dataset(
        workspace_path,
        seed_matrices_file_path=seed_matrices_file_paths
    )
