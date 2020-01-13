import os
import argparse

from .Gene_Expression_Dataset import Gene_Expression_Dataset
from .Data_Mode import Data_Mode
from .Transformation_Method import Transformation_Method
from .Normalization_Method import Normalization_Method


def get_arguments():

    argparser = argparse.ArgumentParser(
        description="Preprocess a SCRAP dataset")
    argparser.add_argument("--pipeline_name", "-p",
                           help="What to name the pipeline", default=None)
    argparser.add_argument("--workspace_path", "-w",
                           help="Path to the workspace", default=None)
    argparser.add_argument("--gene_count_filter", "-g",
                           help="Throw away genes that have no cells >= this",
                           default=5, type=int)
    argparser.add_argument("--transcript_count_filter", "-t",
                           help="Throw away cells that have transcripts >= this",
                           default=None, type=int)

    args = argparser.parse_args()

    return args


def preprocess():

    args = get_arguments()

    if args.workspace_path is None:
        args.workspace_path = os.getcwd()

    workspace_path = args.workspace_path
    pipeline_name = args.pipeline_name
    gene_count_filter = args.gene_count_filter
    transcript_count_filter = args.transcript_count_filter

    gene_expression_dataset = Gene_Expression_Dataset(workspace_path)

    print("Filtering...")
    gene_expression_dataset.filter_low_transcript_counts(
        gene_count_filter)

    if transcript_count_filter is not None:
        gene_expression_dataset.filter_low_transcript_cells(
            transcript_count_filter)

    print("Normalizing cells...")
    gene_expression_dataset.normalize_cells(
        data_mode=Data_Mode.GENE_PROBABILITIES)

    print("Normalizing genes...")
    gene_expression_dataset.normalize_genes(
        Normalization_Method.LOG_PLUS_1,
        use_normalized=True,
        parameters=[5000])

    print("Transforming...")
    gene_expression_dataset.transform(
        Transformation_Method.PCA, num_dimensions=30,
        use_normalized=True)
    gene_expression_dataset.transform(
        Transformation_Method.NMF, num_dimensions=30,
        use_normalized=True)
    gene_expression_dataset.transform(
        Transformation_Method.SVD, num_dimensions=30,
        use_normalized=True)
    gene_expression_dataset.transform(
        Transformation_Method.TSNE, num_dimensions=2,
        use_normalized=True)

    gene_expression_dataset.save(pipeline_name)
