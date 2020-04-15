import pandas
import random
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE

from sklearn import metrics
from sklearn.cluster import AgglomerativeClustering
import os
import math
from scipy import stats
from scipy import optimize
import csv
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.cluster import KMeans
import numpy
from sklearn import mixture
from copy import copy

from sparsedat import Sparse_Data_Table

from .Normalization_Method import Normalization_Method
from .Transformation_Method import Transformation_Method
from .Clustering_Method import Clustering_Method
from .Data_Mode import Data_Mode
from ..utils import fileio


class Gene_Expression_Dataset:

    @staticmethod
    def get_sample_name(cell_name):

        sample_start_index = cell_name.find("_")

        if sample_start_index == -1:
            return None
        else:
            return cell_name[:sample_start_index]

    @staticmethod
    def get_cell_labels_file_path(dataset_path):
        return os.path.join(dataset_path, "labels.csv")

    @staticmethod
    def get_cell_transcript_counts_file_path(dataset_path):
        return os.path.join(dataset_path, "barcode_transcript_counts.sdt")

    @staticmethod
    def initialize_dataset(dataset_path, seed_matrices_file_path):

        if isinstance(seed_matrices_file_path, str):

            seed_matrices_file_path = [seed_matrices_file_path]

        if len(seed_matrices_file_path) > 1:
            raise NotImplementedError("Haven't implemented list support for"
                                      " SDT")

        seed_matrix_file_path = seed_matrices_file_path[0]

        # If the user passed an H5 file, convert it to SDT - we can use this
        # directly as the barcode counts file
        if seed_matrix_file_path.endswith(".h5"):

            cell_transcript_counts_file_path = \
                Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
                    dataset_path
                )

            cell_transcript_counts = fileio.convert_h5_to_sdt(
                seed_matrix_file_path,
                cell_transcript_counts_file_path
            )
        elif seed_matrix_file_path.endswith(".sdt"):
            cell_transcript_counts = Sparse_Data_Table(
                seed_matrix_file_path,
                load_on_demand=False
            )

        # Filter out zero cells/genes - no point in keeping these around
        present_cells = cell_transcript_counts.sum(axis=1) > 0
        present_genes = cell_transcript_counts.sum(axis=0) > 0

        cell_transcript_counts = \
            cell_transcript_counts[present_cells, present_genes]

        cell_transcript_counts.save(
            Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
                dataset_path))

        samples = [
            Gene_Expression_Dataset.get_sample_name(cell_name)
            for cell_name in cell_transcript_counts.row_names]

        samples = set(samples)

        if len(samples) == 1:
            samples = set()

        # If an existing labels file exists, don't overwrite it
        labels_file_path = \
            Gene_Expression_Dataset.get_cell_labels_file_path(dataset_path)

        if os.path.exists(labels_file_path):
            _, label_cells = Gene_Expression_Dataset.get_label_cells_from_file(
                labels_file_path
            )
        else:
            label_cells = {}

        for sample_name in samples:
            label_cells[sample_name] = set()

        for cell_name in cell_transcript_counts.row_names:
            sample_name = Gene_Expression_Dataset.get_sample_name(cell_name)

            if sample_name is not None:
                label_cells[sample_name].add(cell_name)

        Gene_Expression_Dataset.write_label_cells_to_file(
            label_cells,
            Gene_Expression_Dataset.get_cell_labels_file_path(dataset_path))

    @staticmethod
    def get_label_cells_from_file(file_path):

        cell_labels = {}
        labels = []

        cell_labels_file = open(file_path, "r")

        cell_labels_reader = csv.reader(cell_labels_file)

        for row in cell_labels_reader:
            label = row[0]
            cells = row[1:]

            cell_labels[label] = set(cells)
            labels.append(label)

        cell_labels_file.close()

        return labels, cell_labels

    @staticmethod
    def write_label_cells_to_file(label_cells, file_path):

        cell_labels_file = open(file_path, "w", newline="")
        writer = csv.writer(cell_labels_file)
        for label, cells in label_cells.items():
            label_array = [label]
            label_array.extend(list(cells))
            writer.writerow(label_array)
        cell_labels_file.close()

    def __init__(self, dataset_path, name=None, seed_matrix_file_path=None):

        if not os.path.exists(dataset_path):
            os.makedirs(dataset_path)
        elif os.path.isfile(dataset_path):
            raise Exception("Requested to initialize a dataset in '%s' \
                but this is a file, not a folder!" % dataset_path)

        if seed_matrix_file_path is not None:
            Gene_Expression_Dataset.initialize_dataset(
                dataset_path, seed_matrix_file_path)

        self._dataset_path = dataset_path

        # A Sparse Data Table of cells and their transcript counts. Each row
        # is a cell, each column a named transcript
        self._cell_transcript_counts = None

        # A Sparse Data Table of the raw barcodes and transcript counts
        self._barcode_transcript_counts = None

        # A numpy array storing the total transcript count of each cell
        self._cell_total_transcript_counts = None

        self._zero_genes = None
        self._pca = None
        self._NMF_model = None
        self._SVD_model = None
        self._transformed = {}
        self._label_cells = {}
        self._labels = []
        self._transcript_means = None
        self._normalized_cell_transcript_counts = None

        # The mean of each transcript, normalized
        self._normalized_transcript_means = None

        self._gene_metadata = None

        if name is None:
            self._load_dataset_from_path()
        else:
            self.load(name)

        self._transcript_count_threshold = 0

    def reload(self):

        self.__init__(self._dataset_path)

    def get_labels(self):

        return list(self._label_cells.keys())

    def filter_low_transcript_counts(self, threshold):

        self._transcript_count_threshold = threshold

        if threshold <= 0:
            return

        transcripts_above_threshold = \
            (self._cell_transcript_counts.max(axis=0) >= threshold)

        self._cell_transcript_counts = \
            self._cell_transcript_counts[:, transcripts_above_threshold]

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_cell_transcript_counts = \
                self._normalized_cell_transcript_counts[
                :, transcripts_above_threshold]
            self._normalized_transcript_means = \
                self._normalized_transcript_means[transcripts_above_threshold]

        self._transcript_means = self._transcript_means[
            transcripts_above_threshold]

    def filter_genes(self, include_genes=None, exclude_genes=None):
        """
        Filter out genes that either not in the include list, or in the exclude
        list

        :param include_genes: Genes to keep
        :param exclude_genes: Genes to remove
        :return: None
        """

        gene_indices = []
        gene_list = self._cell_transcript_counts.column_names

        if include_genes is not None:
            for gene in include_genes:
                gene_indices.append(gene_list.index(gene))
        if exclude_genes is not None:
            raise NotImplementedError("Some day...")

        self._cell_transcript_counts = \
            self._cell_transcript_counts[:, gene_indices]

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_cell_transcript_counts = \
                self._normalized_cell_transcript_counts[:, gene_indices]

            self._normalized_transcript_means = \
                self._normalized_transcript_means[gene_indices]

        self._transcript_means = self._transcript_means[gene_indices]

    def filter_low_transcript_cells(self, threshold):

        cells_above_threshold = \
            self._cell_transcript_counts.sum(axis=1) >= threshold

        self._cell_transcript_counts = self._cell_transcript_counts[
                                       cells_above_threshold, :]

        self._cell_total_transcript_counts = self._cell_total_transcript_counts[
            cells_above_threshold
        ]

        # Now that we have less cells, we re-apply our transcript filter
        self.filter_low_transcript_counts(self._transcript_count_threshold)

        self._transcript_means = self._cell_transcript_counts.mean(axis=0)

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_cell_transcript_counts = \
                self._normalized_cell_transcript_counts[
                cells_above_threshold, :]
            self._normalized_transcript_means = \
                self._normalized_cell_transcript_counts.mean(axis=0)

    def filter_noise_barcodes(
            self,
            min_num_cells=1000,
            max_num_cells=15000,
            num_sources_noise_expected=1,
            min_num_cell_types_expected=4,
            max_num_clusters=15,
            log_multiplication_factor=5000
    ):
        """
        Removes estimated noise barcodes based on their transcript count
        distribution.

        :param min_num_cells:
        :param max_num_cells:
        :param num_sources_noise_expected:
        :param min_num_cell_types_expected:
        :param max_num_clusters:
        :param log_multiplication_factor
        :return:
        """

        min_num_clusters = num_sources_noise_expected + \
            min_num_cell_types_expected

        sorted_transcript_counts = \
            numpy.sort(self._cell_total_transcript_counts)

        transcript_count_thresholds = int(
            sorted_transcript_counts[-max_num_cells]), int(
            sorted_transcript_counts[-min_num_cells])
        min_transcript_count_threshold = transcript_count_thresholds[0]

        # Get a list of cells above the threshold
        cell_indices_above_threshold = \
            self._cell_total_transcript_counts > min_transcript_count_threshold

        # Get the full gene count matrix and total transcript counts associated
        # with these cells
        filtered_cell_gene_counts = \
            self._cell_transcript_counts[cell_indices_above_threshold, :]
        filtered_total_transcript_counts = \
            self._cell_total_transcript_counts[cell_indices_above_threshold]

        # Filter out zero genes to make it easier on PCA
        non_zero_genes = filtered_cell_gene_counts.sum(axis=0) > 0
        filtered_cell_gene_counts = filtered_cell_gene_counts[:, non_zero_genes]

        # Normalize the cell gene counts - first divide by the sum to get
        # transcripts per cell
        filtered_cell_gene_counts.divide(filtered_cell_gene_counts.sum(axis=1))

        # Multiply by a factor to separate single counts from zero counts
        filtered_cell_gene_counts.multiply(log_multiplication_factor)

        # Add one before taking log
        filtered_cell_gene_counts.add(1)

        # Log scale
        filtered_cell_gene_counts.log10()

        pca = PCA(n_components=50)
        transformed_PCA = pca.fit_transform(
            filtered_cell_gene_counts.to_array())

        cluster_range = range(min_num_clusters, max_num_clusters)

        max_silhouette_score = -numpy.inf
        best_num_clusters = None
        best_clusters = None

        for num_clusters in cluster_range:

            clustering = AgglomerativeClustering(n_clusters=num_clusters)

            clusters = clustering.fit_predict(transformed_PCA)

            silhouette_score = metrics.silhouette_score(
                transformed_PCA, clusters)

            if silhouette_score > max_silhouette_score:
                max_silhouette_score = silhouette_score
                best_num_clusters = num_clusters
                best_clusters = clusters

        num_clusters = best_num_clusters
        clusters = best_clusters

        # Get a range of thresholds to test from the min to the max - this is to
        # establish a slope of cell dropoff per cluster as we increase the
        # transcript count threshold
        fine_grain_transcript_count_thresholds = range(
            min_transcript_count_threshold, transcript_count_thresholds[1], 50)

        # A dataframe that lists the number of cells in each cluster above the
        # thresholds
        cluster_size_from_most_to_least = pandas.DataFrame(
            index=list(range(num_clusters)),
            columns=list(fine_grain_transcript_count_thresholds)
        )

        for threshold in fine_grain_transcript_count_thresholds:

            for cluster in range(num_clusters):
                # Get how many cells are in this cluster at this threshold
                cluster_cell_counts = \
                    filtered_total_transcript_counts[(clusters == cluster) & (
                                filtered_total_transcript_counts > threshold)]

                cluster_size_from_most_to_least.loc[cluster, threshold] = \
                    cluster_cell_counts.shape[0]

        final_cluster_loss = \
            cluster_size_from_most_to_least.iloc[:, -1] / \
            cluster_size_from_most_to_least.iloc[:, 0]
        cluster_losses = \
            cluster_size_from_most_to_least.values[:, 1:] / \
            cluster_size_from_most_to_least.values[:, 0].reshape((-1, 1))

        cluster_losses[numpy.isnan(cluster_losses)] = 0

        noise_clusterer = AgglomerativeClustering(n_clusters=2)
        noise_clusters = noise_clusterer.fit_predict(cluster_losses)

        signal_cluster = noise_clusters[numpy.where(
            final_cluster_loss.values == final_cluster_loss.values.max())[0][0]]

        # Initialize a boolean array with all false
        valid_cells = numpy.array([False] * filtered_cell_gene_counts.num_rows)

        # Loop through all clusters and mark any cells that are part of a valid
        # cluster as valid
        for cluster in range(num_clusters):

            if noise_clusters[cluster] == signal_cluster:
                valid_cells = valid_cells | (clusters == cluster)

        valid_cell_barcodes = \
            filtered_cell_gene_counts[valid_cells].row_names

        self.filter_cells(valid_cell_barcodes, exclude=False)

    def filter_cells(self, cell_barcodes, exclude=True):
        """
        Remove cells from the dataset.

        :param cell_barcodes: The list of barcodes to remove (if exclude is
            True), or keep (if exclude is False)
        :param exclude: Whether to exclude or include the given barcodes
        :return: None
        """

        if exclude:
            all_cells = set(self._cell_transcript_counts.row_names)
            cell_barcodes = all_cells - set(cell_barcodes)
            cell_barcodes = list(cell_barcodes)

        self._cell_transcript_counts = \
            self._cell_transcript_counts[cell_barcodes]

        self._cell_total_transcript_counts = self._barcode_transcript_counts[
            cell_barcodes
        ].sum(axis=1)

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_cell_transcript_counts = \
                self._normalized_cell_transcript_counts[cell_barcodes]

    def filter_unlabeled_cells(self, labels=None):
        """
        Remove any cells that don't have one of the provided labels. If no
        labels are provided, removes cells that have no labels.

        :param labels: A list of labels to filter by

        :return None
        """

        if labels is None:
            raise NotImplementedError("Must specify labels for now")

        labeled_cells = set()
        for label in labels:
            if label not in self._labels:
                continue
            labeled_cells.update(self._label_cells[label])

        labeled_cells = list(labeled_cells)

        self.filter_cells(labeled_cells, exclude=False)

    def normalize_cells(
            self, data_mode=Data_Mode.READS_PER_MILLION_TRANSCRIPTS,
            use_normalized=False
    ):

        if not use_normalized or \
                self._normalized_cell_transcript_counts is None:
            cell_transcript_counts = copy(self._cell_transcript_counts)
        else:
            cell_transcript_counts = self._normalized_cell_transcript_counts

        if data_mode == Data_Mode.READS_PER_MILLION_TRANSCRIPTS:

            cell_transcript_counts.divide(
                self._cell_total_transcript_counts)

            cell_transcript_counts.multiply(1e6)

        elif data_mode == Data_Mode.GENE_PROBABILITIES:

            cell_transcript_counts.divide(
                self._cell_total_transcript_counts)

        self._normalized_cell_transcript_counts = cell_transcript_counts

        self._normalized_transcript_means = \
            self._normalized_cell_transcript_counts.mean(axis=0)

    def label_cells(self, label, cells):

        if label not in self._label_cells:
            self._label_cells[label] = set()

        valid_cells = set(self._cell_transcript_counts.row_names).intersection(
            cells
        )

        self._label_cells[label] = self._label_cells[label].union(valid_cells)

    def delete_label(self, label):

        if label not in self._label_cells:
            return

        del self._label_cells[label]

    def rename_label(self, old_label, new_label):

        if old_label not in self._label_cells:
            raise ValueError("Can't rename label that doesn't exist!")

        self._label_cells[new_label] = self._label_cells[old_label]

        self.delete_label(old_label)

    def transform(self, method=Transformation_Method.PCA, num_dimensions=2,
                  use_normalized=False):

        if use_normalized and self._normalized_cell_transcript_counts is None:
            use_normalized = False

        if use_normalized:
            cell_transcript_counts = self._normalized_cell_transcript_counts
        else:
            cell_transcript_counts = self._cell_transcript_counts

        if method == Transformation_Method.PCA:

            self._pca = PCA(n_components=num_dimensions)

            transformed = self._pca.fit_transform(
                cell_transcript_counts.to_array())

            self._transformed[method] = pandas.DataFrame(transformed)

            self._transformed[method].columns = \
                ["PC_%i" % i for i in range(1, num_dimensions + 1)]

        elif method == Transformation_Method.TSNE:

            if Transformation_Method.PCA in self._transformed:
                transformed = TSNE(
                    verbose=True, perplexity=30, n_components=num_dimensions). \
                    fit_transform(
                    self._transformed[Transformation_Method.PCA])
            else:
                transformed = TSNE(
                    verbose=True, perplexity=30, n_components=num_dimensions). \
                    fit_transform(cell_transcript_counts.to_array())

            self._transformed[Transformation_Method.TSNE] = \
                pandas.DataFrame(transformed)

            self._transformed[Transformation_Method.TSNE].columns = \
                ["tSNE_%i" % i for i in range(1, num_dimensions + 1)]
        elif method == Transformation_Method.NMF:

            self._NMF_model = NMF(
                n_components=num_dimensions, solver="mu", init="random",
                beta_loss="kullback-leibler", max_iter=500, alpha=0.1,
                l1_ratio=0.5)

            transformed = self._NMF_model.fit_transform(
                cell_transcript_counts.to_array())

            self._transformed[method] = pandas.DataFrame(transformed)

            self._transformed[method].columns = \
                ["NMF_%i" % i for i in range(1, num_dimensions + 1)]
        elif method == Transformation_Method.SVD:

            self._SVD_model = TruncatedSVD(n_components=num_dimensions)

            transformed = self._SVD_model.fit_transform(
                cell_transcript_counts.to_array())

            self._transformed[method] = pandas.DataFrame(transformed)

            self._transformed[method].columns = \
                ["NMF_%i" % i for i in range(1, num_dimensions + 1)]

        self._transformed[method].index = cell_transcript_counts.row_names

    @property
    def num_cells(self):

        return self._cell_transcript_counts.shape[0]

    def get_cells(self, labels=None, union=False):
        """
        Get all the cell names that match the given query.

        :param labels: A list of labels
        :param union:
        :return: A set of cells
        """

        # To make dynamic inspection from Pycharm understand this is an empty
        # list...
        if labels is None:
            labels = []

        if len(labels) == 0:
            return set(self._cell_transcript_counts.row_names)

        if isinstance(labels, str):
            return self._label_cells[labels].intersection(
                set(self._cell_transcript_counts.row_names))
        elif not union:
            cells = set(self._cell_transcript_counts.row_names)

            for label in labels:
                cells = cells.intersection(self._label_cells[label])

            return cells
        else:
            cells = set()
            for label in labels:
                label_cells = set(self._cell_transcript_counts.row_names)
                label_cells = label_cells.intersection(
                    self._label_cells[label])

                cells = cells.union(label_cells)
            return cells

    def get_genes(self):

        return self._cell_transcript_counts.column_names

    def get_cell_gene_expression(self, transform=None):

        if transform is None:
            return self._cell_transcript_counts
        else:
            return self._transformed[transform]

    def get_label_cells(self):
        return self._label_cells.copy()

    def get_cell_gene_expression_by_label(self, transform=None):

        label_cells = {}

        for label in self._label_cells:

            cell_names = self.get_cells(label)

            if transform is None:
                label_cells[label] = self._cell_transcript_counts[cell_names]
            else:
                label_cells[label] = \
                    self._transformed[transform].loc[cell_names]

        return label_cells

    def get_label_means(self, random_shuffle=False, is_median=False):

        if random_shuffle:
            label_counts = self.get_label_counts()
            labels = [label_count[0] for label_count in label_counts]
            label_weights = [label_count[1] for label_count in label_counts]

            label_cells = {}

            for label, _ in label_counts:
                label_cells[label] = set()

            for cell_name in self._cell_transcript_counts.row_names:
                sample = random.choices(labels, weights=label_weights)
                label_cells[sample[0]].add(cell_name)

        else:
            label_cells = self._label_cells

        label_means = pandas.DataFrame()

        for label in label_cells:

            cell_names = self.get_cells(label)

            label_SDT = self._cell_transcript_counts[list(cell_names)]

            if is_median:
                label_mean = pandas.Series(
                    label_SDT.median(axis=0),
                    index=self._cell_transcript_counts.column_names
                )
            else:
                label_mean = pandas.Series(
                    label_SDT.mean(axis=0),
                    index=self._cell_transcript_counts.column_names
                )

            label_mean.name = label
            label_means = label_means.append(label_mean)

        return label_means

    def get_cell_gene_differential(self, gene):

        cells_gene_count = self._cell_transcript_counts[:, gene]
        cells_gene_count = cells_gene_count.to_array().astype(numpy.float)

        non_zero_min = cells_gene_count[cells_gene_count > 0].min()

        cells_gene_count[cells_gene_count == 0] = non_zero_min / 2

        gene_index = self._cell_transcript_counts.get_column_index(gene)

        cells_gene_differential = numpy.divide(
            cells_gene_count, self._transcript_means[gene_index])

        cells_gene_differential = numpy.log2(cells_gene_differential)

        return pandas.Series(
            cells_gene_differential.squeeze(),
            index=self._cell_transcript_counts.row_names
        )

    def compare_gene_expression(self, label_1, label_2=None,
                                differential_clusters=None,
                                use_normalized=True):

        if use_normalized and self._normalized_cell_transcript_counts is None:
            use_normalized = False

        label_1_cells = self.get_cells(label_1)

        if label_2 is not None:
            label_2_cells = self.get_cells(label_2)
        else:
            label_2_cells = \
                set(self._cell_transcript_counts.row_names). \
                    difference(label_1_cells)

        if differential_clusters is None or len(differential_clusters) == 0:
            gene_value_scores = self.get_gene_value_scores(
                label_1_cells, label_2_cells, use_normalized)

            gene_DE = pandas.DataFrame.from_dict(
                gene_value_scores, orient="index")

            gene_DE.columns = ["Log2 Change", "p-value", "difference",
                               "Group 1 Mean", "Group 1 SD", "Group 2 Mean",
                               "Group 2 SD"]
        else:

            gene_DE = pandas.DataFrame(
                columns=[
                    "Cluster", "Log2 Change", "p-value", "difference",
                    "Group 1 Mean", "Group 1 SD", "Group 2 Mean", "Group 2 SD"])

            for cluster in differential_clusters:
                label_1_cluster_cells = \
                    label_1_cells.intersection(self.get_cells(cluster))
                label_2_cluster_cells = \
                    label_2_cells.intersection(self.get_cells(cluster))

                if len(label_1_cluster_cells) <= 2 or \
                        len(label_2_cluster_cells) <= 2:
                    continue

                cluster_gene_value_scores = self.get_gene_value_scores(
                    label_1_cluster_cells, label_2_cluster_cells,
                    use_normalized)

                cluster_gene_DE = pandas.DataFrame.from_dict(
                    cluster_gene_value_scores, orient="index")

                cluster_gene_DE.columns = [
                    "Log2 Change", "p-value", "difference", "Group 1 Mean",
                    "Group 1 SD", "Group 2 Mean", "Group 2 SD"
                ]

                cluster_gene_DE["Cluster"] = cluster

                gene_DE = gene_DE.append(cluster_gene_DE)

        p_values = gene_DE["p-value"]

        _, p_values, _, _ = multipletests(p_values, method="fdr_bh")

        gene_DE["p-value"] = p_values

        return gene_DE

    def get_gene_value_scores(self, cells_1, cells_2, use_normalized):

        gene_value_scores = {}

        all_cells = cells_1.union(cells_2)

        if use_normalized:
            cell_gene_counts = self._normalized_cell_transcript_counts[
                               list(all_cells), :
                               ]
        else:
            cell_gene_counts = self._cell_transcript_counts[list(all_cells), :]

        non_zero_genes = cell_gene_counts.sum(axis=0) > 0
        cell_gene_counts = cell_gene_counts[:, non_zero_genes]

        min_value = cell_gene_counts.min()

        cells_1_gene_counts = cell_gene_counts[list(cells_1), :].to_array()
        cells_2_gene_counts = cell_gene_counts[list(cells_2), :].to_array()

        for gene_index, gene in enumerate(cell_gene_counts.column_names):

            sample_1_values = cells_1_gene_counts[:, gene_index].squeeze()
            sample_2_values = cells_2_gene_counts[:, gene_index].squeeze()

            sample_1_mean = sample_1_values.mean()
            sample_1_SD = sample_1_values.std()
            sample_2_mean = sample_2_values.mean()
            sample_2_SD = sample_2_values.std()

            if sample_1_mean == 0:
                sample_1_mean_for_log = min_value / 2
            else:
                sample_1_mean_for_log = sample_1_mean
            if sample_2_mean == 0:
                sample_2_mean_for_log = min_value / 2
            else:
                sample_2_mean_for_log = sample_2_mean

            log_2_fold_change = math.log2(sample_1_mean_for_log /
                                          sample_2_mean_for_log)

            difference, p_value = stats.ks_2samp(sample_1_values,
                                                 sample_2_values)

            gene_value_scores[gene] = (log_2_fold_change, p_value, difference,
                                       sample_1_mean, sample_1_SD,
                                       sample_2_mean, sample_2_SD)

        return gene_value_scores

    def get_gene_expression_for_cell(self, cell):

        cell_gene_counts = self._cell_transcript_counts[cell]
        cell_gene_de = copy(cell_gene_counts)

        gene_mins = self._cell_transcript_counts.min(axis=0)
        cell_gene_de[cell_gene_de == 0] = gene_mins[cell_gene_de == 0]

        cell_gene_de = numpy.log2(cell_gene_de / self._transcript_means)

        cell_gene_expression = numpy.array([
            cell_gene_counts,
            cell_gene_de
        ])

        cell_gene_expression = pandas.DataFrame(
            cell_gene_expression.transpose(),
            columns=["Count", "Log2 Fold Change"],
            index=self._cell_transcript_counts.column_names
        )

        return cell_gene_expression

    def get_gene_counts(self, genes, filter_labels=None, normalized=False):

        cells = self.get_cells(filter_labels)

        if normalized:
            return pandas.Series(
                self._normalized_cell_transcript_counts[list(cells), genes]
                    .to_array().squeeze()
            )
        else:
            return pandas.Series(
                self._cell_transcript_counts[list(cells), genes]
                    .to_array().squeeze()
            )

    def get_label_counts(self, filter_labels=None, union=False):

        label_counts = {}

        cells = self.get_cells(filter_labels, union=union)

        total_cells = len(cells)

        if total_cells == 0 or len(self._label_cells) == 0:
            empty = pandas.DataFrame(columns=["# Cells", "Ratio"])
            return empty

        for label in self._label_cells.keys():

            num_cells = len(cells.intersection(self.get_cells(label)))
            cell_ratio = num_cells / total_cells

            if num_cells > 0:
                label_counts[label] = (num_cells, cell_ratio)

        df = pandas.DataFrame.from_dict(label_counts, orient="index")
        df.columns = ["# Cells", "Ratio"]

        return df

    def normalize_genes(self, method=Normalization_Method.STD, by_label=False,
                        use_normalized=True, parameters=None):

        if by_label:
            raise NotImplementedError("By label abandoned with switch to SDT")

        if not use_normalized or \
                self._normalized_cell_transcript_counts is None:
            self._normalized_cell_transcript_counts = \
                copy(self._cell_transcript_counts)

        if method == Normalization_Method.STD:

            gene_sds = self._normalized_cell_transcript_counts.std(axis=0)

            self._normalized_cell_transcript_counts.divide(
                gene_sds)

        elif method == Normalization_Method.STANDARDIZATION:

            gene_means = self._normalized_cell_transcript_counts.mean(axis=0)
            gene_sds = self._normalized_cell_transcript_counts.std(axis=0)

            self._normalized_cell_transcript_counts.subtract(gene_means)
            self._normalized_cell_transcript_counts.divide(gene_sds)

        elif method == Normalization_Method.LOG_PLUS_1:

            self._normalized_cell_transcript_counts.multiply(parameters[0])
            self._normalized_cell_transcript_counts.add(1)
            self._normalized_cell_transcript_counts.log10()

        elif method == Normalization_Method.SQUARE_ROOT:

            self._normalized_cell_transcript_counts.sqrt()

        elif method == Normalization_Method.L2FC:

            raise NotImplementedError()

        elif method == Normalization_Method.ECDF:

            raise NotImplementedError()

            # gene_index = 0
            #
            # for gene, gene_counts in \
            #         self._normalized_cell_transcript_counts.iterrows():
            #
            #     value_counts = gene_counts.value_counts()
            #     eCDF = value_counts.sort_index().cumsum() * 1. / \
            #            self.num_cells
            #
            #     value_count_map = {}
            #     for i, j in eCDF.iteritems():
            #         value_count_map[i] = j
            #
            #     for cell, gene_count in gene_counts.iteritems():
            #         gene_counts[cell] = value_count_map[gene_count]
            #
            #     gene_index += 1

        self._normalized_transcript_means = \
            self._normalized_cell_transcript_counts.mean(axis=0)

    def save_labels(self):

        Gene_Expression_Dataset.write_label_cells_to_file(
            self._label_cells, self._get_cell_labels_file_path())

    def load(self, name):

        self._barcode_transcript_counts = Sparse_Data_Table(
            self._get_cell_transcript_counts_file_path()
        )

        cell_transcript_counts_path = os.path.join(
            self._dataset_path, "cell_transcript_counts_%s.sdt" % name)

        self._cell_transcript_counts = Sparse_Data_Table(
            cell_transcript_counts_path
        )

        self._cell_total_transcript_counts = \
            self._cell_transcript_counts.sum(axis=1)

        normalized_path = os.path.join(self._dataset_path,
                                       "normalized_%s.sdt" % name)

        if os.path.isfile(normalized_path):
            self._normalized_cell_transcript_counts = \
                Sparse_Data_Table(normalized_path)

            self._normalized_transcript_means = \
                self._normalized_cell_transcript_counts.mean(axis=0)

        for method_name, method in \
                Transformation_Method.__members__.items():

            file_name = "transformed_%s_%s.csv" % (method_name, name)

            file_path = os.path.join(self._dataset_path, file_name)

            if not os.path.isfile(file_path):
                continue

            self._transformed[method] = fileio.read_pandas_csv(
                file_path
            )

        self._initialize_cache()

    def save(self, name):

        self.save_labels()

        cell_transcript_counts_path = os.path.join(
            self._dataset_path, "cell_transcript_counts_%s.sdt" % name)

        self._cell_transcript_counts.save(cell_transcript_counts_path)

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_cell_transcript_counts.save(
                os.path.join(self._dataset_path, "normalized_%s.sdt" % name))

        for method_name, method in \
                Transformation_Method.__members__.items():

            if method not in self._transformed:
                continue

            file_name = "transformed_%s_%s.csv" % (method_name, name)

            file_path = os.path.join(self._dataset_path, file_name)

            fileio.write_pandas_csv(
                self._transformed[method],
                file_path)

    def get_cell_transcript_counts(
            self,
            filter_labels=None,
            normalized=False,
            genes=None
    ):

        if not filter_labels:
            if normalized:
                cell_transcript_counts = self._normalized_cell_transcript_counts
            else:
                cell_transcript_counts =  self._cell_transcript_counts
        else:
            cells = self.get_cells(filter_labels)
            if normalized:
                cell_transcript_counts = \
                    self._normalized_cell_transcript_counts[list(cells), :]
            else:
                cell_transcript_counts = \
                    self._cell_transcript_counts[list(cells), :]

        if genes is not None:
            return cell_transcript_counts[:, genes]
        else:
            return cell_transcript_counts

    def get_cell_total_transcript_counts(self):
        return self._cell_total_transcript_counts

    def auto_cluster(self, num_clusters=20,
                     transformation_method=Transformation_Method.PCA,
                     clustering_method=Clustering_Method.K_MEANS):

        if num_clusters is None:
            return

        data_transformed = self._transformed[transformation_method]

        if clustering_method == Clustering_Method.K_MEANS:
            clusterer = KMeans(n_clusters=num_clusters, random_state=0)
            fitted = clusterer.fit(data_transformed)
            clusters = fitted.labels_
        elif clustering_method == Clustering_Method.GMM:
            clusterer = mixture.GaussianMixture(n_components=num_clusters)
            fitted = clusterer.fit(data_transformed)
            clusters = fitted.predict(data_transformed)
        elif clustering_method == Clustering_Method.MAX_FEATURE:
            data_transformed = pandas.DataFrame(data_transformed)

            num_columns = len(data_transformed.columns)
            columns = data_transformed.columns[0:min(num_clusters, num_columns)]
            data_transformed = data_transformed[columns]

            clusters = data_transformed.idxmax(axis=1)

            for cluster_index, cluster in enumerate(clusters):
                clusters[cluster_index] = \
                    data_transformed.columns.get_loc(cluster)
        else:
            raise ValueError("Invalid clustering method")

        labels_to_delete = []

        for label in self._label_cells.keys():
            if label.find("Auto Cluster") != -1:
                labels_to_delete.append(label)

        for label in labels_to_delete:
            self.delete_label(label)

        for cluster_index in range(num_clusters):
            cluster_cells = data_transformed[
                clusters == cluster_index
                ]

            self.label_cells(
                "Auto Cluster %i" % cluster_index, cluster_cells.index)

    def get_matched_clusters(self, label_1, label_2=None, num_clusters=20,
                             transformation_method=Transformation_Method.PCA,
                             clustering_method=Clustering_Method.K_MEANS):

        label_1_cells = list(self.get_cells(label_1))
        label_2_cells = list(self.get_cells(label_2))

        label_1_transformed = \
            self._transformed[transformation_method].loc[label_1_cells]
        label_2_transformed = \
            self._transformed[transformation_method].loc[label_2_cells]

        if clustering_method == Clustering_Method.K_MEANS:

            label_1_k_means = KMeans(n_clusters=num_clusters, random_state=0)
            label_2_k_means = KMeans(n_clusters=num_clusters, random_state=0)

            label_1_fitted = label_1_k_means.fit(label_1_transformed)
            label_2_fitted = label_2_k_means.fit(label_2_transformed)

            label_1_cluster_centers = label_1_fitted.cluster_centers_
            label_2_cluster_centers = label_2_fitted.cluster_centers_

        elif clustering_method == Clustering_Method.GMM:
            label_1_mixture = mixture.GaussianMixture(n_components=num_clusters)
            label_2_mixture = mixture.GaussianMixture(n_components=num_clusters)

            label_1_fitted = label_1_mixture.fit(label_1_transformed)
            label_2_fitted = label_2_mixture.fit(label_2_transformed)

            label_1_cluster_centers = label_1_fitted.means_
            label_2_cluster_centers = label_2_fitted.means_
        else:
            raise NotImplementedError()

        cluster_distances = numpy.empty((num_clusters, num_clusters))

        labels_to_delete = []

        for label in self._label_cells.keys():
            if label.find("Auto Cluster") != -1:
                labels_to_delete.append(label)

        for label in labels_to_delete:
            self.delete_label(label)

        for label_1_cluster_index in range(num_clusters):
            label_1_cluster_center = \
                label_1_cluster_centers[label_1_cluster_index]
            for label_2_cluster_index in range(0, num_clusters):
                label_2_cluster_center = \
                    label_2_cluster_centers[label_2_cluster_index]
                distance = numpy.linalg.norm(
                    label_1_cluster_center - label_2_cluster_center)
                cluster_distances[label_1_cluster_index][label_2_cluster_index] \
                    = distance

        cluster_assignments = optimize.linear_sum_assignment(
            cluster_distances)

        if clustering_method == Clustering_Method.K_MEANS:
            label_1_clusters = label_1_fitted.labels_
            label_2_clusters = label_2_fitted.labels_
        elif clustering_method == Clustering_Method.GMM:
            label_1_clusters = label_1_fitted.predict(label_1_transformed)
            label_2_clusters = label_2_fitted.predict(label_2_transformed)
        else:
            raise NotImplementedError()

        for cluster_index in range(num_clusters):
            label_2_cluster_index = cluster_assignments[1][cluster_index]

            label_1_cluster_cells = label_1_transformed[
                label_1_clusters == cluster_index
                ]

            label_2_cluster_cells = label_2_transformed[
                label_2_clusters == label_2_cluster_index
                ]

            cluster_cells = list(label_1_cluster_cells.index)
            cluster_cells.extend(list(label_2_cluster_cells.index))

            self.label_cells("Auto Cluster %i" % cluster_index,
                             cluster_cells)

    def _load_dataset_from_path(self):

        self._barcode_transcript_counts = Sparse_Data_Table(
            self._get_cell_transcript_counts_file_path()
        )

        self._cell_transcript_counts = copy(self._barcode_transcript_counts)

        self._cell_total_transcript_counts = \
            self._cell_transcript_counts.sum(axis=1)

        self._initialize_cache()

    def _initialize_cache(self):

        self._transcript_means = self._cell_transcript_counts.mean(axis=0)

        if self._normalized_cell_transcript_counts is not None:
            self._normalized_transcript_means = \
                self._normalized_cell_transcript_counts.mean(axis=0)

        self._labels, self._label_cells = \
            Gene_Expression_Dataset.get_label_cells_from_file(
                self.get_cell_labels_file_path(self._dataset_path))

    def _get_cell_labels_file_path(self):
        return Gene_Expression_Dataset.get_cell_labels_file_path(
            self._dataset_path)

    def _get_cell_transcript_counts_file_path(self):
        return Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
            self._dataset_path)
