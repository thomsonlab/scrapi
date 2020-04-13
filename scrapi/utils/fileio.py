import pandas
import h5py
from sparsedat import wrappers as SDT_wrappers
from sparsedat import Sparse_Data_Table


def write_pandas_csv(data_frame, file_path):
    pandas.DataFrame(data_frame)\
        .to_csv(file_path, sep=',', encoding='utf-8', chunksize=1000)


def read_pandas_csv(file_path):
    return pandas.read_csv(file_path, sep=",", header=0, index_col=0)


def load_mtx(mtx_file_path, features_file_path, barcodes_file_path):

    if features_file_path is not None:
        with open(features_file_path, "r") as row_names_file:
            row_names = [line[:-1].strip() for line in row_names_file]
            row_names = [" ".join(row_name.split("\t")[0:2])
                         for row_name in row_names]
    else:
        row_names = None

    if barcodes_file_path is not None:
        with open(barcodes_file_path, "r") as column_names_file:
            column_names = [line[:-1].strip() for line in column_names_file]

    sdt = SDT_wrappers.load_mtx(mtx_file_path)
    sdt.row_names = row_names
    sdt.column_names = column_names

    return sdt


def convert_h5_to_sdt(
        h5_file_path,
        SDT_file_path,
        cells_as_rows=True
):
    h5_file = h5py.File(h5_file_path)

    cellranger_version = 2

    if "matrix" in h5_file:
        cellranger_version = 3

    matrix_name = None

    if cellranger_version == 2:
        for key, value in h5_file.items():
            matrix_name = key
            break
    else:
        matrix_name = "matrix"

    sdt = Sparse_Data_Table()

    data = h5_file[matrix_name]["data"][()]
    indices = h5_file[matrix_name]["indices"][()]
    indptr = h5_file[matrix_name]["indptr"][()]

    if indices[0] > indices[1]:
        for column_index in range(len(indptr) - 1):
            data[indptr[column_index]:indptr[column_index + 1]] = \
                data[indptr[column_index]:indptr[column_index + 1]][::-1]
            indices[indptr[column_index]:indptr[column_index + 1]] = \
                indices[indptr[column_index]:indptr[column_index + 1]][::-1]
    else:
        for column_index in range(len(indptr) - 1):
            data[indptr[column_index]:indptr[column_index + 1]] = \
                data[indptr[column_index]:indptr[column_index + 1]]
            indices[indptr[column_index]:indptr[column_index + 1]] = \
                indices[indptr[column_index]:indptr[column_index + 1]]

    sdt.from_sparse_column_entries(
        (
            data,
            indices,
            indptr
        ),
        h5_file[matrix_name]["shape"][0],
        h5_file[matrix_name]["shape"][1]
    )

    if cellranger_version == 2:
        gene_names = [x.decode("UTF-8")
                      for x in list(h5_file[matrix_name]["gene_names"])]
        gene_ids = [x.decode("UTF-8")
                    for x in list(h5_file[matrix_name]["gene_ids"])]
    else:
        gene_names = [x.decode("UTF-8")
                      for x in list(h5_file[matrix_name]["features"]["name"])]
        gene_ids = [x.decode("UTF-8")
                    for x in list(h5_file[matrix_name]["features"]["id"])]

    gene_name_indices = {}

    disambiguated_gene_names = []

    for gene_index, gene in enumerate(gene_names):
        if gene not in gene_name_indices:
            gene_name_indices[gene] = [gene_ids[gene_index]]
        else:
            gene_name_indices[gene].append(gene_ids[gene_index])

    for gene_index, gene in enumerate(gene_names):

        if len(gene_name_indices[gene]) > 1:
            # Figure out which gene this is, as sorted by id
            this_gene_id = gene_ids[gene_index]
            duplicate_gene_index = sorted(gene_name_indices[gene]).index(
                this_gene_id)
            disambiguated_gene_name = "%s_%i" % (gene, duplicate_gene_index + 1)
            disambiguated_gene_names.append(disambiguated_gene_name)
        else:
            disambiguated_gene_names.append(gene)

    sdt.column_names = [x.decode("utf-8") for x in
                        h5_file[matrix_name]["barcodes"][()]]
    sdt.row_names = disambiguated_gene_names

    # Transpose so that rows are cells
    if cells_as_rows:
        sdt.transpose()

    sdt.save(SDT_file_path)

    return sdt
