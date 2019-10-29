import numpy
import h5py

from sparsedat import Sparse_Data_Table


def moving_average(x, w):
    return numpy.convolve(x, numpy.ones(w), "valid") / w


def convert_h5_to_sdt(
        h5_file_path,
        SDT_file_path,
        cells_as_rows=True
):
    h5_file = h5py.File(h5_file_path)

    sdt = Sparse_Data_Table()

    data = h5_file["matrix"]["data"][()]
    indices = h5_file["matrix"]["indices"][()]
    indptr = h5_file["matrix"]["indptr"][()]

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
        h5_file["matrix"]["shape"][0],
        h5_file["matrix"]["shape"][1]
    )

    sdt.column_names = [x.decode("utf-8") for x in
                        h5_file["matrix"]["barcodes"][()]]
    sdt.row_names = [x.decode("utf-8") for x in
                     h5_file["matrix"]["features"]["name"][()]]

    # Transpose so that rows are cells
    if cells_as_rows:
        sdt.transpose()

    sdt.save(SDT_file_path)
