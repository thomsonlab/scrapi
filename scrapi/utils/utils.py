import numpy

import warnings

from . import fileio


def moving_average(x, w):
    return numpy.convolve(x, numpy.ones(w), "valid") / w


def convert_h5_to_sdt(
        h5_file_path,
        SDT_file_path,
        cells_as_rows=True
):

    warnings.warn(
        "convert_h5_to_sdt has been moved to utils.fileio",
        DeprecationWarning
    )

    return fileio.convert_h5_to_sdt(
        h5_file_path,
        SDT_file_path,
        cells_as_rows
    )
