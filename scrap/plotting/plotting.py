import numpy

from pepars.plotting import plotting

from .. import utils


def plot_barcode_transcript_counts(
        transcript_counts,
        max_points=1000,
        elbow=True,
        **kwargs):
    """
    Plots a scatter plot of number of barcodes at or above different UMI count
    thresholds. Also known as the "knee" or "elbow" plot.

    :param transcript_counts: An iterable of transcript counts
    :param max_points: The maximum number of points to display
    :param elbow: Whether this is an elbow (True) or knee (False) plot
    :return: A plotly figure object
    """

    step_size = int(numpy.floor(len(transcript_counts)/max_points))

    transcript_counts = numpy.sort(transcript_counts)[::-1]
    # transcript_counts = numpy.log10(transcript_counts)
    transcript_counts = utils.moving_average(transcript_counts, step_size)

    cell_number = [x + 1 for x in list(range(len(transcript_counts)))]

    transcript_counts = transcript_counts[::step_size]
    cell_number = cell_number[::step_size]

    if elbow:
        return plotting.plot_scatter(
            cell_number,
            transcript_counts,
            x_axis_title="# Cells",
            y_axis_title="Transcript Count",
            x_axis_log_scale=True,
            y_axis_log_scale=True,
            **kwargs
        )
    else:
        return plotting.plot_scatter(
            transcript_counts,
            cell_number,
            x_axis_title="Transcript Count",
            y_axis_title="# Cells",
            x_axis_log_scale=True,
            y_axis_log_scale=True,
            **kwargs
        )
