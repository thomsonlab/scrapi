import pandas

from pepars.utils import FASTQ as FASTQ_utils
from pepars.utils import DNA as DNA_utils
from pepars.utils import Sequence_Trie
from pepars.analysis import DNA as DNA_analysis


def get_cell_tag_counts(
        FASTQ_file_set,
        valid_cell_barcode_trie=None,
        valid_tags=None,
        index_sequence=None,
        maximum_tag_distance=0,
        maximum_index_distance=0,
        maximum_cell_barcode_distance=0,
        quality_threshold=None,
        collapse_UMIs=True,
        error_correct=False,
        cell_barcode_range=(0, 16),
        cell_barcode_read_number=1,
        tag_read_number=2,
        tag_range=(0, 8),
        UMI_range=(16, 28),
        UMI_read_number=1,
        output_frequency=1e5
):

    tag_start_index = tag_range[0]
    tag_end_index = tag_range[1]
    tag_length = tag_end_index - tag_start_index

    UMI_start_index = UMI_range[0]
    UMI_end_index = UMI_range[1]
    UMI_length = UMI_end_index - UMI_start_index

    cell_barcode_start_index = cell_barcode_range[0]
    cell_barcode_end_index = cell_barcode_range[1]
    cell_barcode_length = cell_barcode_end_index - cell_barcode_start_index

    if valid_tags:
        valid_tag_trie = Sequence_Trie(DNA_utils.NUCLEOTIDE_INDEX_MAP,
                                       allow_invalid=True)

        for tag in valid_tags:
            valid_tag_trie.add(tag)
    else:
        valid_tag_trie = None

    cell_UMI_tag_counts = {}
    num_reads_passing_quality_threshold = 0
    num_reads_matching_index = 0
    num_reads_matching_tag = 0
    num_tag_reads_matching_cell = 0
    num_reads = 0

    if quality_threshold:
        FASTQ_file_set_iterator = FASTQ_file_set.get_sequence_quality_iterator()
    else:
        FASTQ_file_set_iterator = FASTQ_file_set.get_sequence_iterator()

    for read in FASTQ_file_set_iterator:

        num_reads += 1

        if output_frequency and num_reads % output_frequency == 0:
            print("Read %i reads" % num_reads)

        if quality_threshold:

            passes_quality_score = True

            index_quality_sequence = read[1][0]

            index_quality_scores = FASTQ_utils.convert_quality_string_to_quality_score(
                index_quality_sequence)

            for quality_score in index_quality_scores:
                if quality_score < quality_threshold:
                    passes_quality_score = False
                    break

            if not passes_quality_score:
                continue

            cell_barcode_quality_sequence = read[1][cell_barcode_read_number][
                                            cell_barcode_start_index:cell_barcode_end_index]

            cell_barcode_quality_scores = FASTQ_utils.convert_quality_string_to_quality_score(
                cell_barcode_quality_sequence)

            for quality_score in cell_barcode_quality_scores:
                if quality_score < quality_threshold:
                    passes_quality_score = False
                    break

            if not passes_quality_score:
                continue

            UMI_quality_sequence = read[1][UMI_read_number][
                                   UMI_start_index:UMI_end_index]

            UMI_quality_scores = FASTQ_utils.convert_quality_string_to_quality_score(
                UMI_quality_sequence)

            for quality_score in UMI_quality_scores:
                if quality_score < quality_threshold:
                    passes_quality_score = False
                    break

            if not passes_quality_score:
                continue

            tag_quality_sequence = read[1][tag_read_number][
                                   tag_start_index:tag_end_index]

            tag_quality_scores = FASTQ_utils.convert_quality_string_to_quality_score(
                tag_quality_sequence)

            for quality_score in tag_quality_scores:
                if quality_score < quality_threshold:
                    passes_quality_score = False
                    break

            if not passes_quality_score:
                continue

            sequences = read[0]
        else:
            sequences = read

        num_reads_passing_quality_threshold += 1

        if index_sequence:
            if not DNA_utils.is_template_match(
                    index_sequence,
                    sequences[0],
                    allowable_distance=maximum_index_distance):
                continue

        num_reads_matching_index += 1

        tag_sequence = sequences[tag_read_number][tag_start_index:tag_end_index]

        if error_correct:
            cell_barcode = sequences[cell_barcode_read_number][
                           cell_barcode_start_index:cell_barcode_end_index]
            UMI = sequences[UMI_read_number][UMI_start_index:UMI_end_index]
        else:
            if valid_tag_trie:
                if not valid_tag_trie.find(tag_sequence):
                    if maximum_tag_distance == 0:
                        continue
                    else:
                        possible_tags = DNA_analysis.find_all_sequences_of_distance_n(
                            tag_sequence,
                            valid_tag_trie,
                            n=maximum_tag_distance)

                        if len(possible_tags) != 1:
                            continue

                        tag_sequence = possible_tags[0]

            num_reads_matching_tag += 1

            cell_barcode = sequences[cell_barcode_read_number][
                           cell_barcode_start_index:cell_barcode_end_index]

            if valid_cell_barcode_trie:
                if not valid_cell_barcode_trie.find(cell_barcode):
                    if maximum_cell_barcode_distance == 0:
                        continue
                    else:
                        possible_cell_barcodes = DNA_analysis.find_all_sequences_of_distance_n(
                            cell_barcode,
                            valid_cell_barcode_trie,
                            n=maximum_cell_barcode_distance)

                        if len(possible_cell_barcodes) != 1:
                            continue

                        cell_barcode = possible_cell_barcodes[0]

            num_tag_reads_matching_cell += 1

            UMI = sequences[UMI_read_number][UMI_start_index:UMI_end_index]

        cell_UMI_tag = cell_barcode + UMI + tag_sequence

        if cell_UMI_tag not in cell_UMI_tag_counts:
            cell_UMI_tag_counts[cell_UMI_tag] = 1
        else:
            cell_UMI_tag_counts[cell_UMI_tag] += 1

    cell_UMI_tag_counts_list = [(cell_UMI_tag, count) for cell_UMI_tag, count in
                                cell_UMI_tag_counts.items()]

    if error_correct:

        cell_UMI_tag_counts_list = DNA_analysis.collapse_similar_sequences(
            cell_UMI_tag_counts_list)

        filtered_cell_UMI_tag_counts = {}

        for cell_UMI_tag, count in cell_UMI_tag_counts_list:

            tag_sequence = cell_UMI_tag[cell_barcode_length + UMI_length:]

            if valid_tag_trie:
                if not valid_tag_trie.find(tag_sequence):
                    if maximum_tag_distance == 0:
                        continue
                    else:
                        possible_tags = DNA_analysis.find_all_sequences_of_distance_n(
                            tag_sequence,
                            valid_tag_trie,
                            n=maximum_tag_distance)

                        if len(possible_tags) != 1:
                            continue

                        tag_sequence = possible_tags[0]

            num_reads_matching_tag += count

            cell_barcode = cell_UMI_tag[0:cell_barcode_length]

            if valid_cell_barcode_trie:
                if not valid_cell_barcode_trie.find(cell_barcode):
                    if maximum_cell_barcode_distance == 0:
                        continue
                    else:
                        possible_cell_barcodes = DNA_analysis.find_all_sequences_of_distance_n(
                            cell_barcode,
                            valid_cell_barcode_trie,
                            n=maximum_cell_barcode_distance)

                        if len(possible_cell_barcodes) != 1:
                            continue

                        cell_barcode = possible_cell_barcodes[0]

            num_tag_reads_matching_cell += count

            UMI = cell_UMI_tag[
                  cell_barcode_length:cell_barcode_length + UMI_length]

            cell_UMI_tag = cell_barcode + UMI + tag_sequence

            if cell_UMI_tag not in filtered_cell_UMI_tag_counts:
                filtered_cell_UMI_tag_counts[cell_UMI_tag] = 1
            else:
                filtered_cell_UMI_tag_counts[cell_UMI_tag] += 1

        cell_UMI_tag_counts_list = [(cell_UMI_tag, count) for
                                    cell_UMI_tag, count in
                                    filtered_cell_UMI_tag_counts.items()]

    # Establish our tag list, either from the given list, or as deduced from the data
    if not valid_tags:
        valid_tags = set({})

        for cell_UMI_tag in cell_UMI_tag_counts_list:
            tag = cell_UMI_tag[0][cell_barcode_length + UMI_length:]
            valid_tags.add(tag)

    cell_tag_counts = {}

    tag_indices = {}

    tag_index = 0
    sorted_tag_list = []

    for tag in sorted(valid_tags):
        tag_indices[tag] = tag_index
        tag_index += 1
        sorted_tag_list.append(tag)

    for sequence, count in cell_UMI_tag_counts_list:

        cell_barcode = sequence[0:cell_barcode_length]
        tag = sequence[cell_barcode_length + UMI_length:]

        try:
            tag_index = tag_indices[tag]

            if cell_barcode not in cell_tag_counts:
                cell_tag_counts[cell_barcode] = [0 for tag in valid_tags]

            if collapse_UMIs:
                cell_tag_counts[cell_barcode][tag_index] += 1
            else:
                cell_tag_counts[cell_barcode][tag_index] += count
        except ValueError:
            pass

    cell_tag_counts_df = pandas.DataFrame.from_dict(cell_tag_counts,
                                                    orient="index")
    cell_tag_counts_df.columns = sorted_tag_list

    print("Num reads: %i" % num_reads)
    print("Reads passing quality threshold: %i/%i (%.2f%%)" % (
        num_reads_passing_quality_threshold, num_reads,
        num_reads_passing_quality_threshold / num_reads * 100))
    print("Reads matching index: %i/%i (%.2f%%)" % (
        num_reads_matching_index, num_reads_passing_quality_threshold,
        num_reads_matching_index / num_reads_passing_quality_threshold * 100))
    print("Reads matching tag: %i/%i (%.2f%%)" % (
        num_reads_matching_tag, num_reads_matching_index,
        num_reads_matching_tag / num_reads_matching_index * 100))
    print("Tag reads matching cell: %i/%i (%.2f%%)" % (
        num_tag_reads_matching_cell, num_reads_matching_tag,
        num_tag_reads_matching_cell / num_reads_matching_tag * 100))

    return cell_tag_counts_df
