import os

from pepars.fileio import fileio
from pepars.utils import DNA as DNA_utils
from pepars.utils import Sequence_Trie


def get_chromium_barcodes(chromium_version=3):

    URL_prefix = "https://caltech.box.com/shared/static/"

    if chromium_version == 2:
        file_name = "cell-barcodes.10x.737K-august-2016.txt"
        remote_file_URL = URL_prefix + "du9v3pfdammy3u1iwgal8fin495u6kw9.txt"
    elif chromium_version == 3:
        file_name = "3M-february-2018.txt"
        remote_file_URL = URL_prefix + "evjhzh50rk6zk2jtjpfg3myddxxmnqn5.txt"
    else:
        raise NotImplementedError("Only support versions 2 and 3 of Chromium")

    valid_barcodes_file_path = os.path.join("data", file_name)

    fileio.download_remote_file(
        file_URL=remote_file_URL,
        file_path=valid_barcodes_file_path
    )

    valid_barcode_trie = Sequence_Trie(
        DNA_utils.NUCLEOTIDE_INDEX_MAP,
        allow_invalid=True)

    with open(valid_barcodes_file_path) as valid_barcodes_file:
        for line in valid_barcodes_file.readlines():
            valid_barcode_trie.add(line.strip())

    return valid_barcode_trie
