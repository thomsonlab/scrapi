from enum import Enum


class Normalization_Method(Enum):
    STD = 0
    ECDF = 1
    STANDARDIZATION = 2
    L2FC = 3
    LOG_PLUS_1 = 4
    SQUARE_ROOT = 5
