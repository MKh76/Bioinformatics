def space_sep(inlist):
    """
    :param inlist: an iterable separated not by space
    :return: a space-separated string of the same items
    """
    str_list = [str(el) for el in inlist]
    return " ".join(str_list)


def random_kmer_count(t, le, k):
    """
    :param t: int: number of sequence
    :param le: int: length of each sequence
    :param k: int: length of the desired k-mer
    :return: float: how many times k-mer was found in
             t number of sequences with the length l
    """

    return t * (le - k + 1) * (1/4)**k
