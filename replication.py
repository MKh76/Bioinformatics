from random import randint, uniform
import tools_collection as tools
import itertools
import time
from collections import defaultdict

# ====================================================================72
# ============================================================================80
# ===================================================================================================================120


def pattern_count(text, pattern):
    """takes nucleic acid sequence as :param text and a shorter sequence as :param pattern
        :returns int: count of the pattern within the sequence
    """
    count = 0
    ran = range(len(text)-len(pattern)+1)
    for i in ran:
        if text[i:i+len(pattern)] == pattern:
            count += 1

    return count


def frequency_map(text, k):
    """takes a nucleic acid sequence as :param text
        and a length for a pattern sequence as :param k:int
        searches how many time each pattern repeats within the sequence
        :returns a dictionary {pattern: frequency} """
    freq = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        # add pattern to dict - substring: count
        if pattern in freq.keys():  # if key in dict
            new_value = freq.get(pattern) + 1
            freq[pattern] = new_value

        else:  # if key not in dict
            freq[pattern] = 1

    return freq


def frequent_sequence(text, k, m=2):
    """takes a nucleic acid sequence as :param text, a length for a pattern
        as :param k, the minimum frequency required as :param m
        :returns a list of patterns that repeat :param m times in the sequence
        also, the pattern(s)? with most frequency and its frequency"""

    freq_map_dict = frequency_map(text, k)

    two_plus_list = [key for key in freq_map_dict if freq_map_dict.get(key) >= m]
    find_max = max(list(freq_map_dict.values()))

    # retrieving key using its value
    most_freq_list = [key for key in freq_map_dict if freq_map_dict.get(key) == find_max]

    return two_plus_list, [most_freq_list, find_max]


def pattern_matching(genome, pattern):
    """takes a nucleic acid sequence as :param genome and a shorter sequence as :param pattern
        searches the genome and saves the indices of positions where pattern occurs
        :returns the list of indices as index_list and the pattern's frequency"""

    ran = range(len(genome)-len(pattern)+1)
    index_list = [i for i in ran if genome[i:i+len(pattern)] == pattern]

    return index_list, len(index_list)


def reverse_complement(pattern):
    """takes a sequence from one strand of DNA and returns the other DNA strand"""
    rev_pat_str = pattern[::-1]

    alt_rev_pat_str = ""
    for nucleotide in rev_pat_str:

        if nucleotide == "A":
            alt_rev_pat_str += "T"
        elif nucleotide == "T":
            alt_rev_pat_str += "A"
        elif nucleotide == "G":
            alt_rev_pat_str += "C"
        elif nucleotide == "C":
            alt_rev_pat_str += "G"
        else:
            print("A character other than A, T, C, G is in the input string")
            raise TypeError

    return alt_rev_pat_str


def clump_find(genome, k, le, t):
    """
    :param genome: a string of genome
    :param k: length of the desired k-mer
    :param le: length of sequence within genome to be examined at a time
    :param t: the number of repetition
    :return: a list containing the k-mer(s) that have been repeated t times,
             within the window of le, in the genome
    note: not efficient and practical for big data (use one of the algorithms)
    """
    final_list = []
    for i in range(0, (len(genome) - le)+1):
        text_window = genome[i:le+i]
        indices, length = frequent_sequence(text_window, k, t)
        final_list.extend(indices)

    return list(set(final_list))


def symbol_array(genome, symbol):
    """takes a nucleic acid sequence as :param genome
        and a single character as :param symbol
        adds the first half of the genome sequence to the end of it
        searches half of the genome at each cycle and
        :returns the symbol's frequency"""
    dic = {}
    len_gen = len(genome)
    ext_gen = genome + genome[0: len_gen//2]
    for i in range(len_gen):
        dic[i] = pattern_count(ext_gen[i: i + (len_gen//2)], symbol)

    return dic


def faster_symbol_array(genome, symbol):  # far more efficient than symbol_array()
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n // 2]

    # look at the first half of Genome to compute first array value
    array[0] = pattern_count(symbol, genome[0:n // 2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if extended_genome[i-1] == symbol:
            array[i] = array[i]-1
        if extended_genome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1

    return array


def max_min(dic, max_or_min="max"):
    """it takes a dictionary as :param dic
    :returns the maximum or minimum of its values depending the :param max_or_min"""

    if max_or_min == "max":
        peak = max(dic.values())
        keys_list = [key for key in dic if dic[key] == peak]
        return keys_list, peak

    elif max_or_min == "min":
        valley = min(dic.values())
        keys_list = [key for key in dic if dic[key] == valley]
        return keys_list, valley


def skew_array(sequence):
    """takes a nucleic acid sequence as :param sequence and calculates
        the skew which C and G nucleotides cause.
        this is based on the variance of C to G ratio
        in different half strands of bacteria genome
        :returns  dictionary {index, skew degree} """
    skew_arr = {0: 0}
    ran = range(0, len(sequence))
    for i in ran:
        if sequence[i] == 'A' or sequence[i] == 'T':
            skew_arr[i+1] = skew_arr[i]

        elif sequence[i] == 'G':
            skew_arr[i+1] = skew_arr[i] + 1

        elif sequence[i] == 'C':
            skew_arr[i+1] = skew_arr[i] - 1

    return skew_arr


def min_skew(seq):
    """takes a nucleic acid sequence as :param seq
        runs seq through skew_array() and calculates the minimum
        returns a list of positions of all the low point(s)"""
    skew_dict = skew_array(seq)
    min_val = min(skew_dict.values())
    min_skew_list = [key for key in skew_dict if skew_dict.get(key) == min_val]

    return min_skew_list


def hamming_distance(seq1, seq2):
    """takes two nucleic acid sequences as :param seq1 and :param seq2
        :returns the number of differences between them."""
    ran = range(0, len(seq1))
    ham_dis = 0

    for i in ran:
        if seq1[i] != seq2[i]:
            ham_dis += 1

    return ham_dis


def hamming_distance_thorough(seq1, seq2):
    """a more complete version of the hamming_distance() that can handle unequal sequences
        and also return by how much a sequence is longer"""
    seq = None
    redundant_symbols = 0
    to_output = ""

    if len(seq1) == len(seq2):
        seq = seq1
    elif len(seq1) < len(seq2):
        seq = seq1
        redundant_symbols += (len(seq2) - len(seq1))
        to_output = "seq2"
    elif len(seq1) > len(seq2):
        seq = seq2
        redundant_symbols += (len(seq1) - len(seq2))
        to_output = "seq1"

    ran = range(0, len(seq))
    ham_dis = 0

    for i in ran:
        if seq1[i] != seq2[i]:
            ham_dis += 1

    return ham_dis, [to_output, redundant_symbols]


def approximate_pattern_matching(text, patt, d=1):
    """similar to pattern matching with a specified tolerance for mismatches
        the tolerance can be specified by :param d
        :returns a list of matches"""
    ran = range(0, len(text) - len(patt) + 1)

    index_list = [i for i in ran if hamming_distance(patt, text[i: i + len(patt)]) <= d]
    return index_list


# bioinformatics 1
def approximate_pattern_matching_count(text, patt, d):
    """
    :param text: a string of nucleotides
    :param patt: a pattern to search for
    :param d: number of allowed mismatches
    :return: the highest repeated patterns in the text and their count
    """
    patt_count_dict = {}
    patt_len = len(patt)

    found_seqs_indices = approximate_pattern_matching(text, patt, d)
    found_seqs_string = [text[el:patt_len] for el in found_seqs_indices]
    for i in found_seqs_string:
        try:
            patt_count_dict[i]
        except KeyError:
            patt_count_dict[i] = 1
        else:
            patt_count_dict[i] += 1

    keys, peak = max_min(patt_count_dict, max_or_min="max")

    return keys, peak


def motif_count(motifs, pseudocount=False):
    """takes a list of strings. creates a key value pair in the dictionary
       where the key is each of the 'A', 'C' 'G' 'T' and the value is an empty-
       list in the first for loop, stores int=0 in the lists in the dictionary
       as many time as there are nucleotide in each motif-string
       in the 2nd for loop, iterates over each motif string and within that
       over each nucleotide and
       increases the stored 0 in the lists in the dictionary for that symbol
       in the corresponding position
        returns a dictionary, with count matrix of motifs """
    count_dict = {}
    k = len(motifs[0])
    t = len(motifs)
    initial_placement = 0
    if pseudocount:
        initial_placement = 1
    for symbol in "ACGT":
        count_dict[symbol] = []
        for j in range(k):
            count_dict[symbol].append(initial_placement)

    for i in range(t):
        for j in range(k):
            try:
                symbol = motifs[i][j]
                count_dict[symbol][j] += 1
            except IndexError:
                return count_dict

    return count_dict


def motif_profile(motif_list, pseudocount=False):
    """takes a list of motif strings as input and uses the motif_count function
       to calculate the count(motifs). then divides each number within each list
       in the dictionary by t (the number of rows) to calculate their
       proportion in that column which is the profile(motifs)"""
    profile_motifs = motif_count(motif_list, pseudocount=pseudocount)

    if not pseudocount:
        for key, value in profile_motifs.items():
            value[:] = [x / len(motif_list) for x in value]

    elif pseudocount:
        for i in range(len(motif_list[0])):
            for symbol in "ACGT":
                profile_motifs[symbol][i] /= len(motif_list) + 4

    return profile_motifs


def motif_consensus(motifs):
    """takes a list of motif strings as :param motifs, runs them through
       motif_count() iterates over each position in the value list of each
       nucleotides selects the highest and concatenates them into a string
       :return consensus"""
    k = len(motifs[0])
    count = motif_count(motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol

    return consensus


def motif_score(motifs):
    """takes a list of motif strings, calculates and :returns the total number
       as an int"""
    t = len(motifs)
    k = len(motifs[0])
    consensus = motif_consensus(motifs)
    score = 0
    for i in range(t):
        for j in range(k):
            try:
                if motifs[i][j] != consensus[j]:
                    score += 1
            except IndexError:
                return score
    return score


def motif_score2(motifs):
    """identical to motif_score but uses the hamming_distance function"""
    score = 0
    for i in range(len(motifs)):
        score += hamming_distance(motifs[i], motif_consensus(motifs))

    return score


def motif_probability_profile(k_mer, profile):
    """takes a string as :param k_mer and a dictionary of profile(motifs),
    which contains the probability of each nucleotide at specific positions
    as :param profile and :return the result of calculating the probability
     of that string using the profile dictionary  """
    k = len(k_mer)
    probability = 1

    for j in range(k):
        value = profile[k_mer[j]][j]
        probability *= value

    return probability


def most_probable_kmer(sequence, k, profile):
    """takes a sequence of genome as :param sequence,
    the length of the k_mer as :param k_int
    the dictionary containing four lists of profile(motifs) values as :param profile
    :returns the most probable motif (which has the highest profile score)"""
    ran = range(len(sequence) - k + 1)
    probability_dict = {}
    probability_list = []
    for i in ran:
        k_mer = sequence[i:i + k]

        probability_dict[k_mer] = motif_probability_profile(k_mer, profile)

    highest_value = max(probability_dict.values())
    for key, value in probability_dict.items():
        if value == highest_value:
            probability_list.append(key)

    return probability_list[0]


def greedy_motif_search(dna, k, t, pseudocount=False):
    """
    :param dna: a list of nucleotide strings as
    :param k: the length of the desired k_mers as
    :param t: the number of strings within the dna list as
    :param pseudocount: self-explanatory
    :return: a list containing t * strings which shows the best set of motifs
    """
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])

    n = len(dna[0])
    for i in range(n-k+1):
        motifs = [dna[0][i:i + k]]
        for j in range(1, t):
            p = motif_profile(motifs[0:j], pseudocount=pseudocount)
            motifs.append(most_probable_kmer(dna[j], k, p))

        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs

    return best_motifs


def motif_pseudocount(motifs):
    """same as the motif_count with only one difference
    the initial number is 1 not 0, so we incorporated this function
    into the motif_count() function using a boolean named parameter."""
    count_dict = {}
    k = len(motifs[0])
    t = len(motifs)

    for symbol in "ACGT":
        count_dict[symbol] = []
        for j in range(k):
            count_dict[symbol].append(1)

    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count_dict[symbol][j] += 1

    return count_dict


def motifs_from_profile(profile: dict, dna, k):
    """exactly like the most_probable_kmer_using_profile()
        this just repeats multiple times,
        takes a list of strings as :param dna instead of a string
        :return a list of strings of the best k-mer from each inputted dna string
        """
    dt = len(dna)
    most_probable = [most_probable_kmer(dna[i], k, profile)
                     for i in range(0, dt)]
    return most_probable


def random_motifs(dna, k, t):
    """
    :param dna: takes a list of dna strings
    :param k: length of the desired motifs
    :param t for the number of strings inside dna list (for some reason)
    returns a list of random motif strings from each of dna strings"""
    dk = len(dna[0])
    k_mer_list = []
    try:
        for i in range(0, t):
            rand = randint(0, dk-k)
            k_mer_list.append(dna[i][rand:rand+k])
    except IndexError:
        print(":param t > len(dna) -> only printed result for existing dna strings")

    return k_mer_list


def randomized_motif_search(dna, k, t, n):
    """
    :param dna: a list of dna strings
    :param k: the length of desired motifs
    :param t: the number of dna strings within the dna list (for some reason)
    :param n: the number of runs
    :return: a list of best motifs found from each string of dna
    """
    m = random_motifs(dna, k, t)
    best_motifs = m

    for i in range(n+1):
        while True:
            profile = motif_profile(dna, pseudocount=True)
            m = motifs_from_profile(profile, dna, k)
            if motif_score(m) < motif_score(best_motifs):
                best_motifs = m
            else:
                return best_motifs


def normalize(probabilities: dict):
    """
    :param probabilities: dict - key: value = k-mers: probability
    divides each value by the sum of all values in that probability.value
    :return: the dictionary which the sum of probabilities amount to 1
    """
    value_sum = 0
    normalized_dict = {}
    for value in probabilities.values():
        value_sum += value

    for key, value in probabilities.items():
        normalized_dict[key] = value / value_sum

    return normalized_dict


def weighted_die(probabilities):
    rand = uniform(0, 1)
    for key in probabilities:
        rand -= probabilities[key]
        if rand <= 0:
            return key


def profile_generated_string(text, profile, k):
    """
    :param text: a dna sequence
    :param profile: a profile dictionary
    :param k: to specify k in k-mer
    :return: a randomly generated k-mer of the text using the profile
    """
    n = len(text)
    probabilities = {}

    for i in range(n-k+1):
        probabilities[text[i:i+k]] = \
            motif_probability_profile(text[i:i + k], profile)

    probabilities = normalize(probabilities)
    return weighted_die(probabilities)


def gibbs_sampler(dna, k, t, n):
    """
    :param dna: a list of dna strings
    :param k: length of the desired motif
    :param t: the number of strings in dna list to be considered (useless)
    :param n: the number of iteration
    :return:
    """
    motifs_list = []
    pre_profile_list = []
    profile = {}

    for i in range(len(dna)):
        rand1 = randint(0, len(dna[0])-k)
        motifs_list.append(dna[i][rand1:rand1 + k])
    best_motifs = motifs_list

    for iteration in range(1, n):
        rand2 = randint(1, t-1)

        for motif in range(len(motifs_list)):
            if motifs_list[motif] != rand2:
                pre_profile_list.append(motifs_list[motif])
                profile = motif_profile(pre_profile_list, pseudocount=True)

        motifs_list[rand2] = most_probable_kmer(motifs_list[rand2], 8, profile)

        if motif_score(motifs_list) < motif_score(best_motifs):
            best_motifs = motifs_list

    return best_motifs


# ============================================================================
# =============================bioinformatics 1===============================
# it has some issues
# def neighbors(pattern, d):
#     """
#     :param pattern: a string of nucleotides
#     :param d: the number of permitted mismatches
#     :return: a set of sequences similar to pattern with at most,
#              d number of mismatches
#     """
#     if d == 0:
#         return {pattern}
#     elif len(pattern) == 1:
#         return {"A", "C", "G", "T"}
#
#     neighborhood = set()
#     suffix = pattern[d:]
#     suffix_neighbors = neighbors(suffix, d)
#     for text in suffix_neighbors:
#         if hamming_distance(suffix, text) < d:
#             for n in "ACGT":
#                 neighborhood.add(n + text)
#         else:
#             neighborhood.add(pattern[0] + text)
#
#     return neighborhood

# this one works

def neighbors(pattern, d):
    if d == 0:
        return [pattern]
    elif len(pattern) == 1:
        return "['A','C','G','T']"
    neighborhood = []
    suffix_pattern = pattern[1:]
    first_symbol_pattern = pattern[0]

    suffix_neighbors = neighbors(suffix_pattern, d)
    for text in suffix_neighbors:
        if hamming_distance(suffix_pattern, text) < d:
            for nucleotide in "ACGT":
                neighborhood.append(nucleotide+text)
        else:
            neighborhood.append(first_symbol_pattern + text)
    return neighborhood


# a correct approach copied from coursera (next 34 lines)
def frequent_words_with_mismatches(genome, k, d):
    start = time.process_time()
    approx_frq_words = []
    frequencies = defaultdict(lambda: 0)
    # all existent k-mers with d mismatches of current kmer in genome
    for index in range(len(genome) - k + 1):
        curr_kmer_and_neighbors = permute_motif_distance_times(genome[index: index + k], d)
        for kmer in curr_kmer_and_neighbors:
            frequencies[kmer] += 1

    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            approx_frq_words.append(kmer)
    end = time.process_time()
    print("Time:", end - start)
    return approx_frq_words


def permute_motif_once(motif):
    """
    Gets all strings within hamming distance 1 of motif and returns it as a
    list.
    """
    nucleotides = {"A", "C", "G", "T"}
    return list(set(itertools.chain.from_iterable([[
        motif[:pos] + nucleotide + motif[pos + 1:] for
        nucleotide in nucleotides] for
        pos in range(len(motif))])))


def permute_motif_distance_times(motif, d):
    working_set = {motif}
    for _ in range(d):
        working_set = set(itertools.chain.from_iterable(map(permute_motif_once, working_set)))
    return list(working_set)


def frequent_sequence_mismatches_reverse(text, k, d):
    """runs the frequent_words_with_mismatches twice,
       the latter time for text's reverse"""

    reversed_text = reverse_complement(text)
    freq_set_main = frequent_words_with_mismatches(text, k, d)
    freq_set_reverse = frequent_words_with_mismatches(reversed_text, k, d)
    main_reverse_combo = tools.space_sep(freq_set_main,) + tools.space_sep(freq_set_reverse)
    return tools.space_sep(main_reverse_combo)


def create_random_kmer(k=10):
    """
    :param k: length of the k-mer
    :return: a randomly generated string of k-mer with desired length
    """
    final_kmer = ""
    nuc = "ACTG"
    for i in range(k):
        final_kmer += nuc[randint(0, 3)]

    return final_kmer


# motif enumeration, mine doesn't work but the one after that
# which I got from the course's comments does

# mine
# def motif_enumeration(dna, k, d):
#     patterns = {0: [dna[0][e:e + k] for e in range(len(dna[0]) - k + 1)]}
#     for i in range(1, len(dna)):
#         # # ok up to this point
#         for kmer in patterns[0]:
#             k_index = rep.approximate_pattern_matching(dna[i], kmer, d)
#             # patterns = {dna[i][el:el+k] for el in k_index}
#             for el in k_index:
#
#                 try:
#                     patterns[i].append(dna[i][el:el+k])
#                 except KeyError:
#                     patterns[i] = [dna[i][el:el+k]]
#                 # when adding to the dict, use if-else with operand in,
#                 # to check and see whether they are in all the previous dict entries (use len to loop)
#
#     return patterns

# I didn't code this function myself
def motif_enumeration(dna, k, d):
    patterns = []
    for i in range(0, len(dna[0])-k+1):
        neighbor = neighbors(dna[0][i:i+k], d)

        for j in neighbor:
            count = 0
            for el in dna:
                for e in range(0, len(el)-k+1):
                    if hamming_distance(j, el[e:e+k]) <= d:
                        count += 1
                        break
            if count == len(dna):
                patterns.append(j)
    patts = []

    [patts.append(x) for x in patterns if x not in patts]
    result = ""
    for item in patts:
        result = result + " " + str(item)
    return result


# median string (copied)
def median_string(infile):

    # open infile
    with open(infile, 'r') as file:
        k = int(file.readline())
        dna = file.readlines()

        # iterate through each Dna string
        kmer_list = []
        median = []
        for line in dna:
            string = line.strip('\n')
            for i in range(len(string) - k+1):
                pattern = string[i:i+k]
                if pattern not in kmer_list:
                    kmer_list.append(pattern)

            # pattern that minimizes hamming distance
            distance = float('inf')
            for kmer in kmer_list:
                for i in range(len(string) - k+1):
                    if distance > hamming_distance(kmer, string[i:i+k]):
                        distance = hamming_distance(kmer, string[i:i+k])
                        median.append(kmer)
    return median


# median string problem (copied)

def pattern_string_distance(pattern, dna):
    k = len(pattern)
    distance = 0
    for seq in dna:
        length = len(seq)
        hd = float('inf')
        for i in range(length - k +1):
            hd_curr = hamming_distance(pattern, seq[i:i + k])
            if hd > hd_curr:
                hd = hd_curr
        distance += hd
    return distance
