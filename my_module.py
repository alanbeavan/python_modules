#!/usr/local/bin/python3

#python module:

def get_file_data(filename):
    """
    get_file_data
    
    stores the lines of the program with name filename as a list
    """
    import sys
    try:
        with open(filename) as in_file:
            lines = []
            for line in in_file:
                lines.append(line.rstrip("\n"))
        return(lines)

    except IOError as e:
        print("{}\n Error opening {}. Terminating program.".format(e, filename),
              file=sys.stderr)
        sys.exit(1)
    
    

def generate_sequence(length):
    """
    generate_sequence

    returnes a random sequence DNA of length given as argument
    """
    import random
    sequence = "";
    for i in range(0, length - 1):
        base = ""
        rand = random.randint(0,3)

        if(rand == 0):
            base = "A"
        elif(rand == 1):
            base = "C"
        elif(rand == 2):
            base = "G"
        elif(rand == 3):
            base = "T"

        sequence = sequence + base
    return sequence

def read_fasta(filename):
    """
    read_fasta

    reads a fasta file returning the seqs as a dictionary
    """
    import re
    import collections
    lines = get_file_data(filename)
    seqs = collections.OrderedDict()
    key = ""
    value = ""

    for line in lines:
        line = line.rstrip("\n")
        if re.search(">", line):
            if key:
                seqs[key] = value
                key = line[1:]
            else:
                key = line[1:]
            value = ""
        else:
            value = value + line
    seqs[key] = value
    return(seqs)

def read_phylip_unpartitioned(filename):
    """
    read_phylip_unpartitioned

    Reads a phylip file and return a dictionary of the names and sequences
    Labels and seqs must be separated with whitespace.
    There must not be spaces in the sequence label
    """
    import re
    lines = get_file_data(filename)
    seqs = {}
    key = ""
    lengths = lines[0].strip()
    length = lengths.split()[1]


    for line in lines:
        if re.match('^\S+', line):
            fields = line.split()
            if len(fields) > 1:
                key = fields[0]

                #Sometimes sequences are separated by whitespace - concatinate fields[1-n]
                seq = fields[1]
                if(len(fields) > 2):
                    for i in range(2, len(fields)):
                        seq += fields[i]

                if key in seqs:
                    seqs[key] += seq
                else:
                    seqs[key] = seq

    to_delete = []
    for key in seqs:
        if len(seqs[key]) != int(length):
            to_delete.append(key)


    for key in to_delete:
        del seqs[key]

    return(seqs)

def reverse_complement(string):
    """
    reverse_complement

    returns the reverse complement of a DNA string provided as argument
    """
    i = len(string)-1
    reverse_complement = ""
    while i >= 0 :
        if string[i] == "A":
            reverse_complement = reverse_complement + "T"
        elif string[i] == "C":
            reverse_complement = reverse_complement + "G"
        elif string[i] == "G":
            reverse_complement = reverse_complement + "C"
        elif string[i] == "T":
            reverse_complement = reverse_complement + "A"
        else:
            print("characters that are not A,C,G or T are present in the string\n\n")
            break
        i = i - 1
    return(reverse_complement)

def transcribe(string):
    """
    transcribe

    return the transcribed RNA given the DNA string provided as argument.
    """
    RNA = ""
    for base in string:
        if base == "T":
            RNA = RNA + "U"
        elif base in ["A","C","G"]:
            RNA = RNA + base
        else:
            print("characters that are not A, C, G or T are present in the string\n\n")
            return(0)
    return(RNA)

def reverse_transcribe(string):
    """
    reverse_transcribe

    return the trancscipt of the reverse complement of the DNA string provided.
    """
    DNA = ""
    for base in string:
        if base == "U":
            DNA = DNA + "T"
        elif base in ["A","C","G"]:
            DNA = DNA + base
        else:
            print("characters that are not A, C, G or U are present in the string\n\n")
            return(0)
    return(DNA)

def root_at_midpoint(self):
    """
    root_at_midpoint

    takes a tree and retuns the tree rooted at it's midpoint.
    The tree has to be some module's version but I can't remember which lol.
    I also believe this is not all my original work.
    """
    # Identify the largest pairwise distance
    max_distance = 0.0
    tips = self.get_terminals()
    for tip in tips:
       self.root_with_outgroup(tip)
       new_max = max(self.depths().items(), key=lambda nd: nd[1])
       if new_max[1] > max_distance:
           tip1 = tip
           tip2 = new_max[0]
           max_distance = new_max[1]
    self.root_with_outgroup(tip1)
    # Depth to go from the ingroup tip toward the outgroup tip
    root_remainder = 0.5 * (max_distance - (self.root.branch_length or 0))
    assert root_remainder >= 0
    # Identify the midpoint and reroot there.
    # Trace the path to the outgroup tip until all of the root depth has
    # been traveled/accounted for.
    for node in self.get_path(tip2):
       root_remainder -= node.branch_length
       if root_remainder < 0:
           outgroup_node = node
           outgroup_branch_length = -root_remainder
           break
    else:
       raise ValueError("Somehow, failed to find the midpoint!")
    self.root_with_outgroup(outgroup_node,
                           outgroup_branch_length=outgroup_branch_length)

def hamming_distance(pair):
    """
    hamming_distance

    calculate and return the hamming distance between 2 strings of the same
    length.
    """
    hamming_distance = 0
    i = 0
    while i < len(pair[0]):
        if pair[0][i] != pair[1][i]:
            hamming_distance += 1
        i += 1
    return hamming_distance


#returns the GC content from a DNA string - does not deal with non ACTG characters and doesn't
#give an error. I will change this if necessary
def gc_content(seq):
    """
    gc_content

    returns the 
    """
    gc = (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c"))/len(seq)
    return gc

def translate_DNA(seq):
    """
    translate_DNA

    returns the protein coded for by the first ORF in the provided DNA string.
    """
    map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L","TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S","TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP","TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W","CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L","CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R","ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M","ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T","AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R","GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V","GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    i=0
    protein =""
    while i+2 < len(seq):
        if(map[seq[i:i+3]] == "STOP"):
            break
        else:
            protein = protein + map[seq[i:i+3]]
        i+=3
    return protein

def translate_RNA(seq):
    """
    translate_RNA

    translates RNA to protein in ORF 1 of the provided sequence
    """
    map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L","UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S","UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP","UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W","CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L","CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R","AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M","ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T","AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R","GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V","GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    i=0
    protein =""
    while i+2 < len(seq):
        if(map[seq[i:i+3]] == "STOP"):
            break
        else:
            protein = protein + map[seq[i:i+3]]
        i+=3
    return protein

def occurances(seq, motif):
    """
    occurances
    
    given a seq and motif, return the indices of the motif in the sequence
    """
    i = 0
    occurances = []
    while i < len(seq) - len(motif):
        if (seq[i:i+len(motif)] == motif):
            occurances.append(i)
        i += 1
    return occurances

def counts(seqs):
    """
    counts

    From a set of alignes seqs, return the counts of each residue at each locus.
    Output is a pandas dataframe
    """
    import numpy as np
    import pandas as pd
    chars = ""
    length = 0


    for key in seqs.keys():
        chars = chars + "".join(set(seqs[key]))
        length = len(seqs[key])
        seqs[key] = list(seqs[key])
    chars = sorted(set(chars))


    matrix = np.arange(length * len(chars)).reshape(len(chars), length)
    df = pd.DataFrame(matrix, index = chars)

    counts_at_site = dict.fromkeys(chars, 0)
    alignment = pd.DataFrame.from_dict(seqs)
    for i in range(0, alignment.shape[0]):
        for j in range(0, alignment.shape[1]):
            counts_at_site[alignment.iloc[i,j]] += 1

        for key in counts_at_site.keys():
            df.iloc[df.index.get_loc(key), i] = counts_at_site[key]


        counts_at_site = dict.fromkeys(chars, 0)

    return df

def consensus_seq(df):
    """
    consensus_seq

    Generates and returns a consensus sequence from counts in a pandas dataframe
    """
    import pandas as pd

    seq = ""
    for i in range(0,df.shape[1]):
        seq = seq + df[i].idxmax()
    return seq

def get_all_substrings(input_string):
    """
    get_all_substrings

    returns all the possible substrings for a string
    """
    length = len(input_string)
    return [input_string[i:j+1] for i in range(length) for j in range(i,length)]

def longest_common_substring(data):
    """
    longest_common_substring

    returns the longest common substring from a list of strings.
    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr

def is_substr(find, data):
    """
    is_substr

    return boolean true is find is in data.
    """
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True

def get_uniprot_sequence(id):
    """
    get_uniprot_sequence

    obtain and return the sequence associated with a uniport id.
    """
    from requests import get

    seq = get("http://www.uniprot.org/uniprot/" + id + ".fasta").content.decode("utf-8")
    lines = seq.split("\n")
    protein = {}
    protein[">".join(lines[0].split(">")[1:])] = "".join(lines[1:])

    return protein

def forward_orfs(seq):
    """
    forward_orfs

    Return all possible  forward ORFs in a a sequence. i.e. ending in stop codons
    """
    map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L","TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S","TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP","TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W","CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L","CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R","ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M","ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T","AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R","GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V","GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    proteins = []
    for i in range(0,3):
        j = i
        while j < len(seq) - 2:
            #print(seq[j:j+3])
            if seq[j:j+3] == "ATG":
                orf = "ATG"
                k = j+3
                while k < len(seq) - 2:
                    if map[seq[k:k+3]] == "STOP":
                        orf = orf + seq[k:k+3]
                        proteins.append(translate_DNA(orf))
                        break
                    else:
                        orf = orf + seq[k:k+3]
                    k += 3
            j += 3
    return proteins

def proteins(seq):
    """
    proteins

    Returns all the proteins coded for by ORFs in the forward direction.
    """
    proteins = forward_orfs(seq)
    for protein in forward_orfs(reverse_complement(seq)):
        proteins.append(protein)
    return proteins

def permutations(iterable, r=None):
    """
    permutations

    returns a list of permutations for an iterable item and length
    permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    permutations(range(3)) --> 012 021 102 120 201 210
    """
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

def fasta_to_phylip(infile, outfile):
    """
    fasta_to_phylip

    Reads a fasta file and writes a phylip version.
    """
    from Bio import AlignIO
    alignments = AlignIO.parse(infile, "fasta")
    AlignIO.write(alignments, outfile, 'phylip-relaxed')

def write_matrix(dict, filename = "matrix.phy"):
    """
    write_matrix

    writes a set of seqs as a phylip file.
    """
    f = open(filename, "w")
    f.write(str(len(dict.keys())) + "\t")
    i = 1
    for key in dict.keys():
        if i == 1:
            f.write(str(len(dict[key])) + "\n")
        f.write(key + "\t" + dict[key] + "\n")
        i += 1
    f.close()

def string_permutations(objects, n, acc='', perms=[]):
    """
    string_permutations

    Returns the permutations of lenght n in a string.
    """
    if n == 0:
        perms.append(acc)
    else:
        for char in chars:
            string_permutations(objects, n - 1, acc + char, perms)

    return(perms)

def combinations(iterable, r):
    """
    combinations

    Returns all the combinations of length r from an iterable.
    combinations('ABCD', 2) --> AB AC AD BC BD CD
    combinations(range(4), 3) --> 012 013 023 123
    """
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in list(reversed(range(r))):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def combinations_with_replacement(iterable, r):
    """
    combinations_with_replacement

    Retunrs a list of combinations with replacement.
    combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    """
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)

def product(*args, repeat=1):
    """
    product

    return the list of combinations of some lists.
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    """
    combs = []

    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        combs.append(tuple(prod))
    return(combs)

def TsTvRatio(before, after):
    """
    TsTvRatio

    return the transition transversion ratio of 2 DNA sequences
    """
    purines = ["A", "G"]
    pyrimidines = ["C", "T"]

    if len(before) != len(after):
        print("the sequences are not aligned. Try again with seqs of the same length")
        return 0
    Tv = 0
    Ts = 0
    for i in range(0,len(before)):
        if before[i] != after[i]:
            if before[i] in purines:
                if after[i] in purines:
                    Ts += 1
                elif after[i] in pyrimidines:
                    Tv += 1
            elif before[i] in pyrimidines:
                if after[i] in purines:
                    Tv += 1
                elif after[i] in pyrimidines:
                    Ts += 1
            else:
                print("non - ACTG characters present - check your shit")
                return 0
    return float(Ts)/float(Tv)

def add_branch_lengths(tree_file, table):
    """
    add_branch_lengths

    adds branch lengths to a nwk tree in the format label:length.
    """
    import subprocess
    tree = get_file_data(tree_file)[0]
    table = get_file_data(table)

    for line in table:
        fields = line.split("\t")
        subprocess.call("sed -i .bak \"s/" + fields[0] + ":[[:digit:]]*\.*[[:digit:]]*/" + fields[0] + ":" + fields[1] + "/\" " + tree_file, shell = True)

