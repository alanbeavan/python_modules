#!/usr/local/bin/python3

#python module:

#get_file_data - return a list of lines from a file "filename".
def get_file_data(filename):
    """get_file_data.
    
    
    stores the lines of the program with name filename as a list
    """
    try:
	with open(file) as in_file:
            lines = []
            for line in f:
                lines.append(line.rstrip("\n"))
        return(lines)

    except IOError as e:
	print("{}\n Error opening {}. Terminating program.".format(e, file),
              file=sys.stderr)
        sys.exit(1)
    
    

#Generate a random sequence of user defined length
def generate_sequence(length):
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

#read_fasta - read a fasta file returning a dictionary of names of seq keys and seq values
def read_fasta(filename):
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

#read_phylip - read a phylip file and return a dictionary of the names and sequences - labels and seqs must be separated with whitespace. There must not be spaces in the sequence label
def read_phylip_unpartitioned(filename):
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

#Compute_bic - computes the bayesian information criterion for a set of clusters form 0-X specified by the user.
#input - kmeans - kmeans at k, X = multidementional np array of data points
#I don't think this works!
def compute_bic(kmeans,X):
    from sklearn import cluster
    from scipy.spatial import distance
    import sklearn.datasets
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    # assign centers and labels
    centers = [kmeans.cluster_centers_]
    labels  = kmeans.labels_
    #number of clusters
    m = kmeans.n_clusters
    # size of the clusters
    n = np.bincount(labels)
    #size of data set
    N, d = X.shape

    #compute variance for all clusters beforehand
    cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]],
             'euclidean')**2) for i in range(m)])

    const_term = 0.5 * m * np.log(N) * (d+1)

    BIC = np.sum([n[i] * np.log(n[i]) -
               n[i] * np.log(N) -
             ((n[i] * d) / 2) * np.log(2*np.pi*cl_var) -
             ((n[i] - 1) * d/ 2) for i in range(m)]) - const_term

    return(BIC)

#Return the reverse complement of a DNA string
def reverse_complement(string):
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

#Count the nucleotides of a DNA sequence, returning a dictionary
def count_nucleotides(string):
    nucleotides = ["A","C","G","T"]
    counts = {    N : 0 for N in nucleotides}

    for base in string:
        if base.upper() in nucleotides:
            counts[base.upper()] += 1
        else:
            print("characters that are not A, C, G or T are present in the string\n\n")
            return(0)

    return(counts)

#Count the amino acids of a protein sequence, returning a dictionary
def count_amino_acids(string):
    AAs = ["A","R","N","D","B","C","E","Q","Z","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    counts = {    AA : 0 for AA in AAs}

    for base in string:
        if base.upper() in AAs:
            counts[base.upper()] += 1
        else:
            print("characters that do not represent standard amino acids are present in the string\n\n")
            return(0)

    return(counts)

#Return the transcribed DNA sequence as RNA string
def transcribe(string):
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

#Return the reverse-transcribed RNA sequence as DNA string
def reverse_transcribe(string):
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

#Root the tree at the midpoint of the two most distant taxa in place - i mostly took this from some other module
def root_at_midpoint(self):
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

#Calculate the hamming distance between 2 sequences (which can of course be any strings)
def hamming_distance(pair):
    hamming_distance = 0
    i = 0
    while i < len(pair[0]):
        if pair[0][i] != pair[1][i]:
            hamming_distance += 1
        i += 1
    return hamming_distance

#Take as input a pair of integers from a file. The first is the number of generaions
#the second is the number of new pairs generated per generation
#rabits become sexually mature after 1 month Returns the number of pairs. This doesn't have to
#be rabbits. It is apparently an important style of sequence so might come in handy
def rabbits(starting_young, starting_mature, generations, litter):
    young_pairs = starting_young
    mature_pairs = starting_mature
    i=1
    while i < int(generations):
        new_mature_pairs = young_pairs
        new_young_pairs = mature_pairs * int(litter)
        mature_pairs = mature_pairs + new_mature_pairs
        young_pairs = new_young_pairs
        #print(str(i) + "\tyoung " + str(young_pairs) + "\told" + str(mature_pairs))
        i += 1
    return mature_pairs + young_pairs

#Simulate the birth of population size of rabbits, where 1 pair takes 1 unit of time to become mature
#each pair breeds every unit of time given that they are sexually mature and produces a litter of 1
#the rabbits die every m months. Returns the number of rabbit pairs after m months.
def mortal_rabbits(starting_young, starting_mature, litter, n, m):
    young_pairs = starting_young
    mature_pairs = starting_mature
    to_die = [1]

    i=1
    while i < int(n):
        new_mature_pairs = young_pairs
        new_young_pairs = mature_pairs * int(litter)
        to_die.append(new_young_pairs)
        if i < m:
            mature_pairs = mature_pairs + new_mature_pairs
            young_pairs = new_young_pairs
        else:
            mature_pairs = mature_pairs + new_mature_pairs - to_die[i-m]
            young_pairs = new_young_pairs
        #print(str(i) + "\tyoung " + str(young_pairs) + "\told" + str(mature_pairs))
        i += 1
    return mature_pairs + young_pairs

#returns the GC content from a DNA string - does not deal with non ACTG characters and doesn't
#give an error. I will change this if necessary
def gc_content(seq):
    gc = (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c"))/len(seq)
    return gc

#take a list of homozygous dominant, heterozygous, and homozygous recessive individual counts
#Return the chance of offspring having the phenotype detirmined by the dominant allele
def dominant_probability(freqs):
    total = sum(freqs)
    #freqs[:] = [x/total for x in freqs]
    A1A1, A1A2, A2A2 = freqs

    prob = 1 - (A2A2/total * ((A2A2-1)/(total-1) + (A1A2/(total-1))/2) + ((A1A2/total)/2) * (A2A2/(total-1) + ((A1A2-1)/(total-1))/2))
    return prob

#translates DNA into protein starting from the first codon (ORF1)
def translate_DNA(seq):
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

#translates RNA into protein starting from the first codon (ORF1)
def translate_RNA(seq):
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

#given a sequence and motif, give the indecies of all occurences of the motif in the seq where
#each index correesponds to a sequence the same length as the motif. The position is 0 indexed
def occurances(seq, motif):
    i = 0
    occurances = []
    while i < len(seq) - len(motif):
        if (seq[i:i+len(motif)] == motif):
            occurances.append(i)
        i += 1
    return occurances

#Takes as input a dictionary of seqs (may convert this to including also lists) Seqs must be aligned
#returns a numpy matrix of the frequency of each base at each site in the alignment Currently, this
#cannot deal with gaps
def counts(seqs):
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

#From a pandas df, return the consensus sequence from the counts. Ties result in one or the other base
def consensus(df):
    import pandas as pd

    seq = ""
    for i in range(0,df.shape[1]):
        seq = seq + df[i].idxmax()
    return seq

#find overlaps between sequences of given length. Returns list of space separated pairs connected
#by edge (labels not seqs)
def find_overlaps(seqs, overlap):
    edges = []
    for key in seqs:
        tail = seqs[key][-3:]
        for key2 in seqs:
            if key2 != key:
                if(tail == seqs[key2][:3]):
                    edges.append(" ".join([key,key2]))

    return edges

#get all the possible substrings for a string
def get_all_substrings(input_string):
  length = len(input_string)
  return [input_string[i:j+1] for i in range(length) for j in range(i,length)]

#find the longest common substring from a list of strings
def longest_common_substring(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr

#return True if find is a substring of data
def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True

#n_successes - return the number of possible ways to get k successes in a bernouli distriubtion
#This timesed by the probability of getting there in 1 specific way gives the probability of
#getting this many successes
def n_combinations(n, r):
    #The distribution is symmetrical
    r = min(r, n-r)
    numer = n
    for i in range(n-1, n-r, -1):
        numer = numer * i

    denom = 1
    for i in range(2, r+1):
        denom = denom * i
    return numer/denom

#return the probability of getting K or more successes with a bernouli distriubtion with chance
#of success p in a sample of N
def n_or_more(K, N, p):
    cumprob = 0
    for i in range(K, N+1):
        combinations = n_combinations(N, i)
        prob = combinations * (p ** i) * ((1-p) ** (N-int(i)))
        cumprob += prob
        print(str(i) + "\t" + str(prob))
    return cumprob

#obtains the sequence and header as a key value pair for a uniprot id
def get_uniprot_sequence(id):
    from requests import get

    seq = get("http://www.uniprot.org/uniprot/" + id + ".fasta").content.decode("utf-8")
    lines = seq.split("\n")
    protein = {}
    protein[">".join(lines[0].split(">")[1:])] = "".join(lines[1:])

    return protein

#Obtains positions of a motif (in the form of a regular expression fit for python re) in a seq
#It does not do any computation of the regex. Ensure this is correct when you pass it to the funciton
def find_motif(seq, motif):
    import re
    positions = []
    for i in range(0, len(seq)-len(motif)):
        if re.match(motif, seq[i:i+len(motif)]):
            positions.append(i+1)
    return positions

#
def n_mRNA_possible(seq):

    num_per_AA = {"A":4,
                "R":6,
                "D":2,
                "N":2,
                "C":2,
                "E":2,
                "Q":2,
                "G":4,
                "H":2,
                "I":3,
                "L":6,
                "K":2,
                "M":1,
                "F":2,
                "P":4,
                "S":6,
                "T":4,
                "W":1,
                "Y":2,
                "V":4}
    n = 1
    x = 0
    for AA in seq:
        n = n * num_per_AA[AA]
        x = x * num_per_AA[AA]
        if n > 1000000:
            x += int(n/1000000)
            n = n % 1000000

    #Also 3 STOP codons so the number of sequences needs to be times by 3 at the end
    n = n*3
    x = x*3
    if n > 1000000:
        x += int(n/1000000)
        n = n % 1000000

    return( x, n)

#return all the possible proteins from a DNA sequence in all 6 ORFs. They must end in stop codons
def forward_orfs(seq):
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

#get the protein sequences of all the ORF of a DNA string and it's reverse compliment
def proteins(seq):
    proteins = forward_orfs(seq)
    for protein in forward_orfs(reverse_complement(seq)):
        proteins.append(protein)
    return proteins

#return the number of permutations and a list of those permutations of length n where the contents is n
#unique identifiers, here numbers from 1 to n. The results of this could be used with a dictionary
#to be translated into the data we are looking at. The return value is a list of strings of space separated
#lists. One could also just use itertools which is in the standard libraries. This is the Heap's algorithm
def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
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

#From a protein sequence, return its monoisotonic mass
def protein_mass(seq):
    masses = {"A":71.03711,"C":103.00919,"D":115.02694,"E":129.04259,"F":147.06841,
    "G":57.02146,"H":137.05891,"I":113.08406,"K":128.09496,"L":113.08406,
    "M":131.04049,"N":114.04293,"P":97.05276,"Q":128.05858,"R":156.10111,
    "S":87.03203,"T":101.04768,"V":99.06841,"W":186.07931,"Y":163.06333}

    mass = 0
    for aa in seq:
        mass += masses[aa]
    return mass

#Writes a phylip file from a fasta file
def fasta_to_phylip(infile, outfile):
    from Bio import AlignIO
    alignments = AlignIO.parse(infile, "fasta")
    AlignIO.write(alignments, outfile, 'phylip-relaxed')

#from a dictionary of species (or anything really) with their values for each site (or trait or whatever)
#write a matrix in phylip format
def write_matrix(dict, filename = "matrix.phy"):
    f = open(filename, "w")
    f.write(str(len(dict.keys())) + "\t")
    i = 1
    for key in dict.keys():
        if i == 1:
            f.write(str(len(dict[key])) + "\n")
        f.write(key + "\t" + dict[key] + "\n")
        i += 1
    f.close()

#Return a string that is string without pattern - basically just re.sub
def splice_string(string, pattern):
    import re
    string = re.sub(pattern, "", string)
    return string

#return the combinations of chars of length n
def string_permutations(objects, n, acc='', perms=[]):
    #Recusiveboi
    if n == 0:
        perms.append(acc)
    else:
        for char in chars:
            string_permutations(objects, n - 1, acc + char, perms)

    return(perms)

#Return the number of combinations of r elements from a list - iterable. Slightly modified from itertools
def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
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
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
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
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    combs = []

    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        combs.append(tuple(prod))
    return(combs)

#from a list of numbers (or alphabetical strings actually), return the longest ascending or
#descending orderred list that can be obtained by removing elements from the list
def longest_orderred_list(seq, reverse = False):
    P = [None] * len(seq)
    M = [None] * len(seq)

    L = 1
    M[0] = 0
    for i in range(1, len(seq)):
        if(reverse):
            lo = 0
            hi = L
            if seq[M[hi - 1]] > seq[i]:
                j = hi
            else:
                while hi - lo > 1:
                    mid = (hi + lo) // 2
                    if seq[M[mid - 1]] > seq[i]:
                        lo = mid
                    else:
                        hi = mid

                j = lo
            P[i] = M[j - 1]
            if j == L or seq[i] > seq[M[j]]:
                M[j] = i
                L = max(L, j + 1)

        else:
            lo = 0
            hi = L
            if seq[M[hi - 1]] < seq[i]:
                j = hi
            else:
                while hi - lo > 1:
                    mid = (hi + lo) // 2
                    if seq[M[mid - 1]] < seq[i]:
                        lo = mid
                    else:
                        hi = mid

                j = lo
            P[i] = M[j - 1]
            if j == L or seq[i] < seq[M[j]]:
                M[j] = i
                L = max(L, j + 1)

    result = []
    pos = M[L - 1]
    for k in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return (result[::-1])

#from a set of sequences, find at print the shortest subsequence that contains them all
#This will sometimes be wrong - i did it to get better at the programming
#This now returns a list of contigs if just 1 cannot be constructed
def shortest_superstring_greedy(seqs):
    while(len(seqs) > 1):
        largest = len(sorted(seqs, key=len)[-2])
        overlaps = {}
        for combination in combinations(range(0, len(seqs)),2):
            if largest_overlap(seqs[combination[0]], seqs[combination[1]]):
                overlaps.update(largest_overlap(seqs[combination[0]], seqs[combination[1]]))

        #If no more contigs can be constructed
        if len(overlaps) == 0:
            break

        strings = max(overlaps, key = overlaps.get).split("\t")
        seqs.remove(strings[0])
        #if the 2 strings are identical, we need to check strings[1] exists
        if strings[1] in seqs:
            seqs.remove(strings[1])
        seqs.append(strings[2])

    if len(seqs) == 1:
        return seqs[0]
    else:
        return seqs

#find the largest overlap and return thos seqs concatinated with the others.
#return a dictionary with key being "\t" separated values of the 2 seqs and the concatinated seq
#and value of the length of the overlap
#return 0 if there is no overlap
def largest_overlap(seq0, seq1):
    cat = {}
    len0 = len(seq0)
    len1 = len(seq1)
    i =  min(len0, len1)
    #The maximum overlap is the shortest of the 2 lengths
    while i > 0:
        #compare the last i characters in seq 0 with the first i characters in seq 1
        if seq0[len0-i:len0] == seq1[0:i]:
            if i == len(seq0):
                cat = {"\t".join([seq0, seq1, max([seq0, seq1], key = len)]) : i}
            else:
                cat = {"\t".join([seq0, seq1, seq0 + seq1[i:]]) : i}
            break
        #else if the last i characters in seq 1 with the first characters in seq 0
        elif seq1[len1-i:len1] == seq0[0:i]:
            if i == len(seq1):
                cat = {"\t".join([seq0, seq1, max([seq0, seq1], key = len)]) : i}
            else:
                cat = {"\t".join([seq0, seq1, seq1 + seq0[i:]]) : i}
            break
        i -= 1
    if(len(cat) > 0):
        return cat
    else:
        return 0

#return the indicies of matches for char in string
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

#Provided an equal numbers of A and U, and of C and G, return a list of the
#perfect graphs in the form of lists of tupples representing the edges, with
#the numbers being the index in the sequence This has O(n^2) so takes ages with
#large datasets
def perfect_RNA_graphs(RNA):
    import re
    if re.search("[^ACGU]", RNA):
        print("Not a valid RNA string due to non-RNA nucleotide characters")
        return 0
    indicies = {"A":find(RNA, "A"), "C":find(RNA,"C"), "G":find(RNA,"G"),"U":find(RNA,"U")}
    if len(indicies["A"]) != len(indicies["U"]) or len(indicies["C"]) != len(indicies["G"]):
        print("No perfect graphs can be made from this RNA string")
        return 0
    matches = {"A":"U", "C":"G", "G":"C", "U":"A"}

    AtoU = all_possible_tuples(indicies["A"], indicies["U"])
    CtoG = all_possible_tuples(indicies["C"], indicies["G"])

    graphs = list(product(AtoU, CtoG))

    for i in range(len(graphs)):
        graphs[i] = graphs[i][0] + graphs[i][1]

    return graphs

#return all possible ways of matching equally sized lists without replacement
def all_possible_tuples(list0, list1):
    zips = []
    for permutation in permutations(list1, len(list1)):
        zips.append(list(zip(permutation, list0)))
    return list(zips)

#return just the number of perfect graphs
def n_perfect_graphs(RNA):
    from math import factorial
    import re
    if re.search("[^ACGU]", RNA):
        print("Not a valid RNA string due to non-RNA nucleotide characters")
        return 0
    indicies = {"A":find(RNA, "A"), "C":find(RNA,"C"), "G":find(RNA,"G"),"U":find(RNA,"U")}
    if len(indicies["A"]) != len(indicies["U"]) or len(indicies["C"]) != len(indicies["G"]):
        print("No perfect graphs can be made from this RNA string")
        return 0

    AU = 0
    GC = 0
    for nt in RNA:
        if nt == 'A':
            AU += 1
        elif nt == 'G':
            GC += 1
    return factorial(AU) * factorial(GC)

#return the number of partial permutations of a set of n containing k
#return the number of millions and the modulo value 1,000,000
def n_partial_permutaions(n,k):
    mils = 0
    mod = 1
    #the number is n * n-1 * n-2 * ... * n-k i think
    for i in range(0,k):
        mils = mils * (n-i)
        mod = mod * (n-i)

        mils += int(mod/1000000)
        mod = mod%1000000
    return mils,mod

#return the likelihood of getting a sequence from a randomly modelled seq of the
#same length given only the probability that a base is g or c.
def likelihood_from_gc(gc, seq):
    import math
    probability = {"A":(1-gc)/2,
                    "C":gc/2,
                    "G":gc/2,
                    "T":(1-gc)/2}

    likelihood = 1
    for base in seq:
        likelihood = likelihood * probability[base]
    return math.log(likelihood, 10)

#return the signed permutations of a sequence of 1:n I might modify this to take any
#list
def signed(n):
    positive = list(range(1,n+1))
    perms = permutations(positive, len(positive))
    print(perms)
    signed = []
    for perm in perms:
        perm = list(perm)
        for i in range(0, len(positive)+1):
            combs = combinations(perm, i)
            for comb in combs:
                copy = perm.copy()
                for k in comb:
                    copy[k-1] = -copy[k-1]
                signed.append(list(copy))
    return signed

#return a list of indicies of substring in seq if substring is a subsequence of seq
def subsequence_indicies(seq, substring):
    i = 0
    indicies = []
    for j in range(0, len(seq)):
        if substring[i] == seq[j]:
            indicies.append(j+1)
            i += 1
            if i == len(substring):
                break
    if len(indicies) == len(substring):
        return indicies
    else:
        print("substring not a subsequence of seq")
        return 0

#return the transition transversion ratio of 2 DNA sequences
def TsTvRatio(before, after):
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

#adds branch lengths to a nwk tree in the format label:length
def add_branch_lengths(tree_file, table):
    import subprocess
    tree = get_file_data(tree_file)[0]
    table = get_file_data(table)

    for line in table:
        fields = line.split("\t")
        subprocess.call("sed -i .bak \"s/" + fields[0] + ":[[:digit:]]*\.*[[:digit:]]*/" + fields[0] + ":" + fields[1] + "/\" " + tree_file, shell = True)
