# Implement PatternToNumber
dna=input()
dna_re=dna[::-1]
result=0
i=0
while i<len(dna):
    if dna_re[i]=='A':
        pass
    elif dna_re[i]=='C':
        result+=4**i
    elif dna_re[i]=='G':
        result+=2*4**i
    else:
        result+=3*4**i
    i+=1
print(result)

# Generate the Frequency Array of a String
from itertools import *
seq = "TATACTGTGGGGCAGTGATTACCCCCAGCAAGGGGGTCCAGTGGTATGCAGCTCCACCGGCAAACAAGGAAGATCAGATCGACAGGCGTAACAATGCACCTCGTGGCGACTAGACATCCCATCGAGCGTATTGTTATCTTGTTCACATGAGACTCAGGTGGACTACTCCTACTGCAGGCGTGCAAGTTGTCTATCTGCCCTTGAGGCGGGCAGTAGTAGCCATGCCGCGTCTCCGAGTTATCAGGTCACATGACCTGGTTTTCCGTACCGCTAGGACACTCTTTAATAGAATATGAAGCCAATTCAGTATCCCACCAAAATCATGATGGTGCCGGCGCCGATGTTTTAGGACGTGCAGATACAACCAGCATGAATTCGATTTAACTAACTCGAGCTTAACCAGAAAGCGGTCGACTGATAAGTGATACAACTTGTTAACCACCCCGTTCGCCGTTGCTTATCAAAACATTGCTTAATGCGTTCTGTGGGCTTAATCAGCCAGAAGGCAATCGAGTAGATCGGTGCTGGAACGGACTGCAAAGAGTGAATAGAACGTGTACTCAAGCACAGAGCACCTTAAGTCAGTGTTGGCCACTCTGTAGGCGAAAACTCTAAACCGACCGTAGGTATTTTTAAGGAATTACGCCCATCGGGCCCCCAAAGCGACCGCGGCTGTTGCATGCATGCTATCGCGAGGGTACAAGGGAT"
n = 7
def gen_all_k_mer(n):
    list_result = []
    list = [i for i in product(["A","T","C","G"], repeat=n)]
    str_list = []
    for tuple in list:
        string = ""
        for i in tuple:
            string += i
        str_list.append(string)
    return sorted(str_list)
kmers = gen_all_k_mer(n)
freq_list = []
for kmer in kmers:
    counter = 0
    for i in range(len(seq)):
        if seq[i:i+n] == kmer:
            counter += 1
    freq_list.append(counter)
for i in freq_list:
    print(i,end=" ")

# Find the Most Frequent Words in a String
from itertools import *
seq = "GTCTCTCTAAAGATGGTACGCCCCTTAGAACACAAGTCTCTCTAGTCTCTCTAAAGATGGTAAGATGGTAGAACACAAAAGATGGTACGCCCCTTAAGATGGTAAGATGGTAAGATGGTAGAACACAAAAGATGGTAGAACACAAAAGATGGTGTCTCTCTAGTCTCTCTAAGAACACAAGTCTCTCTAACGCCCCTTAAGATGGTAGAACACAAGTCTCTCTAGTCTCTCTAGATAACATGAAGATGGTAGAACACAAAAGATGGTAAGATGGTAAGATGGTGTCTCTCTAAGAACACAAACGCCCCTTGATAACATGAGAACACAAGTCTCTCTAAGAACACAAAGAACACAAAAGATGGTGATAACATGGTCTCTCTAACGCCCCTTGTCTCTCTAGTCTCTCTAGTCTCTCTAACGCCCCTTGATAACATGAGAACACAAGTCTCTCTAACGCCCCTTAAGATGGTAGAACACAAGATAACATGAGAACACAAGATAACATGGATAACATGGATAACATGGTCTCTCTAACGCCCCTTAAGATGGTAAGATGGTGATAACATGGTCTCTCTAACGCCCCTTGTCTCTCTAGATAACATGGTCTCTCTAGTCTCTCTAACGCCCCTTGATAACATGGATAACATGACGCCCCTTAGAACACAAAGAACACAAAAGATGGTACGCCCCTTAGAACACAAAAGATGGTGATAACATGACGCCCCTTACGCCCCTTAAGATGGTAAGATGGTACGCCCCTTACGCCCCTTGATAACATGAGAACACAAAAGATGGTGTCTCTCTAGATAACATGGTCTCTCTAAAGATGGTAGAACACAAAAGATGGTAAGATGGTAAGATGGTAGAACACAAAGAACACAAGTCTCTCTAGTCTCTCTA"
n = 12
dict = {}
for i in range(len(seq)-n+1):
    if seq[i:i+n] not in dict:
        dict[seq[i:i+n]] = 1
    else:
        dict[seq[i:i+n]] += 1
for i in dict:
    if dict[i] == max(dict.values()):
        print(i,end=" ")

# Find the Most Frequent Words with Mismatches in a String
seq = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k = 4
d = 1
dict = {}

# Compute the Probability of a Hidden Path: UNSOLVED
path = "AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB"
prob = 1
for i in path:
    if i == "A":
        prob *= (0.194+0.273)*0.5
    if i == "B":
        prob *= (0.806+0.727)*0.5
prob

# Generate the Theoretical Spectrum of a Cyclic Peptide: UNSOLVED
from itertools import *
def read_data(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split()[0])
    infile.close()
    return data
data = read_data("test")
protein = data[0]
list = []
for i in protein:
    list.append(i)

sublists = []

def read_data_list(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split())
    infile.close()
    return data
protein_mass = read_data_list("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[amino[0]] = eval(amino[1])

mass_list = []
for sublist in sublists:
    mass = 0
    for amino in sublist:
        mass += mass_dict[amino]
    mass_list.append(mass)
mass_list.sort()
for i in mass_list:
    print(i,end=" ")

# Compute the Number of Times a Pattern Appears in a Text
data = read_data("ba1a")
seq = data[0]
pattern = data[1]
count = 0
for i in range(len(seq)):
    if seq[i:i+len(pattern)] == pattern:
        count += 1
print(count)

# Compute the Probability of an Outcome Given a Hidden Path
path = "BAABAABBABBAABBABAAAABAAAAAABBAAAABBABAABAAAABABAA"
outcome = "yyzxxxxzxyyzxxzxzyyxyxxxyxxxzzxxzxxxzzxxzyzxyzxzzx"
prob = 1
for i in range(len(path)):
    if path[i] == "A":
        if outcome[i] == "x":
            prob *= 0.337
        elif outcome[i] == "y":
            prob *= 0.335
        else:
            prob *= 0.328
    else:
        if outcome[i] == "x":
            prob *= 0.637
        elif outcome[i] == "y":
            prob *= 0.319
        else:
            prob *= 0.044
prob

# Generate the k-mer Composition of a String
data = read_data("ba3a (1)")
k = int(data[0])
seq = data[1]
list = []
for i in range(len(seq)-k+1):
    if seq[i:i+k] not in list:
        list.append(seq[i:i+k])
for i in list:
    print(i)

# Find All Occurrences of a Pattern in a String
data = read_data("ba1d")
pattern = data[0]
seq = data[1]
positions = []
for i in range(len(seq)-len(pattern)+1):
    if seq[i:i+len(pattern)] == pattern:
        positions.append(i)
for i in positions:
    print(i,end=" ")

# Find a Position in a Genome Minimizing the Skew
data = read_data("ba1f")
seq = data[0]
skew = [0]
C = 0
G = 0
for i in seq:
    if i == "C":
        C += 1
    elif i == "G":
        G += 1
    skew.append(G-C)
for i, x in enumerate(skew):
    if x == min(skew):
        print(i,end=" ")

# Find All Occurrences of a Collection of Patterns in a String
data = read_data("ba9n")
seq = data[0]
patterns = data[1:]
positions = []
for pattern in patterns:
    for i in range(len(seq)-len(pattern)+1):
        if seq[i:i+len(pattern)] == pattern:
            positions.append(i)
for i in positions:
    print(i,end=" ")

# Find All Approximate Occurrences of a Pattern in a String
data = read_data("ba1h")
pattern = data[0]
seq = data[1]
d = int(data[2])
result = []
for i in range(len(seq)-len(pattern)+1):
    mis_counter = 0
    for j in range(len(pattern)):
        if seq[i+j] != pattern[j]:
            mis_counter += 1
    if mis_counter <= d:
        result.append(i)
for i in result:
    print(i,end=" ")

# Find All Approximate Occurrences of a Collection of Patterns in a String
data = read_data_list("ba9o")
seq = data[0][0]
patterns = data[1]
d = int(data[2][0])
result = []
for pattern in patterns:
    for i in range(len(seq)-len(pattern)+1):
        mis_counter = 0
        for j in range(len(pattern)):
            if seq[i+j] != pattern[j]:
                mis_counter += 1
                if mis_counter > d:
                    break
        if mis_counter <= d:
            result.append(i)
result.sort()
for i in result:
    print(i,end=" ")

# Translate an RNA String into an Amino Acid String
from Bio.Seq import translate
data = read_data("ba4a")
seq = data[0]
translate(seq)

# Find the Minimum Number of Coins Needed to Make Change: UNSOLVED
data = read_data_list("test")
money = [int(data[0][0])]
str_coins = data[1][0].split(",")
coins = []
for coin in str_coins:
    coins.append(int(coin))
coins.reverse()

counter = 0
for coin in coins:
    if money[-1] // coin != 0:
        counter += money[-1]//coin
        money.append(money[-1] % coin)

# Find Frequent Words with Mismatches and Reverse Complements
from Bio.Seq import reverse_complement
data = read_data_list("ba1j (1)")
seq = data[0][0]
reverse = reverse_complement(seq)
list = [seq,reverse]
k = int(data[1][0])
d = int(data[1][1])

def gen_all_k_mer(n):
    list_result = []
    list = [i for i in product(["A","T","C","G"], repeat=n)]
    str_list = []
    for tuple in list:
        string = ""
        for i in tuple:
            string += i
        str_list.append(string)
    return sorted(str_list)
kmers = {}
kmers_list = gen_all_k_mer(k)
for kmer in kmers_list:
    counter = 0
    for dna in list:
        for i in range(len(dna)-len(kmer)+1):
            if dna[i:i+len(kmer)] == kmer:
                counter += 1
    kmers[kmer] = counter
result = {}
for i in kmers:
    result[i] = kmers[i]
for kmer in kmers:
    for xkmer in kmers:
        if kmer != xkmer:
            mis_counter = 0
            for i in range(len(kmer)):
                if kmer[i] != xkmer[i]:
                    mis_counter += 1
                    if mis_counter > d:
                        break
            if mis_counter <= d:
                result[kmer] += kmers[xkmer]
keys = [i for i in result.keys()]
values = [i for i in result.values()]
for i, x in enumerate(values):
    if x == max(values):
        print(keys[i],end=" ")

# Find a Longest Common Subsequence of Two Strings
def lcs(s1, s2):
    matrix = [["" for x in range(len(s2))] for x in range(len(s1))]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                if i == 0 or j == 0:
                    matrix[i][j] = s1[i]
                else:
                    matrix[i][j] = matrix[i-1][j-1] + s1[i]
            else:
                matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1], key=len)
    cs = matrix[-1][-1]
    return len(cs), cs
data = read_data("ba5c")
s1 = data[0]
s2 = data[1]
lcs(s1,s2)

# Find the Longest Repeat in a String
def more_than_once(seq, str):
    counter = 0
    for i in range(len(seq)-len(str)+1):
        if seq[i:i+len(str)] == str:
            counter += 1
    if counter >= 2:
        return True
    else:
        return False
data = read_data("ba9d")
seq = data[0]
substr = []
for i in range(len(seq)):
    for j in range(len(seq)-i):
        if more_than_once(seq, seq[i:i+j+1]) and seq[i:i+j+1] not in substr:
            substr.append(seq[i:i+j+1])
length = []
for i in substr:
    length.append(len(i))
for i,x in enumerate(length):
    if x == max(length):
        print(substr[i])

# Generate the Theoretical Spectrum of a Linear Peptide
data = read_data("ba4j")
seq = data[0]
protein_mass = read_data_list("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[amino[0]] = int(eval(amino[1]))
sub = []
for i in range(len(seq)):
    for j in range(len(seq)-i):
        sub.append(seq[i:i+j+1])
mass_list = [0]
for subseq in sub:
    mass = 0
    for amino in subseq:
        mass += mass_dict[amino]
    mass_list.append(mass)
mass_list.sort()
for i in mass_list:
    print(i,end=" ")

# Find the Most Frequent Words with Mismatches in a String: SLOW!
data = read_data_list("ba1i")
seq = data[0][0]
k = int(data[1][0])
d = int(data[1][1])
def gen_all_k_mer(n):
    list_result = []
    list = [i for i in product(["A","T","C","G"], repeat=n)]
    str_list = []
    for tuple in list:
        string = ""
        for i in tuple:
            string += i
        str_list.append(string)
    return sorted(str_list)
list = [seq]
kmers = {}
kmers_list = gen_all_k_mer(k)
for kmer in kmers_list:
    counter = 0
    for dna in list:
        for i in range(len(dna)-len(kmer)+1):
            if dna[i:i+len(kmer)] == kmer:
                counter += 1
    kmers[kmer] = counter
result = {}
for i in kmers:
    result[i] = kmers[i]
for kmer in kmers:
    for xkmer in kmers:
        if kmer != xkmer:
            mis_counter = 0
            for i in range(len(kmer)):
                if kmer[i] != xkmer[i]:
                    mis_counter += 1
                    if mis_counter > d:
                        break
            if mis_counter <= d:
                result[kmer] += kmers[xkmer]
keys = [i for i in result.keys()]
values = [i for i in result.values()]
for i, x in enumerate(values):
    if x == max(values):
        print(keys[i],end=" ")

# Convert a Peptide into a Peptide Vector
data = read_data("ba11c")
seq = data[0]
protein_mass = read_data_list("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[amino[0]] = int(eval(amino[1]))
sub = []
for i in range(len(seq)):
    sub.append(seq[:i+1])
mass_list = []
for subseq in sub:
    mass = 0
    for amino in subseq:
        mass += mass_dict[amino]
    mass_list.append(mass)
mass_list.sort()
result_list = []
for i in range(1,max(mass_list)+1):
    if i in mass_list:
        result_list.append(1)
    else:
        result_list.append(0)
for i in result_list:
    print(i, end=" ")

# Construct the Overlap Graph of a Collection of k-mers
data = read_data("ba3c")
for seq in data:
    for xseq in data:
        if seq != xseq:
            if seq[1:] == xseq[:-1]:
                print(seq, " -> ",xseq)

# Construct the Graph of a Spectrum
data = read_data_list("ba11a")
str_spectrum = data[0]
spectrum = [0]
for i in str_spectrum:
    spectrum.append(int(i))
protein_mass = read_data_list("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[int(eval(amino[1]))] = amino[0]
for i in range(len(spectrum)-1):
    for j in range(1,len(spectrum)-i):
        if spectrum[i+j] - spectrum[i] < 200:
            if spectrum[i+j] - spectrum[i] in mass_dict and spectrum[i+j] - spectrum[i] > 5:
                print(str(spectrum[i])+"->"+str(spectrum[i+j])+":"+
                      str(mass_dict[spectrum[j+i]-spectrum[i]]))

# Construct the Suffix Array of a String
data = read_data("ba9g")
seq = data[0]
suffixes = []
for i in range(len(seq)):
    suffixes.append(seq[i:]+str(i))
suffixes.sort()
for suffix in suffixes:
    for i in range(len(suffix)):
        if suffix[i] == "$":
            print(suffix[i+1:]+",",end=" ")

# Convert a Peptide Vector into a Peptide
data = read_data_list("ba11d")
str_vector = data[0]
vector = []
for str in str_vector:
    vector.append(int(str))
mass_list = [0]
for i in range(len(vector)):
    if vector[i] == 1:
        mass_list.append(i+1)
protein_mass = read_data_list("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[int(eval(amino[1]))] = amino[0]
protein = ""
for i in range(len(mass_list)-1):
    if mass_list[i+1] - mass_list[i] in mass_dict:
        protein += mass_dict[mass_list[i+1] - mass_list[i]]
print(protein)

# Implement TrieMatching
data = read_data("ba9b")
seq = data[0]
patterns = data[1:]
positions = []
for pattern in patterns:
    for i in range(len(seq)-len(pattern)+1):
        if seq[i:i+len(pattern)] == pattern:
            positions.append(i)
for i in positions:
    print(i,end=" ")

# Compute the Number of Breakpoints in a Permutation: UNSOLVED
data = read_data_list("test")
r_permutation = data[0]
permutation = []
for i in r_permutation:
    if "(" in i:
        permutation.append(eval(i.replace("(","")))
    elif ")" in i:
        permutation.append(eval(i.replace(")","")))
    else:
        permutation.append(eval(i))
breaks = 0
for i in range(1,len(permutation)-1):
    if permutation[i-1] - permutation[i] != permutation[i] - permutation[i+1]:
        breaks += 1
breaks

# Pattern Matching with the Suffix Array
data = read_data("ba9h")
seq = data[0]
patterns = data[1:]
positions = []
for pattern in patterns:
    for i in range(len(seq)-len(pattern)+1):
        if seq[i:i+len(pattern)] == pattern:
            positions.append(i)
for i in positions:
    print(i,end=" ")

# Find the Shortest Non-Shared Substring of Two Strings
data = read_data("ba9f")
seq1 = data[0]
seq2 = data[1]
substrings = []
for i in range(1,len(seq1)):
    for j in range(len(seq1)-i+1):
        count = 0
        for m in range(len(seq2)-i+1):
            if seq2[m:m+i] == seq1[j:j+i]:
                count += 1
                if count != 0:
                    break
        if count == 0:
            substrings.append(seq1[j:j+i])
    if substrings:
        break
print(substrings[0])

# Find the Longest Substring Shared by Two Strings: UNSOLVED
data = read_data("test")
seq1 = data[0]
seq2 = data[1]
result = []
for i in range(1,len(seq1)):
    for j in range(len(seq1)-i+1):
        for m in range(len(seq2)-i+1):
            if seq2[m:m+i] == seq1[j:j+i]:
                result.append(seq1[j:j+i])
                break
        if len(result[-1]) == i:
            break
print(result[-1])

# Generate the d-Neighborhood of a String
data = read_data("ba1n")
target = data[0]
d = int(data[1])
kmers = gen_all_k_mer(len(target))
result = []
for kmer in kmers:
    mis_counter = 0
    for i in range(len(kmer)):
        if kmer[i] != target[i]:
            mis_counter += 1
            if mis_counter > d:
                break
    if mis_counter <= d:
        result.append(kmer)
for i in result:
    print(i)

# 



