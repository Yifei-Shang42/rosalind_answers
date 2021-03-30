# Counting DNA Nucleotides
infile = open("C:\\Users\\surface\\Downloads\\rosalind_dna.txt", "r")
lines = infile.readlines()
for line in lines:
    bases = line.split()
infile.close()
G = 0
T = 0
A = 0
C = 0
for base in bases[0]:
    if base == "A":
        A += 1
    elif base == "G":
        G += 1
    elif base == "T":
        T += 1
    else:
        C += 1
print(A, T, G, C)

# Transcribing DNA into RNA
infile = open("C:\\Users\\surface\\Downloads\\rosalind_rna.txt", "r")
lines = infile.readlines()
for line in lines:
    data = line.split()
infile.close()

output = []
for base in data[0]:
    if base == "A":
        output.append("A")
    elif base == "T":
        output.append("U")
    elif base == "C":
        output.append("C")
    else:
        output.append("G")
string = ""
for i in output:
    string = string + i
print(string)

# Complementing a Strand of DNA
infile = open("C:\\Users\\surface\\Downloads\\rosalind_revc.txt", "r")
lines = infile.readlines()
for line in lines:
    data = line.split()
infile.close()
complement = []
    # generate complement strand
for base in data[0]:
    if base == "A":
        complement.append("T")
    elif base == "T":
        complement.append("A")
    elif base == "C":
        complement.append("G")
    else:
        complement.append("C")
    # generate reverse strand
string = ""
reverse = reversed(complement)
for i in reverse:
    string = string + i
print(string)

# Mendel's First Law
infile = open("C:\\Users\\surface\\Downloads\\rosalind_iprb.txt", "r")
lines = infile.readlines()
for line in lines:
    data = line.split()
infile.close()

k = int(data[0])
m = int(data[1])
n = int(data[2])
N = k + n + m
ratio = (k*(k-1) + 2*k*(n+m) + 0.75*m*(m-1) + m*n)/(N*(N-1))
print(ratio)

# Translating RNA into Protein
    # create codon dict
infile = open("C:\\Users\\surface\\Downloads\\codon.txt", "r")
codons = infile.read()
codon_split = codons.split()
infile.close()
codon_base = []
codon_amino = []
for i in range(len(codon_split)):
    if len(codon_split[i]) > 1 and codon_split[i] != "Stop":
        codon_base.append(codon_split[i])
    else:
        codon_amino.append(codon_split[i])
codon_dict = dict(zip(codon_base, codon_amino))

infile = open("C:\\Users\\surface\\Downloads\\rosalind_prot.txt", "r")
lines = infile.readlines()
for line in lines:
    data = line.split()
infile.close()
protein = []
for i in range(0, len(data[0]), 3):
    if codon_dict[data[0][i:i+3]] == "Stop":
        break
    else:
        protein.append(codon_dict[data[0][i:i+3]])

string = ""
for i in protein:
    string = string + i
print(string)

# Rabbits and Recurrence Relations
infile = open("C:\\Users\\surface\\Downloads\\rosalind_fib.txt", "r")
lines = infile.readlines()
for line in lines:
    data = line.split()
infile.close()
k = int(data[1])
n = int(data[0])
count = [1, 1]
for i in range(n-2):
    count.append(count[i+1] + count[i]*k)
print(count[-1])

# Computing GC Content
def read_data(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split()[0])
    infile.close()
    return data
data = read_data("gc")
    # create dict of bases
def read_fasta_as_dict(data):
    dict_data = {}
    for i in range(len(data)):
        if data[i][0] == ">":
            base = ""
            for j in range(len(data[i:-1])):
                if data[i+j+1][0] != ">":
                    base += data[i+j+1]
                else:
                    break
            dict_data[data[i]] = base
    return dict_data
dict_data = read_fasta_as_dict(data)
dict_result = {}
for key in list(dict_data.keys()):
    counter = 0
    for base in dict_data[key]:
        if base == "C" or base == "G":
            counter += 1
    dict_result[key] = 100*counter/len(dict_data[key])
print(dict_result)

# Counting Point Mutations
infile = open("C:\\Users\\surface\\Downloads\\rosalind_ba1g.txt", "r")
lines = infile.readlines()
data = []
for line in lines:
    data.append(line.split()[0])
infile.close()
gene_1 = data[0]
gene_2 = data[1]
hamming_distance = 0
for i in range(len(gene_1)):
    if gene_1[i] != gene_2[i]:
        hamming_distance += 1
print(hamming_distance)

# Finding a Motif in DNA
infile = open("C:\\Users\\surface\\Downloads\\rosalind_subs.txt", "r")
lines = infile.readlines()
data = []
for line in lines:
    data.append(line.split()[0])
infile.close()

s = data[0]
t = data[1]

for i in range(len(s) - len(t) - 1):
    if s[i:i+len(t)] == t:
        print(i+1)

# Consensus and Profile
def read_data(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split()[0])
    infile.close()
    return data
data = read_data("cons")
def read_fasta_as_dict(data):
    dict_data = {}
    for i in range(len(data)):
        if data[i][0] == ">":
            base = ""
            for j in range(len(data[i:-1])):
                if data[i+j+1][0] != ">":
                    base += data[i+j+1]
                else:
                    break
            dict_data[data[i]] = base
    return dict_data
dict_data = read_fasta_as_dict(data)
list_data = list(dict_data.values())
A_list = []
T_list = []
C_list = []
G_list = []
for i in range(len(list_data[0])):
    A = 0
    C = 0
    T = 0
    G = 0
    for j in range(len(list_data)):
        if list_data[j][i] == "A":
            A += 1
        elif list_data[j][i] == "T":
            T += 1
        elif list_data[j][i] == "C":
            C += 1
        else:
            G += 1
    A_list.append(A)
    G_list.append(G)
    C_list.append(C)
    T_list.append(T)
print("A: ", end=" ")
for i in A_list:
    print(i, end=" ")
print("C: ", end=" ")
for i in C_list:
    print(i, end=" ")
print("G: ", end=" ")
for i in G_list:
    print(i, end=" ")
print("T: ", end=" ")
for i in T_list:
    print(i, end=" ")
consensus = ""
for i in range(len(list_data[0])):
    if max(T_list[i], G_list[i], A_list[i], C_list[i]) == T_list[i]:
        consensus += "T"
    elif max(T_list[i], G_list[i], A_list[i], C_list[i]) == G_list[i]:
        consensus += "G"
    elif max(T_list[i], G_list[i], A_list[i], C_list[i]) == A_list[i]:
        consensus += "A"
    else:
        consensus += "C"
print(consensus)

# Mortal Fibonacci Rabbits: UNSOLVED!!!!!
#run for n months, rabbits die after m months.
n, m = input("Enter months to run, and how many months rabbits live, separated by a space ").split()
n, m = int(n), int(m)
generations = [1, 1] #Seed the sequence with the 1 pair, then in their reproductive month.
def fib(i, j):
    count = 2
    while (count < i):
        if (count < j):
            generations.append(generations[-2] + generations[-1]) #recurrence relation before rabbits start dying (simply fib seq Fn = Fn-2 + Fn-1)
        elif (count == j or count == j+1):
            print ("in base cases for newborns (1st+2nd gen. deaths)") #Base cases for subtracting rabbit deaths (1 death in first 2 death gens)
            generations.append((generations[-2] + generations[-1]) - 1)#Fn = Fn-2 + Fn-1 - 1
        else:
            generations.append((generations[-2] + generations[-1]) - (generations[-(j+1)])) #Our recurrence relation here is Fn-2 + Fn-1 - Fn-(j+1)
        count += 1
    return (generations[-1])


print (fib(n, m))
print ("Here's how the total population looks by generation: \n" + str(generations))

# Calculating Expected Offspring
data = read_data_list("iev")
a = int(data[0][0])
b = int(data[0][1])
c = int(data[0][2])
d = int(data[0][3])
e = int(data[0][4])
f = int(data[0][5])
total = 2*a + 2*b + 2*c + 2*0.75*d + 2*0.5*e
print(total)

# Inferring mRNA from Protein
from numpy import *
data = read_data("mrna")
protein = str(data[0][0])
def gen_codon_dict():
    infile = open("C:\\Users\\surface\\Downloads\\codon.txt", "r")
    codons = infile.read()
    codon_split = codons.split()
    infile.close()
    codon_base = []
    codon_amino = []
    for i in range(len(codon_split)):
        if len(codon_split[i]) > 1 and codon_split[i] != "Stop":
            codon_base.append(codon_split[i])
        else:
            codon_amino.append(codon_split[i])
    codon_dict = dict(zip(codon_base, codon_amino))
    return codon_dict
codon_dict = gen_codon_dict()
    # for each amino:
        # count all keys with value == amino
        # put count into a list
counter_list = []
for amino in data[0][0]:
    counter = 0
    for i in codon_dict:
        if codon_dict[i] == amino:
            counter += 1
    counter_list.append(counter)
result = 1
for i in counter_list:
    result = result*i
print((result*3)%1000000)

# Overlap Graphs
data = read_data_list("grph")
dict_data = {}
for i in range(len(data)):
    if data[i][0][0] == ">":
        base = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0][0] != ">":
                base += data[i+j+1][0]
            else:
                break
        dict_data[data[i][0][1:]] = base
for string in dict_data:
    for gene in dict_data:
        if gene != string and dict_data[string][-3:] == dict_data[gene][0:3]:
            print(string + " " + gene)

# Enumerating Gene Orders
from itertools import permutations
from math import *
data = read_data_list("perm")
n = int(data[0][0])
print(factorial(n))
list_result = list(permutations(range(1, n+1), n))
for i in range(len(list_result)):
    for j in range(len(list_result[0])):
        print(str(list_result[i][j]),end= " ")
    print("")

# Independent Alleles
from math import *
data = read_data("lia")
k = int(data[0][0])
N = int(data[0][1])
total_pseudo = 0
for i in range(N):
    total_pseudo += comb(2**k, i)*((0.25)**i)*((0.75)**(2**k-i))
print(1-total_pseudo)

# Locating Restriction Sites
data = read_data_list("revp")
for i in range(len(data)):
    if data[i][0][0] == ">":
        base = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0][0] != ">":
                base += data[i+j+1][0]
            else:
                break

def reverse_pali_location_length(data):
    def gen_reverse_com(string):
        reverse_string = string[::-1]
        reverse_compliment = ""
        for i in reverse_string:
            if i == "A":
                reverse_compliment += "T"
            elif i == "T":
                reverse_compliment += "A"
            elif i == "C":
                reverse_compliment += "G"
            else:
                reverse_compliment += "C"
        return reverse_compliment
    for i in range(0, len(data)-3):
        for j in [4,6,8,10,12]:
            reverse = gen_reverse_com(data[i:i+j])
            if reverse == data[i:i+j] and i+j-1 < len(base):
                print(i+1, j)
reverse_pali_location_length(base)

# Calculating Protein Mass
protein_mass = read_data("protein_mass")
mass_dict = {}
for amino in protein_mass:
    for i in amino:
        mass_dict[amino[0]] = eval(amino[1])
data = read_data("prtm")
protein_sequence = data[0][0]
mass_total = 0
for amino in protein_sequence:
    mass_total += mass_dict[amino]
print(mass_total)

# Open Reading Frames: NOT SOLVED
f = open("\\Users\\surface\\Downloads\\rosalind_xxxx.txt", "r")
lines = f.readlines()
data = []
for line in lines:
    data.append(line)
for line in data:
    if len(line) > 1:
        if line[-2] == "*":
            print(line[:-2])

data = read_data_list("orf")

for i in range(len(data)):
    if data[i][0][0] == ">":
        bases = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0][0] != ">":
                bases += data[i+j+1][0]
            else:
                break
codon_dict = gen_codon_dict()

def gen_reverse_com(bases):
    reverse_string = bases[::-1]
    reverse_compliment = ""
    for i in range(len(reverse_string)):
        if reverse_string[i] == "A":
            reverse_compliment += "T"
        elif reverse_string[i] == "T":
            reverse_compliment += "A"
        elif reverse_string[i] == "C":
            reverse_compliment += "G"
        else:
            reverse_compliment += "C"
    return reverse_compliment
reverse = gen_reverse_com(bases)

def transcribe_template(data):
    rna = ""
    for base in data:
        if base == "A":
            rna += "A"
        elif base == "T":
            rna += "U"
        elif base == "C":
            rna += "C"
        else:
            rna += "G"
    return rna

rna_bases = transcribe_template(bases)
rna_reverse = transcribe_template(reverse)

def translate_many(mrna):
    protein_list = []
    for i in range(len(mrna)-2):
        protein = ""
        if mrna[i:i+3] == "AUG":
            for j in range(0, len(mrna)-2-i, 3):
                if codon_dict[mrna[i+j:i+j+3]] == "Stop":
                    break
                else:
                    protein += codon_dict[mrna[i+j:i+j+3]]
            protein_list.append(protein)
    return protein_list

translate_many(rna_bases)
translate_many(rna_reverse)

# Enumerating k-mers Lexicographically
import itertools
data = read_data_list("lexf")
n = int(data[1][0])
alphabets = data[0]
keywords = [''.join(i) for
            i in itertools.product(alphabets, repeat = 3)]
for i in keywords:
    print(i)

# Finding a Shared Motif
data = read_data("lcsm")
genes = []
for i in range(len(data)):
    if data[i][0] == ">":
        bases = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0] != ">":
                bases += data[i+j+1]
            else:
                break
    genes.append(bases)

def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr

def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True

print(long_substr(genes))

# RNA Splicing
data = read_data("splc")
data_edit = []
for i in range(len(data)):
    if data[i][0] == ">":
        bases = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0] != ">":
                bases += data[i+j+1]
            else:
                break
        data_edit.append(bases)

sequence = data_edit[0]
introns = data_edit[1:]

for i in range(len(introns)):
    for j in range(len(sequence)-len(introns[i])):
        if sequence[j:j+len(introns[i])] == introns[i]:
            sequence = sequence.replace(introns[i], "")

new_gene = sequence

def transcribe_template(data):
    rna = ""
    for base in data:
        if base == "A":
            rna += "A"
        elif base == "T":
            rna += "U"
        elif base == "C":
            rna += "C"
        else:
            rna += "G"
    return rna
exons = transcribe_template(new_gene)

def translate(mrna):
    protein = ""
    for i in range(len(mrna)-2):
        if mrna[i:i+3] == "AUG":
            for j in range(0, len(mrna)-2-i, 3):
                if codon_dict[mrna[i+j:i+j+3]] == "Stop":
                    break
                else:
                    protein += codon_dict[mrna[i+j:i+j+3]]
            break
    return protein

translate(exons)

# Transitions and Transversions
data = read_data("tran")
data_edit = []
for i in range(len(data)):
    if data[i][0] == ">":
        bases = ""
        for j in range(len(data[i:-1])):
            if data[i+j+1][0] != ">":
                bases += data[i+j+1]
            else:
                break
        data_edit.append(bases)

s1 = data_edit[0]
s2 = data_edit[1]

transitions = 0
transversions = 0
for i in range(len(s1)):
    if s1[i] != s2[i]:
        if s1[i] == "A":
            if s2[i] == "G":
                transitions += 1
            else:
                transversions += 1
        elif s1[i] == "T":
            if s2[i] == "C":
                transitions += 1
            else:
                transversions += 1
        elif s1[i] == "G":
            if s2[i] == "A":
                transitions += 1
            else:
                transversions += 1
        else:
            if s2[i] == "T":
                transitions += 1
            else:
                transversions += 1

print(transitions/transversions)

# Introduction to Random Strings
from math import *
data = read_data_list("prob")
sequence = data[0][0]
list_gc_content = data[1]

prob_list = []
for str_gc_content in list_gc_content:
    gc_content = eval(str_gc_content)
    probability = 1
    for base in sequence:
        if base == "A" or base == "T":
            probability *= (1-gc_content)/2
        if base == "C" or base == "G":
            probability *= gc_content/2
    prob_list.append(probability)
for prob in prob_list:
    print(log(prob, 10), end=" ")

# Finding a Protein Motif
from urllib.request import urlopen
data = read_data("mprt")
protein_dict = {}
for protein_id in data:
    url = "https://www.uniprot.org/uniprot/" + protein_id + ".fasta"
    page = urlopen(url)
    html_bytes = page.read()
    html = html_bytes.decode("utf-8")
    for i in range(len(html)):
        if html[i] == "\n":
            protein = html[i:]
            break
    protein = protein.replace("\n", "")
    protein_dict[protein_id] = protein

location_dict = {}
for id in protein_dict:
    location = []
    for i in range(len(protein_dict[id])-3):
        if protein_dict[id][i] == "N":
            if protein_dict[id][i+1] != "P":
                if protein_dict[id][i+2] == "S" or protein_dict[id][i+2] == "T":
                    if protein_dict[id][i+3] != "P":
                        location.append(i+1)
    location_dict[id] = location

for id in location_dict:
    if location_dict[id] != []:
        print(id)
        for location in location_dict[id]:
            print(location, end=" ")
    print("")

# Perfect Matchings and RNA Secondary Structures
from math import *
data = read_data("pmch")
mrna = ""
for i in data[1:]:
    mrna += i
A = 0
G = 0
for base in mrna:
    if base == "A":
        A += 1
    elif base == "G":
        G += 1
factorial(A)*factorial(G)

# Partial Permutations
from math import *
data = read_data_list("pper")
n = int(data[0][0])
k = int(data[0][1])
result = (factorial(n)/factorial(n-k)) % 1000000
result

# Finding a Spliced Motif
def read_fasta(filename):
    data = read_data(filename)
    data_edit = []
    for i in range(len(data)):
        if data[i][0] == ">":
            bases = ""
            for j in range(len(data[i:-1])):
                if data[i+j+1][0] != ">":
                    bases += data[i+j+1]
                else:
                    break
            data_edit.append(bases)
    return data_edit

data = read_fasta("sseq")
s = data[0]
t = data[1]

list_result = [0]
for base in t:
    for i in range(len(s)):
        if s[i] == base and i+1 not in list_result and i > list_result[-1]:
            list_result.append(i+1)
            break
list_result.remove(0)
for location in list_result:
    print(location, end=" ")

# Enumerating Oriented Gene Orderings
from itertools import permutations
from math import *
data = read_data("sign")
n = int(data[0])
print(factorial(n)*2**n)
list_result = []
for i in permutations(range(1, n+1), n):
    list_result.append(i)

alphabets = [1,-1]
keywords = [i for i in product(alphabets, repeat = n)]

result = []
for i in range(len(list_result)):
    for m in range(len(keywords)):
        list = []
        for j in range(len(list_result[i])):
            list.append(list_result[i][j]*keywords[m][j])
        result.append(list)
for i in range(len(result)):
    for j in range(len(result[0])):
        print(str(result[i][j]),end= " ")
    print("")

# Completing a Tree
data = read_data_list("tree")
n = int(data[0][0])
l = data[1:]

import networkx
from networkx.algorithms.components.connected import connected_components


def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """
        treat `l` as a Graph and returns it's edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current

G = to_graph(l)
a = connected_components(G)
result = []
for i in a:
    result.append(i)

list_len = []
for list in l:
    for number in list:
        if number not in list_len:
            list_len.append(number)
print(n-len(list_len)+len(result)-1)

# Counting Phylogenetic Ancestors
data = read_data("inod")
n = int(data[0])
N = 0
while n > 1:
    if n % 2 == 1:
        N += (n - 1) / 2
        n = (n + 1) / 2
    else:
        N += n/2
        n = n/2
print(N-1)

# Ordering Strings of Varying Length Lexicographically: NOT SOLVED
from itertools import permutations
def read_data_list(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split())
    infile.close()
    return data
data = read_data_list("lexv")
n = int(data[1][0])
alphabet = data[0]
permutation_list = []
for i in range(n):
    permutation_list.append([i for i in product(alphabet, repeat=i+1)])

string_list = []
for length in permutation_list:
    sub_list = []
    for tuple in length:
        string = ""
        for letter in tuple:
            string += letter
        sub_list.append(string)
    string_list.append(sub_list)


alphabet_dict = {}
for i in range(len(alphabet)):
    alphabet_dict[alphabet[i]] = i+1

# Counting Subsets
2**959 % 1000000

# Introduction to Set Operations: UNSOLVED
data = read_data_list("test")
n = int(data[0][0])
list_1 = []
list_2 = []
list_1.append(data[1][0][1:-1])
list_2.append(data[2][0][1:-1])
for i in data[1][1:]:
    list_1.append(i[:-1])
for j in data[2][1:]:
    list_2.append(j[:-1])
    # union
union = []
for i in list_1:
    union.append(i)
for number in list_2:
    if number not in list_1:
        union.append(number)
print("{",end="")
for i in range(len(union)-1):
    print(union[i],end=", ")
print(union[-1],end="")
print("}")
    # intersect
intersect = []
if len(list_1) >= len(list_2):
    for i in list_1:
        if i in list_2:
            intersect.append(i)
else:
    for i in list_2:
        if i in list_1:
            intersect.append(i)
print("{",end="")
for i in range(len(intersect)-1):
    print(intersect[i],end=", ")
print(intersect[-1],end="")
print("}")
    # A-B
copy_1 = []
for i in list_1:
    copy_1.append(i)
for i in copy_1:
    if i in list_2:
        copy_1.remove(i)
print("{",end="")
for i in range(len(copy_1)-1):
    print(copy_1[i],end=", ")
print(copy_1[-1],end="")
print("}")
    # B-A
copy_2 = []
for i in list_2:
    copy_2.append(i)
for i in copy_2:
    if i in list_1:
        copy_2.remove(i)
print("{",end="")
for i in range(len(copy_2)-1):
    print(copy_2[i],end=", ")
print(copy_2[-1],end="")
print("}")
    # complement
total_list = [str(i) for i in range(1, n+1)]
for i in list_1:
    if i in total_list:
        total_list.remove(i)
print("{",end="")
for i in range(len(total_list)-1):
    print(total_list[i],end=", ")
print(total_list[-1],end="")
print("}")

total_list = [str(i) for i in range(1, n+1)]
for i in list_2:
    if i in total_list:
        total_list.remove(i)
print("{",end="")
for i in range(len(total_list)-1):
    print(total_list[i],end=", ")
print(total_list[-1],end="")
print("}")

# k-Mer Composition
from itertools import *
data = read_data("kmer")
data = read_fasta_as_dict(data)
sequence = ""
for i in data:
    sequence += data[i]

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

kmers = gen_all_k_mer(4)
kmer_dict = {}
for kmer in kmers:
    counter = 0
    for i in range(0, len(sequence)-3):
        if sequence[i:i+4] == kmer:
            counter += 1
    kmer_dict[kmer] = counter

for i in kmer_dict:
    print(kmer_dict[i], end=" ")

# Creating a Distance Matrix
from numpy import *
data = read_data("pdst")
data = read_fasta_as_dict(data)
data_list = []
for i in data:
    data_list.append(data[i])

result_list = []
for i in range(len(data_list)):
    list = []
    for j in range(len(data_list)):
        counter = 0
        for m in range(len(data_list[i])):
            if data_list[i][m] != data_list[j][m]:
                counter += 1
        list.append(counter/len(data_list[j]))
    result_list.append(list)
for row in result_list:
    for column in row:
        print(column, end=" ")
    print("")

# Introduction to Alternative Splicing
from math import *
n = 1830
m = 724
total = 0
for i in range(m,n+1):
    total += comb(n,i)
total % 1000000

# Longest Increasing Subsequence: UNSOLVED
data = read_data_list("test")
n = int(data[0][0])
number_list = []
for number in data[1]:
    number_list.append(eval(number))

total_increasing_list = []
for i in range(n):
    increasing_list = []
    increasing_list.append(number_list[i])
    for j in range(i+1, n):
        if increasing_list[-1] < number_list[j]:
            increasing_list.append(number_list[j])
    total_increasing_list.append(increasing_list)

result1 = max(total_increasing_list)

# Speeding Up Motif Finding: UNSOLVED
data = read_fasta("kmp")
sequence = data[0]
failure_array = [0]
for i in range(1, len(sequence)):
    array = []
    for k in range(1,int(i/2)+1):
        if sequence[i-k+1] == sequence[0]:
            if sequence[:k] == sequence[i-k+1:i+1]:
                array.append(k)
    if array:
        failure_array.append(max(array))
    else:
        failure_array.append(0)
for i in failure_array:
    print(i,end=" ")

# Independent Segregation of Chromosomes
from math import *
n = 41
list = []
for k in range(2*n):
    list.append(comb(2*n,2*n-k)*0.5**(2*n))
result_list = []
for i in range(len(list)):
    result_list.append(sum(list[:i+1]))
result_list.reverse()
for i in result_list:
    print(log10(i),end=" ")

# Counting Disease Carriers
from math import *
data = read_data_list("afrq")
list_prob = []
for i in data[0]:
    list_prob.append(eval(i))
list_result = []
for j in list_prob:
    list_result.append(j+2*sqrt(j)*(1-sqrt(j)))
for i in list_result:
    print(i,end=" ")

# Finding a Shared Spliced Motif
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

data = read_fasta("lcsq")
lcs(data[0],data[1])
product
