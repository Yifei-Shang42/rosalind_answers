# Introduction to the Bioinformatics Armory
from Bio.Seq import Seq
def read_data(filename):
    infile = open("C:\\Users\\surface\\Downloads\\rosalind_" + filename + ".txt", "r")
    lines = infile.readlines()
    data = []
    for line in lines:
        data.append(line.split()[0])
    infile.close()
    return data
my_seq = read_data("ini")[0]
my_seq = Seq(my_seq)
my_seq.count("A")
my_seq.count("C")
my_seq.count("G")
my_seq.count("T")

# Introduction to Protein Databases: BAD CODE
from Bio import ExPASy
from Bio import SwissProt
handle = ExPASy.get_sprot_raw('B5ZC00')
record = SwissProt.read(handle)

# GenBank Introduction
from Bio import Entrez
Entrez.email = "yifeis.20@intl.zju.edu.cn"
handle = Entrez.esearch(db="nucleotide",
                        term='"Soriculus"[Organism] AND ("2000/03/31"[Publication Date] : "2010/09/29"[Publication Date])')
record = Entrez.read(handle)
record["Count"]

# Data Formats
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "yifeis.20@intl.zju.edu.cn"
handle = Entrez.efetch(db="nucleotide",
                       id=["JX469983"],
                       rettype="fasta")
result_list = []
for record in SeqIO.parse(handle, "fasta"):
    result_list.append(len(record.seq))
    print(record.id,end=" ")
result_list

# FASTQ format introduction
    # convert between formats
"http://sequenceconversion.bugaco.com/converter/biology/sequences/index.php"
from Bio import SeqIO
with open("cor6_6.gb", "rU") as input_handle:
    with open("cor6_6.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")

# Pairwise Global Alignment
"https://www.ebi.ac.uk/Tools/psa/emboss_needle/"
from Bio import Entrez
from Bio import SeqIO, pairwise2
Entrez.email = "yifeis.20@intl.zju.edu.cn"
genbank_ids = "JX398977.1, JQ342169.1"
handle = Entrez.efetch(db="nucleotide", id=[genbank_ids], rettype="fasta")
records = [i.seq for i in SeqIO.parse(handle, "fasta")]
print(pairwise2.align.globalms(records[0], records[1], 5, -4, -10, -1)[0][2])

# Protein Translation
from Bio.Seq import translate
data = read_data("ptra")
dna = data[0]
protein = data[1]
protein_dict = {}
for i in [1,2,3,4,5,6,9,10,11,12,13,14,15]:
    protein_dict[i] = translate(dna, table=i)
for i in protein_dict:
    if protein == protein_dict[i][:-1]:
        print(i, end=" ")

# Complementing a Strand of DNA
from Bio.Seq import Seq
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
data = read_fasta("rvco")
counter = 0
for seq in data:
    my_seq = Seq(seq)
    if my_seq == my_seq.reverse_complement():
        counter += 1
print(counter)

# Finding Genes with ORFs
"http://www.bioinformatics.org/sms2/orf_find.html"

# Read Quality Distribution
from Bio import SeqIO
from statistics import *
record = SeqIO.parse("C:\\Users\\surface\\Downloads\\rosalind_phre.txt",
            "fastq")
list = []
for i in record:
    list.append(mean(i.letter_annotations['phred_quality']))
counter = 0
for mean in list:
    if mean < 22:
        counter += 1
print(counter)

# Read Filtration by Quality
"https://usegalaxy.org/"
    # filter by quality

# Base Quality Distribution
from Bio import SeqIO
from statistics import *
def bphr(data):
    count = 0
    with open(data, "r") as f:
        threshold = int(f.readline())
        qualities = []
        for record in SeqIO.parse(f, "fastq"):
            quality = record.letter_annotations["phred_quality"]
            qualities.append(quality)
    for i in range(len(qualities[0])):
        if sum([q[i] for q in qualities])/len(qualities) < threshold:
            count += 1
    return count
bphr("C:\\Users\\surface\\Downloads\\rosalind_bphr.txt")

# Base Filtration by Quality
"https://usegalaxy.org/root?tool_id=fastq_quality_trimmer"




