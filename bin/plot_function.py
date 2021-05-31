#!/usr/bin/env python
#script to plot miRNAs expression data
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        if ID.count("_") > 1:
            id1 = ID.split("_")
            ID = id1[0]+"_"+id1[1]
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _Alignment_(seqs):
    temp = open("/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/plotTest/temp.txt","w")
    count = 1
    for i in seqs:
        id = str(count)+"_"+str(i[1])
        seq = i[0]
        temp.write(">"+id+"\n"+seq+"\n")
        count += 1
    temp.close()
    os.system("mafft /home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/plotTest/temp.txt > /home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/plotTest/temp_aligned.txt")
    alignment = AlignIO.read(open("/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/plotTest/temp_aligned.txt"), "fasta")
    final = []
    for record in alignment:
        exp = str(record.id).split("_")[1]
        seq = str(record.seq).replace("-",".")
        final.append([seq, exp])
    return final

def _ParseCDHIT_(cdhit):
    clusters = {}
    stars = {}
    expT = {}
    a = open(cdhit,"r")
    for line in a:
        if line.startswith(">"):
            mir = line.strip().replace(">","")
            clusters.setdefault(mir, [])
            stars.setdefault(mir, [])
            expT.setdefault(mir, [0, 0])
        if line.startswith("\t"):
            line1 = line.strip().split(" ")
            seq = line1[1]
            exp = int(line1[2])
            if line1[3] == "canonical":
                clusters[mir].append([seq, exp])
                expT[mir][0] += exp
            if line1[3] == "star":
                stars[mir].append([seq, exp])
                expT[mir][1] += exp
    a.close()
    #print(clusters["id_1003"], stars["id_1003"], expT["id_1003"])
    return clusters, stars, expT

def _PlotS_(AC, AS, id, type, output):

    N = 70

    A = [0]*70
    T = [0]*70
    C = [0]*70
    G = [0]*70

    for i in AC:
        seq = i[0]
        exp = int(i[1])
        for n in range(0, len(seq)):
            if seq[n] == "a":
                A[n] += exp
            if seq[n] == "t":
                T[n] += exp
            if seq[n] == "c":
                C[n] += exp
            if seq[n] == "g":
                G[n] += exp

    for i in AS:
        seq = i[0]
        exp = int(i[1])
        for n in range(0, len(seq)):
            if seq[n] == "a":
                A[n+30] += exp
            if seq[n] == "t":
                T[n+30] += exp
            if seq[n] == "c":
                C[n+30] += exp
            if seq[n] == "g":
                G[n+30] += exp

    r = np.arange(N)    # The position of the bars on the x-axis
    width = 1       # the width of the bars: can also be len(x) sequence

    # Heights of bars1 + bars2
    barsAT = np.add(A, T).tolist()
    barsATC = np.add(barsAT, C).tolist()

    fig = plt.figure()
    # Create brown bars
    plt.bar(r, A, color='green', edgecolor='white', width=width)
    plt.bar(r, T, bottom=A, color='lightgreen', edgecolor='white', width=width)
    plt.bar(r, C, bottom=barsAT, color='blue', edgecolor='white', width=width)
    plt.bar(r, G, bottom=barsATC, color="lightblue", edgecolor='white', width=width)

    # Custom X axis
    plt.ylabel("number of reads")
    plt.xlabel("nucleotide position")
    plt.legend(("A","T","C","G"))
    plt.title(id)

    # Show graphic
    nf = id.split(", ")[0]
    plt.savefig(output+"plots/"+type+"/"+nf+".pdf")

def _PlotC_(AC, id, type, output):

    # textstr = ""
    # for i in AC:
    #     seq = i[0]
    #     exp = str(i[1])
    #     textstr += seq.center(30,".")+exp.rjust(10)+"\n"

    N = 30

    A = [0]*30
    T = [0]*30
    C = [0]*30
    G = [0]*30
    #a - green; t - red; c - blue; g - orange

    for i in AC:
        seq = i[0]
        exp = int(i[1])
        for n in range(0, len(seq)):
            if seq[n] == "a":
                A[n] += exp
            if seq[n] == "t":
                T[n] += exp
            if seq[n] == "c":
                C[n] += exp
            if seq[n] == "g":
                G[n] += exp

    r = np.arange(N)    # The position of the bars on the x-axis
    width = 1       # the width of the bars: can also be len(x) sequence

    # Heights of bars1 + bars2
    barsAT = np.add(A, T).tolist()
    barsATC = np.add(barsAT, C).tolist()

    fig = plt.figure()
    # Create brown bars
    plt.bar(r, A, color='green', edgecolor='white', width=width)
    plt.bar(r, T, bottom=A, color='lightgreen', edgecolor='white', width=width)
    plt.bar(r, C, bottom=barsAT, color='blue', edgecolor='white', width=width)
    plt.bar(r, G, bottom=barsATC, color="lightblue", edgecolor='white', width=width)

    # Custom X axis
    #plt.xticks(r, names, fontweight='bold')
    plt.ylabel("number of reads")
    plt.xlabel("nucleotide position")
    plt.legend(("A","T","C","G"))
    plt.title(id)
    #plt.legend((p1[0], p2[0]), ("A","T","C","G"))
    #plt.gcf().text(0.02, 0.5, "hey ho lets go", fontsize=14)
    #plt.figtext(0.5, 0.01, textstr, ha="center", fontsize=18)
    #fig.text(.1,.1,textstr)

    # Show graphic
    #plt.show()
    nf = id.split(", ")[0]
    plt.savefig(output+"plots/"+type+"/"+nf+".pdf")

# known = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/SB0147_bjussu/known_mirnas.fa"
# K = _ParseFasta_(known)
# novel = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/SB0147_bjussu/novel_mirnas.fa"
# N = _ParseFasta_(novel)
# CD = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/my_pipeline_denovo/SB0147_bjussu/mirna_clusters.txt"
# canonical, stars, expT = _ParseCDHIT_(CD)

def _GenOutput_(known, novel, CD, output):

    if os.path.isdir(output+"plots/") == False:
        os.mkdir(output+"plots/")
    if os.path.isdir(output+"plots/known/") == False:
        os.mkdir(output+"plots/known/")
    if os.path.isdir(output+"plots/novel/") == False:
        os.mkdir(output+"plots/novel/")

    K = _ParseFasta_(known)
    N = _ParseFasta_(novel)
    canonical, stars, expT = _ParseCDHIT_(CD)

    for k in canonical.keys():
        if k in K.keys():
            id = k+", known, "
            AC = _Alignment_(canonical[k])
            if stars[k] != []:
                id += "ExpMature = "+str(expT[k][0])+", ExpStar = "+str(expT[k][1])
                AS = _Alignment_(stars[k])
                _PlotS_(AC, AS, id, "known", output)
            if stars[k] == []:
                id += "ExpMature = "+str(expT[k][0])
                _PlotC_(AC, id, "known", output)
        if k in N.keys():
            id = k+", novel, "
            AC = _Alignment_(canonical[k])
            if stars[k] != []:
                id += "ExpMature = "+str(expT[k][0])+", ExpStar = "+str(expT[k][1])
                AS = _Alignment_(stars[k])
                _PlotS_(AC, AS, id, "novel", output)
            if stars[k] == []:
                id += "ExpMature = "+str(expT[k][0])
                _PlotC_(AC, id, "novel", output)

def _main_():

    if len (sys.argv) != 5:
        print("Basic usage: plot_function.py known_mirnas.fa novel_mirnas.fa mirna_clusters.txt output_folder")
        print("\t> known_mirnas.fa: fasta with known miRNAs identified")
        print("\t> novel_mirnas.fa: fasta with novel miRNAs identified")
        print("\t> mirna_clusters.txt: clusters of miRNA generated in the pipeline")
        print("\t> output_folder: path/to/output/folder (use \".\" to set the current directory)")
        quit()

    known = sys.argv[1]
    novel = sys.argv[2]
    CD = sys.argv[3]
    output = sys.argv[4]
    if output == ".":
        output = str(os.getcwd())+"/"
    _GenOutput_(known, novel, CD, output)

_main_()

#END