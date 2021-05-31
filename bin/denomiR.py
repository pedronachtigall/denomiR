#!/usr/bin/env python
#script designed to perform de novo miRNA identification
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import os
import datetime as dt
from optparse import OptionParser
from difflib import SequenceMatcher
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
except:
    print("Biopython module is not installed.\nInstall it using \"pip install biopython\" or check https://biopython.org for more details.")

def _ReverseSeq_(string):
    my_seq = Seq(string)
    return str(my_seq.reverse_complement())

def _ParseFasta_(fasta, minL, maxL, reads):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        if int(ID.split("_x")[1]) > int(reads) and len(SEQ) >= int(minL) and len(SEQ) <= int(maxL) and "N" not in SEQ:
            final[ID] = SEQ
    return final

def _CollapseReads_(fasta, output):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Collapsing reads...")
    os.system("collapse_reads.py "+fasta+" "+output)

def _FilterRFAM_(fasta, output, db):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Filtering tRNAs and rRNAs reads...")
    os.system("filter_rfam.py "+fasta+" "+db+"RFAM "+output)

def _RunCDHIT_(output, db, similarity):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Refinement of clustering miRNAs and identification...")
    os.system("cd-hit-est-2d -d 120 -i "+output+"mirna_seq.fa -i2 "+db+"metazoa_mature_MBMG.fa -c 0.98 -o "+output+"cdhit")
    os.system("parse_cdhit.py "+output+"mirna_seq.fa "+output+"cdhit.clstr "+output)

def _ClusterMIR_(fastaIN, output, db, collapse, rfam, similarity, minL, maxL, reads):
    FastaFiltered = fastaIN
    if collapse == "True" or collapse == "TRUE" or collapse == "true":
        _CollapseReads_(FastaFiltered, output)
        FastaFiltered = output+"reads_collapsed.fa"
    if rfam == "True" or rfam == "TRUE" or rfam == "true":
        _FilterRFAM_(FastaFiltered, output, db)
        FastaFiltered = output+"reads_collapsed_filtered.fa"
    F = _ParseFasta_(FastaFiltered, minL, maxL, reads)
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Clustering step...")
    reads = list(F.keys())
    computed = []
    clusters = {}
    for i in reads:
        expi = int(i.split("_x")[1])
        #use only SEQs with more than 100 reads
        if i not in computed and expi > 100:
            for j in reads:
                expj = int(j.split("_x")[1])
                if j != i and j not in computed and expj > 10:
                    seq1 = F[i]
                    seq2 = F[j]
                    seq2rev = _ReverseSeq_(seq2)
                    #same direction
                    if SequenceMatcher(a=seq1,b=seq2).ratio() >= 0.90:
                        alignments = pairwise2.align.globalms(seq1, seq2,  5, -4, -4, -1)
                        gaps = alignments[0][0].count("-")+alignments[0][1].count("-")
                        matches = len(alignments[0][0])-gaps
                        percid = (matches*100)/len(alignments[0][0])
                        percidS = (matches*100)/min([len(seq1), len(seq2)])
                        #if percidS >= 98.0:
                        if percidS >= float(similarity)*100:
                            #print(i,j,"same strand")
                            #print(percid, percidS, matches, gaps, len(alignments[0][0]), min([len(seq1), len(seq2)]), max([len(seq1), len(seq2)]))
                            #print(alignments[0])
                            #print(SequenceMatcher(a=seq1,b=seq2).ratio())
                            computed.append(i)
                            computed.append(j)
                            exp1 = i.split("_x")[1]
                            exp2 = j.split("_x")[1]
                            clusters.setdefault(i, [[i,seq1,exp1,"canonical"]])
                            clusters[i].append([j,seq2,exp2,"canonical"])
                            continue
                    #reverse complement
                    if SequenceMatcher(a=seq1,b=seq2rev).ratio() >= 0.80:
                        alignments = pairwise2.align.globalms(seq1, seq2rev,  5, -4, -4, -1)
                        #alignments = pairwise2.align.localxx(seq1, seq2rev)
                        gaps = alignments[0][0].count("-")+alignments[0][1].count("-")
                        matches = len(alignments[0][0])-gaps
                        percid = (matches*100)/len(alignments[0][0])
                        percidS = (matches*100)/min([len(seq1), len(seq2)])
                        #identity = matches
                        #print(percid, percidS, matches, gaps, len(alignments[0][0]), min([len(seq1), len(seq2)]), max([len(seq1), len(seq2)]))
                        #print(alignments[0])
                        #if percidS >= 98.0:
                        if percidS >= float(similarity)*100:
                            #print(i,j,"reversecomplement <<<<<<<<<<<<<")
                            #print(percid, percidS, matches, gaps, len(alignments[0][0]), min([len(seq1), len(seq2)]), max([len(seq1), len(seq2)]))
                            #print(alignments[0])
                            #print(SequenceMatcher(a=seq1,b=seq2rev).ratio())
                            #print(percid, matches, gaps, len(alignments[0][0]))
                            #print(alignments[0])
                            computed.append(i)
                            computed.append(j)
                            exp1 = i.split("_x")[1]
                            exp2 = j.split("_x")[1]
                            clusters.setdefault(i, [[i,seq1,exp1,"canonical"]])
                            clusters[i].append([j,seq2,exp2,"star"])
                            continue
        computed.append(i)

    OUT1 = open(output+"mirna_seq.fa","w")
    OUT2 = open(output+"mirna_clusters.txt", "w")
    count = 1
    for k in clusters.keys():
        stars = []
        expT = 0
        expTS = 0
        for i in clusters[k]:
            if i[-1] == "canonical":
                expT += int(i[2])
            if i[-1] == "star":
                stars.append(i)
                expTS += int(i[2])
        OUT1.write(">id_"+str(count)+"_canonical_x"+str(expT)+"\n"+clusters[k][0][1]+"\n")
        if stars != []:
            OUT1.write(">id_"+str(count)+"_star_x"+str(expTS)+"\n"+stars[0][1]+"\n")
        OUT2.write(">id_"+str(count)+"\n")
        for i in clusters[k]:
            OUT2.write("\t"+" ".join(i)+"\n")
        count += 1
    OUT1.close()
    OUT2.close()

    _RunCDHIT_(output, db, similarity)

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Mandatory - reads clipped in fq format (it can be compressed in gz format also)", metavar="string", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - path to output folder [default = \"dnm_output\"]", metavar="path", default=None)
    parser.add_option("-d", "--db", dest="db", help="Mandatory - path/to/folder/with/DBs, folder containing the RFAMdb and miRNAs from miRBase and MirGeneDB", metavar="path", default=None)
    parser.add_option("-c", "--CollapseReads", dest="collapse", help="Optional - Turn On/Off the collpase reads step with True/False. If your reads are in fasta format and were collapsed elsewhere [default = True]", metavar="boolean", default="True")
    parser.add_option("-r", "--rfam", dest="rfam", help="Optional - Turn On/Off the rfam filtering step with True/False. [default = True]", metavar="boolean", default="True")
    parser.add_option("-l", "--MinLength", dest="minL", help="Optional - Minimum length of reads to used in the analysis [default = 19]", metavar="int", default="19")
    parser.add_option("-L", "--MaxLength", dest="maxL", help="Optional - Maximum length of reads to used in the analysis [default = 25]", metavar="int", default="25")
    parser.add_option("-R", "--ReadsDepth", dest="reads", help="Optional - Minimum number of reads to be kept in the analysis [default = 10]", metavar="int", default="10")
    parser.add_option("-s", "--similarity", dest="similarity", help="Optional - Similarity presented by reads to be clustered together. The number must be between 0.1 and 1.0 [default = 0.98]", metavar="float", default="0.98")

    (options, args) = parser.parse_args()

    if float(options.similarity) > 1.0:
        print("Error: You set the similarity threshold in a bad way, it must be between 0.1 and 1.0 (e.g., 0.98 or 0.97 or 0.90 ...).")

    if options.input == None or options.db == None:
        print(
        """
>>>> De novo analysis of miRNA seq <<<<
      ****Use -h for help!****

USAGE:
denomiR.py -i reads_clipped.fq(.gz) -o output_folder -db path/to/rfam&mirnaDB
""")
        quit()

    if options.input != None and options.db != None:
        print("""

>>>> denomiR - De novo analysis of miRNA-seq <<<<

        """)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Starting denomiR...")
        if options.output != None:
            if not options.output.endswith("/"):
                options.output += "/"
        if options.output == None:
            options.output = os.getcwd()+"/dnm_output/"
        if os.path.isdir(options.output) == False:
            os.mkdir(options.output)
        print("\tInput file -> "+options.input)
        print("\tOutput folder -> "+options.output)
        print("\tRFAM/miRNA DBs -> "+options.db)
        print("\tMinimum length -> "+options.minL)
        print("\tMaximum length -> "+options.maxL)
        print("\tMinimum reads number -> "+options.reads)
        print("\tSimilarity threshold -> "+options.similarity)

        _ClusterMIR_(options.input,
                    options.output,
                    options.db,
                    options.collapse,
                    options.rfam,
                    options.similarity,
                    options.minL,
                    options.maxL,
                    options.reads)

        print("\n\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Finished!!!")

__main__()



#END
