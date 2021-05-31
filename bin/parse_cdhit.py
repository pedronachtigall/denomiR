#!/usr/bin/env python
#Script designed to parse the .clstr output file from cd-hit-est-2d
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys, os
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id).split("_x")[0]
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseCDHIT_(cdhit):
    clusters = {}
    mir = ""
    hit = ""
    a = open(cdhit,"r")
    for line in a:
        if line.startswith("0\t"):
            line1 = line.strip().replace("... *","").split(", >")
            line2 = line1[1].split("_x")
            mir = line2[0]
            exp = line2[1]
            clusters.setdefault(mir, [exp])
        if line.startswith("1\t"):
            line1 = line.strip().split(", >")
            line2 = line1[1].split("...")
            hit = line2[0].split("_")[0]
        if line.startswith(">"):
            clusters.setdefault(mir, [])
            if hit != "":
                clusters[mir].append(hit)
            mir = ""
            hit = ""
    a.close()
    del clusters['']
    return clusters

def _GenOutput_(fasta, cdhit, output):
    C = _ParseCDHIT_(cdhit)
    F = _ParseFasta_(fasta)
    WithHit = {}
    NoHitWithStar = {}
    NoHitNoStar = {}
    for k in C.keys():
        if "canonical" in k:
            id = k.replace("_canonical", "")
            if len(C[k]) > 1:
                exp = C[k][0]
                WithHit[id] = [C[k][1], exp, F[k]]
                if k.replace("canonical", "star") in C.keys():
                    expS = C[k.replace("canonical", "star")][0]
                    WithHit[id].append(expS)
                    WithHit[id].append(F[k.replace("canonical", "star")])
                if k.replace("canonical", "star") not in C.keys():
                    WithHit[id].append("-")
                    WithHit[id].append("-")
            if len(C[k]) == 1:
                exp = C[k][0]
                if k.replace("canonical", "star") in C.keys():
                    expS = C[k.replace("canonical", "star")][0]
                    seqS = F[k.replace("canonical", "star")]
                    NoHitWithStar[id] = ["-",exp,F[k], expS, seqS]
                if k.replace("canonical", "star") not in C.keys():
                    NoHitNoStar[id] = ["-", exp, F[k],"-","-"]

    OUT = open(output+"final_results.txt","w")
    OUT.write("Cluster of reads with hits at miRBase and/or MirGeneDB (known miRNAs)\n")
    OUT.write("ID\tBestHit\tCanonicalRead\tCanonicalSeq\tStarRead\tStarSeq\n")
    OUTF = open(output+"known_mirnas.fa","w")
    for k, v in WithHit.items():
        OUT.write(k+"\t"+"\t".join([str(x) for x in v])+"\n")
        OUTF.write(">"+k+"_"+v[0]+"\n"+v[2]+"\n")
    OUTF.close()
    OUT.write("\nCluster of reads with no hits, but presents star reads also (probability of being novel miRNAs)\n")
    OUT.write("ID\tBestHit\tCanonicalRead\tCanonicalSeq\tStarRead\tStarSeq\n")
    OUTF = open(output+"novel_mirnas.fa","w")
    for k, v in NoHitWithStar.items():
        OUT.write(k+"\t"+"\t".join([str(x) for x in v])+"\n")
        OUTF.write(">"+k+"\n"+v[2]+"\n")
    OUTF.close()
    OUT.write("\nCluster of reads with no hits and no star reads (low-probability of being miRNAs) - Other seqs\n")
    OUT.write("ID\tBestHit\tCanonicalRead\tCanonicalSeq\tStarRead\tStarSeq\n")
    OUTF = open(output+"other_seqs.fa","w")
    for k, v in NoHitNoStar.items():
        OUT.write(k+"\t"+"\t".join([str(x) for x in v])+"\n")
        OUTF.write(">"+k+"\n"+v[2]+"\n")
    OUTF.close()
    OUT.close()

def _main_():

    if len (sys.argv) != 4:
        print("Basic usage: parse_cdhit.py mirna_seq.fa cdhit.clstr output_folder")
        print("\t> mirna_seq.fa: mirna_seq.fa")
        print("\t> cdhit.clstr: output from cd-hit-est-2d")
        print("\t> output_folder: path/to/output/folder (use \".\" to set the current directory)")
        quit()

    fasta = sys.argv[1]
    cdhit = sys.argv[2]
    output = sys.argv[3]
    if output == ".":
        output = str(os.getcwd())+"/"
    _GenOutput_(fasta, cdhit, output)

_main_()

#END
