#!/usr/bin/env python
#script designed to BLAST search the RFAM and filter out tRNAs and rRNAs from the small RNA-seq dataset in reads_collapsed.fa format
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys, os
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _GenOutput_(reads, rfam, output):
    #parse reads
    READS = _ParseFasta_(reads)
    #run blast search
    os.system("blastn -db "+rfam+" -outfmt 6 -query "+reads+" -out "+output+"rfam_blast.out -max_target_seqs 5 -num_threads 4")
    hit = set([])
    a = open(output+"rfam_blast.out","r")
    for line in a:
        line1 = line.strip().split("\t")
        hit.add(line1[0])
    a.close()
    #write output
    OUT = open(output+"reads_collapsed_filtered.fa","w")
    for k in READS.keys():
        if k not in hit:
            OUT.write(">"+k+"\n"+READS[k]+"\n")
    OUT.close()

def _main_():

    if len (sys.argv) != 4:
        print("Basic usage: filter_rfam.py reads_collapsed.fa RFAMdb output_folder")
        print("\t> reads_collapsed.fa: output from collapse_reads.py")
        print("\t> RFAMdb: path/to/RFAMdb")
        print("\t> output_folder: path/to/output/folder (use \".\" to set the current directory)")
        quit()

    reads = sys.argv[1]
    rfam = sys.argv[2]
    output = sys.argv[3]
    if output == ".":
        output = str(os.getcwd())+"/"
    _GenOutput_(reads, rfam, output)

_main_()

#END