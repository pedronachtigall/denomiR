#!/usr/bin/env python
#script designed to collapse reads and convert from fq to fasta
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

import sys, os
import gzip
from Bio import SeqIO
from operator import itemgetter

def _ParseFastq_(fastq, output):
    final = {}
    for record in SeqIO.parse(fastq, "fastq"):
        SEQ = str(record.seq)
        final.setdefault(SEQ, 0)
        final[SEQ] += 1

    OUT = open(output,"w")
    count = 1
    finalS = sorted(final.items(), key=itemgetter(1), reverse=True)
    for k in finalS:
        OUT.write(">Seq_"+str(count)+"_x"+str(final[k[0]])+"\n"+k[0]+"\n")
        count += 1
    OUT.close()

def _ParseFastqGZ_(fastq, output):
    final = {}
    with gzip.open(fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            SEQ = str(record.seq)
            final.setdefault(SEQ, 0)
            final[SEQ] += 1

    OUT = open(output,"w")
    count = 1
    finalS = sorted(final.items(), key=itemgetter(1), reverse=True)
    for k in finalS:
        OUT.write(">Seq_"+str(count)+"_x"+str(final[k[0]])+"\n"+k[0]+"\n")
        count += 1
    OUT.close()

def _GenOutput_(fastq, output):
    out = output+"reads_collapsed.fa"
    if fastq.endswith(".fastq.gz"):
        _ParseFastqGZ_(fastq, out)
    if fastq.endswith(".fastq"):
        _ParseFastq_(fastq, out)

def _main_():

    if len (sys.argv) != 3:
        print("Basic usage: collapse_reads.py reads_clipped.fastq output_folder")
        print("\t> reads_clipped.fastq: reads with adapters clipped and low-complexity reads filtered (it can be in fastq or fastq.gz format)")
        print("\t> output_folder: path/to/output/folder (use \".\" to set the current directory)")
        quit()

    reads = sys.argv[1]
    output = sys.argv[2]
    if output == ".":
        output = str(os.getcwd())+"/"
    _GenOutput_(reads, output)

_main_()

#END
