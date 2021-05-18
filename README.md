![denomiR_logo](denomiR_logo.png)


# denomiR
**De no**vo identification of **miR**NAs. It is an experimental approach designed to identify miRNAs in small RNA-seq data from species with no genome available to be used as reference.

## Pipeline

- collapse 100% identical reads
- filter out tRNA and rRNA reads
- miRNA identification
	- cluster similar reads (isoforms with 98% identity)
	- detect the star reads (with 90% identity)
	- identify known and putative novel miRNAs

## Requirements

- [Python3](https://www.python.org/)
    - [Biopython](https://biopython.org/wiki/Download)
    - [NumPy](https://numpy.org/)
    - [Matplotlib](https://matplotlib.org/2.0.2/index.html)
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Installation

If the user has all requirements installed an working properly properly, you just need to do the following steps:
```
git clone https://github.com/pedronachtigall/denomiR.git
echo "export PATH=$PATH:$(pwd)/denomiR/bin/" >> ~/.bash_profile
source ~/.bash_profile
```

Alternatively, the user can install all requirements through conda manager as follow:
```
conda create --name denomiR_env -c bioconda python=3.7 biopython numpy matplotlib cd-hit blast
git clone https://github.com/pedronachtigall/denomiR.git
echo "export PATH=$PATH:$(pwd)/denomiR/bin/" >> ~/.bash_profile
source ~/.bash_profile
conda activate denomiR_env
```

It may be needed to apply "execution permission" to all bin executables:
  - ```chmod 777 path/to/denomiR/bin/*```

## Usage
```
Usage: denomiR.py [options]

Options:
  -h, --help            show this help message and exit
  -i string, --input=string
                        Mandatory - reads clipped in fq format (it can be
                        compressed in gz format also)
  -o path, --output=path
                        Optional - path to output folder [default =
                        "dnm_output"]
  -d path, --db=path    Mandatory - path/to/folder/with/DBs, folder containing
                        the RFAMdb and miRNAs from miRBase and MirGeneDB
  -c boolean, --CollapseReads=boolean
                        Optional - Turn On/Off the collpase reads step with
                        True/False. If your reads are in fasta format and were
                        collapsed elsewhere [default = True]
  -r boolean, --rfam=boolean
                        Optional - Turn On/Off the rfam filtering step with
                        True/False. [default = True]
  -l int, --MinLength=int
                        Optional - Minimum length of reads to used in the
                        analysis [default = 19]
  -L int, --MaxLength=int
                        Optional - Maximum length of reads to used in the
                        analysis [default = 25]
  -R int, --ReadsDepth=int
                        Optional - Minimum number of reads to be kept in the
                        analysis [default = 10]
  -s float, --similarity=float
                        Optional - Similarity presented by reads to be
                        clustered together. The number must be between 0.1 and
                        1.0 [default = 0.98]
```

Basic usage:
```
denomiR.py -i reads_clipped.fq(.gz) -o output_folder -db path/to/rfam&mirnaDB
```

To plot charts for each miRNA, please use the `plot_function.py` script as follow:
```
plot_function.py known_mirnas.fa novel_mirnas.fa mirna_clusters.txt output_folder
```

## Input

The denomiR was designed to work with the clipped reads in fastq format as input. The file can be compressed in gzip format.
Please notice, that the adapter must be trimmed. The user must use any tool available to perform this task.

The user must also indicate the PATH to the RFAM and miRNAdb `-db path/to/denomiR/mirDB/`.
Please, notice that the current DB is available in a zip file format. Decompress this file before running ```unzip mirDB.zip```
The current DB available is designed for Metazoa species, but the user can design a DB specific to their lineage and replace the files in the mirDB folder.

## Output
By default, all files generated during the analysis are kept. The `final_results.txt` is the final summary of results and containing read counts and the mature and star sequences identified for known and putative novel miRNAs. The sequences of miRNAs identified can be found in the fasta files `known_mirnas.fa` and `novel_mirnas.fa`. If the user also uses the `plot_function.py` a folder named "plots" will contain charts with the expression data in each position of the miRNA for each known and putative novel miRNA identified in the analysis.

```
denomiR_output/
├──  cdhit.clstr
├──  final_results.txt
├──  known_mirnas.fa
├──  mirna_clusters.txt
├──  mirna_seq.fa
├──  novel_mirnas.fa
├──  other_seqs.fa
├──  rfam_blast.out
└──  plots/
    ├──  known/
    |    ├── id_1.pdf
    |    ├── ...
    |    └── id_N.pdf
    └──  novel/
         ├── id_1.pdf
         ├── ...
         └── id_N.pdf
```

## Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

## Cite

