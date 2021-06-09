# Botany 2021 workshop tutorial

This is a step-by-step tutorial to prossess genome skimming data for phylogeny reconstruction. This tutorial will use the example dataset hosted on the PhyloHerb repository to reconstruct a plastome-based phylogeny for five species.

## I. Prerequisites

Install the following software

1. [GetOrganelle v1.7.0+](https://github.com/Kinggerm/GetOrganelle), [Bowtie2 v2.2.2+](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Install Anaconda or miniconda (see https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

Install GetOrganelle and its dependencies under a conda environment
```
conda create --name getorg
#answer promt questions 
#proceed ([y]/n)?
#y
conda install -c bioconda getorganelle
```
2. [Bandage](https://rrwick.github.io/Bandage/) has a graphical user interface. Can be downloaded and installed based on the instruction .

3. Optional [Geneious](https://www.geneious.com/): the free trial lasts for 14 days, after that some functions are restricted. Some institutes purchased licenses.

### Alignment
6. [Biopython](https://biopython.org/)
7. Aligner: 

	[Pasta](https://github.com/smirarab/pasta) for highly variable regions such as the ITS sequences
	
	[MAFFT](https://mafft.cbrc.jp/alignment/software/) for less variable regions or long alignments (>5 kb) that pasta may not be able to handle when the number of species is high (>500 sp)
8. Manual assembly examination: [Geneious](https://www.geneious.com/) (the licensed version are required)

### Phylogeny
9. [IQ-TREE](http://www.iqtree.org/) or [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) or [ExaML](https://cme.h-its.org/exelixis/web/software/examl/index.html)

## II. Assembly

We will use [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) to assemble the plastome, mitochondrial genome, and ribosomal regions. 

### 1. Input:

Illumina FASTQ reads for each species, single-ended or pair-ended, zipped or unzipped. Do not filter the reads or trim adapters, GetOrganelle will take care of it.

### 2. How to:

Load dependencies (assuming installing GetOrganelle under the conda environment)
```
module
```
Create a working directory for assembly and enter it

```
mkidr 1_getorg
cd 1_getorg
```
After loading GetOrganelle to your environment, the basic commands for running assembly with pair end data is as follows:

```
#To assemble plant plastome
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <plastome_output> -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <nr_output> -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <mito_output> -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
```

### 3. Large dataset and batch submission to cluster

If you are dealing with large number of species, running them one by one is too tedious. Here, we will submit individual assembly task to the cluster and run them simultaneously. An example bash file is provided in `/utilities/getorg.sh`. We will also use short and informative output prefix for each species. You can submit your job by typing

```
sbatch getorg.sh <forward.fq> <backward.fq> <output prefix>
```
*IMPORTANT*: Make sure you load the correct environment and provide absolute path to the input data if they are not in the current directory by modifying relavant variables in `getorg.sh`. Instructions for single-end data can also be found in `getorg.sh`.

### 4. Output

The batch submission will generate three subdirectories `chl/`, `ITS/`, and `mito/`, each containing Getorganelle output directories named after sample-specific prefixes.

## VI. Annotation and organelle structure variarion

## V. Alignment generation

