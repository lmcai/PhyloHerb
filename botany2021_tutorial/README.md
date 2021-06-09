# Botany 2021 workshop tutorial

This is a step-by-step tutorial to prossess genome skimming data for phylogeny reconstruction. This tutorial will use the example dataset hosted on the PhyloHerb repository to reconstruct a plastome-based phylogeny for five species.

## I. Prerequisites

Install the following software

1. Install [GetOrganelle v1.7.0+](https://github.com/Kinggerm/GetOrganelle), [Bowtie2 v2.2.2+](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) all at once with conda

Install Anaconda or miniconda for Linux, Mac, or PC (see https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

Create a conda environment and install GetOrganelle and its dependencies (may take a while)
```
conda create --name getorg
#answer promt questions 
#proceed ([y]/n)?
#y
conda install -c bioconda getorganelle
conda install -c conda-forge biopython
```
2. [Bandage](https://rrwick.github.io/Bandage/): Assembly graph viewer with a graphical user interface. Can be downloaded and installed based on the instruction online.

3. Optional [Geneious](https://www.geneious.com/): the free trial version lasts for 14 days. After that some functions are restricted. Your institute may have purchased a license for you.

4. [Pasta](https://github.com/smirarab/pasta) or [MAFFT](https://mafft.cbrc.jp/alignment/software/). Pasta is highly recommended for distantly related taxa and large datasets. I recommend installing it with Docker image. Installation from source can be somewhat difficult and may require root permissions. If you plan to use it during the workshop, please make sure you have successfully installed it prior to the workshop.

5. [IQ-TREE](http://www.iqtree.org/) or [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) or [ExaML](https://cme.h-its.org/exelixis/web/software/examl/index.html). Choose your favorite one!

## II. Assembly

We will use [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) to assemble the plastome, mitochondrial genome, and ribosomal regions. 

### 1. Input:

Illumina FASTQ reads for each species, single-ended or pair-ended, zipped or unzipped. Do not filter the reads or trim adapters, GetOrganelle will take care of it.

### 2. How to:

Load dependencies (assuming installing GetOrganelle under the conda environment named 'getorg')
```
#load Anaconda
module load Anaconda
#load GetOrganelle
source activate getorg
```
Create a working directory for assembly and enter it

```
mkidr 1_getorg
cd 1_getorg
mkdir chl
mkdir ITS
mkdir mito
```
Download example datasets
```
wget 
```

The basic commands for running assembly with pair end data is as follows:

```
#To assemble plant plastome
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o chl/<output prefix> -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o ITS/<output prefix> -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o mito/<output prefix> -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
```

For this dataset, each plastome assembly takes ~600MB memory, and ~60s CPU time.


### 3. Large dataset and batch submission to cluster

If you are working with a high performance computing cluster with slurm workload manager. You can modify the [bash file](/utilities/getorg.sh) and submit jobs to your cluster simultaneously.

The bash job looks like this:
```
#!/bin/bash
#
#SBATCH -n 8                 # Number of cores
#SBATCH -N 1                 # Number of nodes for the cores
#SBATCH -t 0-10:00           # Runtime in D-HH:MM format
#SBATCH -p <name of the partition>    # Partition to submit to
#SBATCH --mem=8000            # Memory pool for all CPUs
#SBATCH -o getorg.out      # File to which standard out will be written
#SBATCH -e getorg.err      # File to which standard err will be written


#Activate conda environment, assuming the name of the conda environment is 'getorg'
#module load Anaconda
#source activate getorg

#export DATA_DIR=<absolute path to input data>

#To assemble plant plastome
get_organelle_from_reads.py -1 $1 -2 $2 -o chl/$3 -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 $1 -2 $2 -o ITS/$3 -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 $1 -2 $2 -o mito/$3 -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt

```
*IMPORTANT*: Make sure you modify the above bash file to provide correct path to the data, request correct resources, etc (marked by "<>").

To submit the job, you can provide a tab-delimited text file [sample_sheet.tsv](/example/sample_sheet.tsv). The `sp_prefix` will be used throughout and  in the final phylogeny.
```
sp_prefix	Forward_reads	Reverse_reads
sp1	sp1.100m.R1.fq.gz	sp1.100m.R2.fq.gz
sp2	sp2.100m.R1.fq.gz	sp2.100m.R2.fq.gz
sp3	sp3.100m.R1.fq.gz	sp3.100m.R2.fq.gz
sp4	sp4.100m.R1.fq.gz	sp4.100m.R2.fq.gz
sp5	sp5.100m.R1.fq.gz	sp5.100m.R2.fq.gz
```

Then generate a batch file using the submission function of phyloherb
```
#python phyloherb.py -a submision -b <batch file name> -s <sample sheet> -o <output>
python phyloherb.py -a submision -b getorg.sh -s sample_sheet.tsv -o submitter.sh

#submit jobs
sh submitter.sh
```

Or you can submit your job one by one:

```
sbatch getorg.sh <forward.fq> <backward.fq> <output prefix>

```
### 4. Output

The batch submission will generate three subdirectories `chl/`, `ITS/`, and `mito/`, each containing Getorganelle output directories named after sample-specific prefixes.

### 5. Assembly QC and visualization with Bandage

After the assemblies are completed, you can summarize the results using the QC function of phyloherb. For each species, it will extract the following information: the number of total input reads, the number of reads used for assembly, the total length of the assembly, GC%, average base coverage, and whether the assembly is circularized. 

```
python phyloherb.py -a qc -s sample_sheet.tsv -d ./ -o assembly_sum.tsv
```


## VI. Annotation and organelle structure variarion

## V. Alignment generation
