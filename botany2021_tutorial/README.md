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
conda install -c etetoolkit ete3
```
2. [Bandage](https://rrwick.github.io/Bandage/): Assembly graph viewer with a graphical user interface. Can be downloaded and installed based on the instruction online.

3. Optional [Geneious](https://www.geneious.com/): the free trial version lasts for 14 days. After that some functions are restricted. Your institute may have purchased a license for you.

4. [Pasta](https://github.com/smirarab/pasta) or [MAFFT](https://mafft.cbrc.jp/alignment/software/). Pasta is highly recommended for distantly related taxa and large datasets. I recommend installing it with Docker image. Installation from source can be somewhat difficult and may require root permissions. If you plan to use it during the workshop, please make sure you have successfully installed it prior to the workshop.

5. [IQ-TREE](http://www.iqtree.org/) or [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) or [ExaML](https://cme.h-its.org/exelixis/web/software/examl/index.html). Choose your favorite one!

## II. Assembly

We will use [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) to assemble the plastome, mitochondrial genome, and ribosomal regions. 

### 1. Input preparation:

Download PhyloHerb and example datasets from Github
```
git clone https://github.com/lmcai/PhyloHerb.git
mkdir herbariomics_workshop
cd herbariomics_workshop
export $PH=../PhyloHerb
```

Move and rename the example dataset from PhyloHerb to the current working directory
```
mkdir 0_fastq
mv ../PhyloHerb/example/*.gz 0_fastq
```

### 2. How to:

Load dependencies (assuming installing GetOrganelle under the conda environment named 'getorg')
```
#load Anaconda and GetOrganelle
module load Anaconda
source activate getorg
```
Create a working directory for assembly

```
mkdir 1_getorg
cd 1_getorg

mkdir chl
mkdir ITS
mkdir mito
```

The basic commands for running assembly with pair end data is as follows:

```
export DATA_DIR=../0_fastq
#To assemble plant plastome
get_organelle_from_reads.py -1 $DATA_DIR/SP1_R1.100m.1.fastq.gz -2 $DATA_DIR/SP1_R2.100m.1.fastq.gz -o chl/sp1 -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 $DATA_DIR/SP1_R1.100m.1.fastq.gz -2 $DATA_DIR/SP1_R2.100m.1.fastq.gz  -o ITS/sp1 -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 $DATA_DIR/SP1_R1.100m.1.fastq.gz -2 $DATA_DIR/SP1_R2.100m.1.fastq.gz -o mito/sp1 -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
```

For this dataset, each plastome assembly takes ~5GB memory, and ~120s CPU time.


### 3. Large dataset and batch submission to cluster

If you are working with a high performance computing cluster with slurm workload manager. You can modify the [bash file getorg.sh](/phyloherbLib/getorg.sh) and submit jobs to your cluster simultaneously.

Copy the example `getorg.sh` to current directory
```
cp $PH/phyloherbLib/getorg.sh .
```

The bash job looks like this:
```
#!/bin/bash
#
#SBATCH -n 8                 # Number of cores
#SBATCH -N 1                 # Number of nodes for the cores
#SBATCH -t 0-10:00           # Runtime in D-HH:MM format
#SBATCH -p <partition>    # Partition to submit to
#SBATCH --mem=10000            # Memory pool for all CPUs
#SBATCH -o getorg.out      # File to which standard out will be written
#SBATCH -e getorg.err      # File to which standard err will be written


#Activate conda environment, assuming the name of the conda environment is 'getorg'
#module load Anaconda
#source activate getorg

export DATA_DIR=<absolute path to input data>

#To assemble plant plastome
get_organelle_from_reads.py -1 $DATA_DIR/$1 -2 $DATA_DIR/$2 -o chl/$3 -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 $DATA_DIR/$1 -2 $DATA_DIR/$2 -o ITS/$3 -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 $DATA_DIR/$1 -2 $DATA_DIR/$2 -o mito/$3 -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt

```
*IMPORTANT*: Make sure you modify the above bash file to load correct environment, provide full path to the data, request correct resources, etc (marked by "<>").

To submit the job, you can provide a tab-delimited text file [sample_sheet.tsv](/example/sample_sheet.tsv). The `sp_prefix` will be used throughout and  in the final phylogeny.
```
sp_prefix	Forward_reads	Reverse_reads
sp1	SP1_R1.100m.1.fastq.gz	SP1_R2.100m.1.fastq.gz
sp2	SP2_R1.100m.1.fastq.gz	SP2_R2.100m.1.fastq.gz
sp3	SP3_R1.100m.1.fastq.gz	SP3_R2.100m.1.fastq.gz
sp4	SP4_R1.100m.1.fastq.gz	SP4_R2.100m.1.fastq.gz
sp5	SP5_R1.100m.1.fastq.gz	SP5_R2.100m.1.fastq.gz
```

Then generate a batch file using the submission function of phyloherb
```
cp $DATA_DIR/sample_sheet.tsv .
python $PH/phyloherb.py -m submision -b getorg.sh -s sample_sheet.tsv -o submitter.sh

#submit jobs
sh submitter.sh
```

### 4. Output

The batch submission will generate three subdirectories `chl/`, `ITS/`, and `mito/`, each containing Getorganelle output directories named after sample-specific prefixes. For detailed descriptions of the output, see [Getorganellel instructions](https://github.com/Kinggerm/GetOrganelle#Instruction)

### 5. Assembly visualization with Bandage


### 6. Assembly QC 

After the assemblies are completed, you can summarize the results using the QC function of phyloherb. For each species, it will extract the following information: the number of total input reads, the number of reads used for assembly, average base coverage, the total length of the assembly, GC%, and whether the assembly is circularized. 

```
python $PH/phyloherb.py -m qc -s sample_sheet.tsv -i ./chl -o ../2_assemblies/chl
```
This command will copy all of the assemblies under the input directory `./chl` to a new directory `../2_assemblies/chl` and rename the files based on their species prefixes. In the output directory, you will also find a summary spreadsheet `assembly_sum.tsv` that looks like this
```
sp_prefix       Total_reads     Reads_in_target_region  Average_base_coverage   Length  GC%     Circularized
sp1     666666  30400.0 31.4    159658  0.36851895927545125     Yes
sp2     666666  48728.0 -246.9  133603  0.3567809106082947      No
```

For nuclear ribosomal regions and mitochondrial assemblies:
```
python phyloherb.py -m qc -s sample_sheet.tsv -i ./ITS -o ../2_assemblies/ITS
python phyloherb.py -m qc -s sample_sheet.tsv -i ./mito -o ../2_assemblies/mito
```

## VI. Annotation and organelle structure variations

We will upload the circularized assembly to the web-based annnotator [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html). The annotations will be returned as genbank files, graphics, and alignments.

The genbank files can be imported to Geneious for visualization and manual curation.

## V. Ortholog identification
Phyloherb will identify the best-matching region of each gene/intergenic region in the assemblly using BLAST. We provide a build-in database of plastid genes from 100 angiosperm species. This database is sufficient for getting genes for most angiosperms 

