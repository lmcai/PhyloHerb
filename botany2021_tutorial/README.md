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
export PH=../PhyloHerb
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
sp3     666666  35628.0 43.1    133223  0.35511135464596955     No
sp4     666666  24956.0 25.4    159410  0.3685026033498526      No
sp5     666666  20474.0 35.3    131162  0.357679815800308       No
```

For nuclear ribosomal regions and mitochondrial assemblies:
```
python $PH/phyloherb.py -m qc -s sample_sheet.tsv -i ./ITS -o ../2_assemblies/ITS
python $PH/phyloherb.py -m qc -s sample_sheet.tsv -i ./mito -o ../2_assemblies/mito
```

## VI. Annotation and organelle structure variations

We will upload the circularized assembly to the web-based annnotator [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html). The annotations will be returned as genbank files, graphics, and alignments.

The genbank files can be imported to Geneious for visualization and manual curation.

## V. Ortholog identification and alignment
1. Ortholog identification

We will use the build-in references to extract plastid genes. The list of our reference species is [here](/database/plastid_reference_sp.csv). The list of the genes in the database is [here](/database/plastid_gene.list). We will only use blast hits longer than 120 bp.

Makesure executable `blast` commands are loaded in your environment. In the `herbariomics_workshop` directory, type the following commands
```
python $PH/phyloherb.py -m ortho -i 2_assemblies/chl -o 3_alignments/chl/ -l 120
```

In the output directory `3_alignments/chl`, orthologous genes will be written to separate fasta files and the headers will be the species prefixes.

3. Alignment

I like to use the `--adjustdirection` function from `MAFFT` to correct reverse complimentary sequences. Then I will use `pasta` to more accurately align high variable sequences such as the intergenic regions and the ITS regions. `pasta` first generates a guidance tree, then align among closely-related species, finally merge the alignments to produce the output.

This is a potentially time consuming step so I recommend running it on the cluster using the example batch file [mafft_pasta.sh](phyloherbLib/mafft_pasta.sh).

Copy `mafft_pasta.sh` to the same directory where the gene sequences are located. Modify the file to include appropriate environmental parameters. Then the batch job can be submitted to the cluster by typing
```
sbatch mafft_pasta.sh <gene_1>
sbatch mafft_pasta.sh <gene_2>
```

4. Intergenic regions


5. Nuclear ribosomal and mitochondrial regions

The nuclear ribosomal data requires a slightly different curation strategy. The highly variable sequence requires more manual curation than the plastome. The nuclear ribosomal region exists as tandem repeats on multiple chromosomes.

<img src="/images/ITS.png" width="400" height="400">

Based on our experiences, NTS is not alignable even between closely related taxa. The entire rDNA region (18S+ITS1+5.8S+ITS2+25S) and some portion of ETS can be aligned at family level. 

We will try to align the entire region using MAFFT first.   

6. mitochondrial regions

For most plant groups, mitochondria are not phylogenetically informative because the genes evolve too slowly, but the intergenic regions are highly variable. Moreover, the qualities of mitochondrial genomes are usually not as good as plastomes. So we will only extract mitochondrial genes for comparative purposes. The methods is similar to plastid genes.


7. Manual curation in Geneious

At this point, it is recommended that you take a initial look at your alignments. **Initial** means be prepared to complete the "alignment-manual check-phylogeny" cycle for at least two rounds to get publication quality data.

The purpose of the initial check is to remove obvious low-quality sequences. Do not conduct any site-based filtering yet! For example, the two sequences highlighted in red below contain too many SNPs (marked in black). They should be removed.

Geneious provides a nice interface to work with alignments. You can view alignment statistics, edit sequences, and concatenate alignments. Alternative automatic tools include [trimAL](http://trimal.cgenomics.org/getting_started_with_trimal_v1.2) and [Seaview](http://doua.prabi.fr/software/seaview).


## VI. Phylogeny reconstruction

1. Concatenation

Many tools are available for concatenating alignments. I recommend the `conc` function of phyloherb or Geneious. I have applied both tools to dataset with 1000 sp x 100 genes. The `conc` function of phyloherb will also output a gene delineation file required by `PartitionFinder`.

To use the `conc` function of phyloherb, use the following commands
```
python phyloherb.py -m conc -i <directory containing alignments> -o <output directory> -suffix <suffix>
```
This command will concatenate all of the fasta sequences in the input directory with the specified suffix. Again, if you only want to use a subset of the genes or want the genes to appear in a specific order, you can supply a gene order file by adding `-g gene_subset.txt`.

2. Maximum likehood phylogeny

For an initial quick and dirty analysis, I recommend ExaML with unpartitioned alignment. IQTREE or RAxML generates more accurate estimations of the phylogeny and substitution paramters, but may not accomodate thousands of species with millions of sites. The commands I use for ExaML analysis is as follows.
```
#generate a parsimony tree with RAxML
raxmlHPC-SSE3 -y -m GTRGAMMA -s DNA_algnment.fas -n test -p 3256179
parse-examl -s DNA_aln.phy -n test -m GTRGAMMA
examl-AVX -S -s test.binary -m GAMMA -n testML -t RAxML_parsimonyTree.test
```

3. Second round of manual alignment curation

It can be a quite satisfying experience when you browse through a well-curated alignment. To get us there, we need to conduct a second round of alignment curation and remove spurious regions arising from assembly errors or false positive BLAST hits. 

First, using a reference ExaML species tree (newick format), we will order the sequences based on their phylogenetic affinity. This will facilitate manual curation of the alignments in Geneious because you can see shared mutations in closely related species.

The `order` function of phyloherb takes a reference tree and reorders all alignments in the input directory based on the phylogeny. If you want to additionally filter sequences based on missing data, using the optional `-missing` flag. A float number from 0 to 1 is required for `-missing` to indicate the maximum proportion of ambiguous sites allowed for each sequence.
```
python phyloherb.py -m order -t <reference.tre> -i <directory containing alignments> -o <output directory> -suffix <suffix> -missing <float number 0 to 1>
```

This will generate an ordered alignment `*.ordered.fas` and a companion tree file `*.pasta_ref.tre` for each gene. You will need this tree for the second round of pasta alignment after manual curation.

Now let's load the ordered alignments to Geneious for some fine tuning. This time we will delete blocks of problematic sequences. They usually appears as a cluster of SNPs highlighted in red below. These SNPs are not conserved in their close relatives, so they are phylogenetically uninformative autapomorphies. Regardless of the causes, we can safely delete them. 

<img src="/images/Geneious2.png" width="400" height="250">

After a second manual check, your alignments is ready for re-alignment in `pasta`. This time we will use a reference tree `*.pasta_ref.tre` to guide the alignment for each gene.
```
run_pasta.py -i <input sequence> -a -t <reference tree> -o <output directory>
```

To submit batch job to the cluster, you can modify [this](phyloherbLib/pasta_rnd2.sh) batch file. Make sure you have loaded the correct environment.

When the alignment is done, you can use them to produce your final species tree.
