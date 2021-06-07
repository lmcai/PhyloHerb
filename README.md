# PhyloHerb
**Phylo**genomic Analysis Pipeline for **Herb**arium Specimens

This bioinformatic pipeline provides detailed guidance to process **genome skimming** data collected from herbarium specimens. The outcomes include plastid genome assemblies, mitochondrial genome assemblies, nuclear 35S ribosomal DNAs (NTS+ETS+18S+ITS1+5.8S+ITS2+25S), alignments of gene and intergenic regions, and a species tree. Combined with the morphological and distribution data from herbarium specimens, this approach provides an unparalleled opportunity to study **taxonomy, biogeography, and macroevolution with nearly complete taxon sampling**.

We have tested this pipeline in the Barbados Cherry family Malpighiaceae, Clusiaceae, and several groups of algae. Each of these datasets contains hundreds to thousands of species and our pipeline extracts ample data to resolve both recent radiations (e.g., *Bunchosia*, Malpighiaceae >135 sp within 10 Myr) and ancient divergences (e.g., the divergence of red algea at hundreds of millions of years ago). 

## Prerequisites
To process large datasets (>20 sp), high performance cluster is recommended. Mac and PC can suffer from insufficient memory during the assembly, alignment, and phylogenetic reconstruction.

### Assembly
1. [GetOrganelle v1.7.0+](https://github.com/Kinggerm/GetOrganelle)
2. [Bowtie2 v2.2.2+](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
3. [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
4. Assembly viewer: [Bandage](https://rrwick.github.io/Bandage/)
5. Optional assembly viewer: [Geneious](https://www.geneious.com/) (the complementary version is sufficient)

### Alignment
6. [Biopython](https://biopython.org/)
7. Aligner: 
	[Pasta](https://github.com/smirarab/pasta) for highly variable regions such as the ITS sequences
	
	[MAFFT](https://mafft.cbrc.jp/alignment/software/) for less variable regions or long alignments (>5 kb) that pasta may not be able to handle when the number of species is high (>500 sp)

### Phylogeny
8. [IQ-TREE](http://www.iqtree.org/)

## General guidelines for genome skimming data collection
About 1-3% of the reads from genome skimming (low-coverage genome sequencing) are from plastids. Theoretically this value should vary with the size of the nuclear genome and the abundance of plastids within a cell. 

![plastic_perc_equation](/images/plastid_perc.png){:height="50%" width="50%"}

## Step 1: Filter adapters with blastn
Choose a taxonID for each data set. This taxonID will be used throughout the analysis. Use short taxonIDs with no special characters.

Remove adapters for paired end reads:

	python filter_fastq.py taxonID_1.fq taxonID_2.fq adapter_file num_cores

Alternatively, if sequences are single end reads:
	
	python filter_fastq.py taxonID.fq adapter_file num_cores

The output files are taxonID_1.fq.filtered and taxonID_2.fq.filtered for paired end reads, and taxonID.fq.filtered for single end reads. I use filter_fastq.py only for sequences downloaded from SRA with adapter sequences and phred score offset unknown. For data sets with known adapter and phred scores I use Trimmomatic instead (see below). Trimmomatic is much faster but would not give an error message if you the adapter or phred score offset were misspedified.

##Step 2: 

