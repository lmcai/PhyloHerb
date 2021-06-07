# PhyloHerb
**Phylo**genomic Analysis Pipeline for **Herb**arium Specimens

This bioinformatic pipeline provides detailed guidance to process genome skimming data collected from herbarium specimens. The outcomes include plastid genome assemblies, mitochondrial genome assemblies, nuclear 35S ribosomal DNAs (NTS+ETS+18S+ITS1+5.8S+ITS2+25S), alignments of gene and intergenic regions, and a species tree. On the other hand, herbarium specimens provides unparalleled opportunity to collect molecular data, morphology, and distribution with nearly complete taxon sampling. Thus this approach has tremendous potential in taxonomy, biogeography, and macroevolution.

Using the Barbados Cherry family Malpighiaceae as an example, our genome skimming approach provides ample data to resolve both the recent radiation in *Bunchosia* (>135 sp within 10 Myr) as well as the early divergence history of the family in the late Cretaceous. 

Citation: ###

## General guidelines for genome skimming data collection

## Step 1: Filter adapters with blastn
Choose a taxonID for each data set. This taxonID will be used throughout the analysis. Use short taxonIDs with no special characters.

Remove adapters for paired end reads:

	python filter_fastq.py taxonID_1.fq taxonID_2.fq adapter_file num_cores

Alternatively, if sequences are single end reads:
	
	python filter_fastq.py taxonID.fq adapter_file num_cores

The output files are taxonID_1.fq.filtered and taxonID_2.fq.filtered for paired end reads, and taxonID.fq.filtered for single end reads. I use filter_fastq.py only for sequences downloaded from SRA with adapter sequences and phred score offset unknown. For data sets with known adapter and phred scores I use Trimmomatic instead (see below). Trimmomatic is much faster but would not give an error message if you the adapter or phred score offset were misspedified.

##Step 2: 

