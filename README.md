# PhyloHerb
**Phylo**genomic Analysis Pipeline for **Herb**arium Specimens

This bioinformatic tutorial provides detailed guidance to process **genome skimming** data collected from herbarium specimens. The outcomes include plastid genome assemblies, mitochondrial genome assemblies, nuclear 35S ribosomal DNAs (NTS+ETS+18S+ITS1+5.8S+ITS2+25S), alignments of gene and intergenic regions, and a species tree. Combined with the morphological and distribution data from herbarium specimens, this approach provides an unparalleled opportunity to study **taxonomy, biogeography, and macroevolution with nearly complete taxon sampling**.

We have tested this pipeline in the Barbados Cherry family Malpighiaceae, Clusiaceae, and several groups of algae. Each of these datasets contains hundreds to thousands of species and our pipeline extracts ample data to resolve both recent radiations (e.g., *Bunchosia*, Malpighiaceae >135 sp within 10 Myr) and ancient divergences (e.g., the divergence of red algea at hundreds of millions of years ago). 

## I. Prerequisites
To process large datasets (>20 sp), high performance cluster is recommended. Mac and PC can suffer from insufficient memory during the assembly, alignment, and phylogenetic reconstruction. Installation instructions for some of the following software can be found [here](/botany2021_tutorial/README.md).

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
8. Manual assembly examination: [Geneious](https://www.geneious.com/) (the licensed version are required)

### Phylogeny
9. [IQ-TREE](http://www.iqtree.org/) or [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) or [ExaML](https://cme.h-its.org/exelixis/web/software/examl/index.html)

## II. General guidelines for genome skimming data collection

**For the impatient:**

If interested in phylogeny alone, up to 384 samples (4 plates * 96 samples/plate) can be multiplexed on a single Illumina HiSeq 2500 lane for most flowering plants. Using the NovaSeq plastform can generate more complete genomes due to its larger output, but currently we cannot put more than 384 multiplexed samples due to the barcode limitation. If circularized plastid genomes are needed, >2 Gb data per species can usually get you there, which translates to ~60 samples per lane.

*IMPORTANT*: If your species have fewer-than-usual plastids per cell or exceptionally large genomes, you need to reduce the number of multiplexed species per sequencing lane. Use the following equation to calculate the expected coverage of plastid genome:

<img src="/images/plastid_perc.png" width="600" height="100">

Minimally, you want the plastid coverage to be larger than 10%.

**FAQ**

1. DNA extraction from herbarium specimens? How?!

We have successfully extracted DNA from 200-year-old specimens. Age matters less than the preservation methods (see [this paper](https://www.frontiersin.org/articles/10.3389/fevo.2019.00439/full)). Standard commercial DNA extraction kits are frequently used to obtain DNA (e.g, Tiangen DNAsecure Plant Kit, Qiagen DNeasy Plant Mini Kit). We used a [Promega Maxwell](https://www.promega.com/products/lab-automation/maxwell-instruments/maxwell-rsc-instrument/?catNum=AS4500) instrument that can process 16 DNA samples simultaneously and extract their DNAs within an hour. This automatic approach is certainly more labour efficient, but manual extractions have more guaranteed yields for delicate precious samples.

2. Where can I find the genome sizes of my species?

In addition to searching through the literature or conduct your own flow cytometry experiments, you could also check the [Plant DNA C-value database](https://cvalues.science.kew.org/) put together by Kew.

3. NGS library preparation and multiplexing

We used the [KAPA HyperPlus Kit](https://sequencing.roche.com/en/products-solutions/products/sample-preparation/dna-reagents/library-preparation/kapa-hyperplus.html) for NGS library. Many institutes provided services for NGS library preparation with robots. We have used quarter reaction (1/4 of all reagents) for our NGS libraries, and it works just fine.

4. Where are the limits?

About 1-3% of the reads from genome skimming are from plastomes. Theoretically this value vary with the size of the nuclear genome and the abundance of plastids within a cell, but we found it to be relatively consistent across flowering plant species despite the dramatic difference in their genome sizes (200 Mb to 3Gb). Below is a **very rough** guidance of what you may expect from certain amount of input data.

<img src="/images/coverage.png" width="400" height="130">

## III. Assembly

We will use [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) to assemble the plastome, mitochondrial genome, and ribosomal regions. It requires minimal tweak for various types of data. I highly recommend [installing it using conda](https://github.com/Kinggerm/GetOrganelle#installation--initialization) so that all its dependencies are in your environment.

### 1. Input:

Illumina FASTQ reads for each species, single-ended or pair-ended, zipped or unzipped. Do not filter the reads or trim adapters, GetOrganelle will take care of it.

### 2. How to:

After loading GetOrganelle to your environment, the basic commands for running assembly with pair end data is as follows:

```
#To assemble plant plastome
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <plastome_output> -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <nr_output> -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 <forward.fq> -2 <reverse.fq> -o <mito_output> -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
```

If you want to use your own reference sequences for assembly, you can provide the seed fasta file by adding `-s <reference.fas>`.

### 3. Large dataset and batch submission to cluster

If you are dealing with large number of species, running them one by one is too tedious. Here, we will submit individual assembly task to the cluster and run them simultaneously. An example bash file is provided in `/phyloherbLib/getorg.sh`. We will also use short and informative output prefix for each species. You can submit your job by typing

```
sbatch getorg.sh <forward.fq> <backward.fq> <output prefix>
```
*IMPORTANT*: Make sure you load the correct environment and provide absolute path to the input data if they are not in the current directory by modifying relavant variables in `getorg.sh`. Instructions for single-end data can also be found in `getorg.sh`.

### 4. Output

The batch submission will generate three subdirectories `chl/`, `ITS/`, and `mito/`, each containing Getorganelle output directories named after sample-specific prefixes.

### 5. Assembly visualization with Bandage


### 6. Assembly QC 

After the assemblies are completed, you can summarize the results using the QC function of phyloherb. For each species, it will extract the following information: the number of total input reads, the number of reads used for assembly, average base coverage, the total length of the assembly, GC%, and whether the assembly is circularized. 

```
python phyloherb.py -a qc -s sample_sheet.tsv -i <directory containing assemblies> -o <output directory>
```
This command will copy all of the assemblies under the input directory to a new directory and rename the files based on their species prefixes. In the output directory, you will also find a summary spreadsheet `assembly_sum.tsv` with the following header:
```
sp_prefix	Total_reads	Reads_in_target_region	Average_base_coverage	Length	GC%	Circularized
```

## VI. Annotation and organelle structure variarion
Annotation is not necessary if you are interested in phylogeny alone, but if you want to submit your circularized assemblies to GenBank or extract intergenic regions from your spcecies, it is a must.

The most convenient tool I have used is the web-based tool [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html). I have concatenated 100 plastomes and annotated them all at once on GeSeq. But if you are annotating hundreds of plastomes, command-line based tools like [PGA](https://github.com/quxiaojian/PGA) might be a better option.

## V. Alignment generation
Phyloherb will identify the best-matching region of each gene/intergenic region in the assemblly using BLAST. We provide a build-in database of plastid genes from 100 angiosperm species. This database is sufficient for getting genes from species that are not too distantly related to our reference species (species list [here](/phyloherbLib/reference_sp.list)). But you can also supply your own reference in a fasta file following the instructions below.

1. Generate reference gene sequences.
You can obtain your gene references from GenBank or the alignment files generated by GeSeq above. All references can be included in a single fasta file `gene_ref.fas`. The headers should begin by gene names, an underscore, and then a species name (which will be ignored):

```
>gene1_sp1
atcg...
>gene1_sp2
atcg...
>gene2_sp1
atcg...
>gene2_sp2
```

If you are using our build-in plastid gene database, a list of the genes in the database can be found [here](/phyloherbLib/gene.list). You can use a subset of these genes by supplying a `gene_list.txt` file in the next step.

2. Extract
 

