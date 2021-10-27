# PhyloHerb		<img src="/images/logo.png" width="80" height="80">
**Phylo**genomic Analysis Pipeline for **Herb**arium Specimens

PhyloHerb is a wrapper program to process **genome skimming** data collected from herbarium specimens. The outcomes include the plastid genome (plastome) assemblies, mitochondrial genome assemblies, nuclear 35S ribosomal DNAs (NTS+ETS+18S+ITS1+5.8S+ITS2+28S), alignments of gene and intergenic regions, and a species tree. Combined with the morphological and distribution data from herbarium specimens, this approach provides an unparalleled opportunity to study **taxonomy, biogeography, and macroevolution with nearly complete taxon sampling**.

We have tested this pipeline in the Barbados Cherry family Malpighiaceae, Clusiaceae, and several groups of algae. Each of these datasets contains hundreds to thousands of species and our pipeline extracts ample data to resolve both recent radiations (e.g., *Bunchosia*, Malpighiaceae >135 sp within 10 Myr) and ancient divergences (e.g., the divergence of red algea hundreds of millions of years ago). 

**License**: GNU General Public License

**Citation**: TBD; PhyloHerb relies on various dependancies, please cite these works as well.

GetOrganelle: Jin, Jian-Jun, Wen-Bin Yu, Jun-Bo Yang, Yu Song, Claude W. Depamphilis, Ting-Shuang Yi, and De-Zhu Li. "GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes." Genome biology 21, no. 1 (2020): 1-31. 

Bandage: Wick, Ryan R., Mark B. Schultz, Justin Zobel, and Kathryn E. Holt. "Bandage: interactive visualization of de novo genome assemblies." Bioinformatics 31, no. 20 (2015): 3350-3352.

pasta: Mirarab, Siavash, Nam Nguyen, Sheng Guo, Li-San Wang, Junhyong Kim, and Tandy Warnow. "PASTA: ultra-large multiple sequence alignment for nucleotide and amino-acid sequences." Journal of Computational Biology 22, no. 5 (2015): 377-386.

MAFFT: Katoh, Kazutaka, and Daron M. Standley. "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." Molecular biology and evolution 30, no. 4 (2013): 772-780.

## I. Prerequisites and installation

To process large datasets (>20 sp), high performance cluster is recommended. Mac and PC can suffer from insufficient memory during the assembly, alignment, or phylogenetic reconstruction. Installation instructions for some of the following software can be found [here](/botany2021_tutorial/README.md).

### PhyloHerb
PhyloHerb is designed to extract orthologous genetic regions from preexisting assemblies. It also provides utility functions to assist assembly and phylogenomic data curation. 

**Installation steps:**

1. Install [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): Make sure BLAST is callable in your current environment.

2. Install two python modules: [Biopython](https://biopython.org/) and [ete3](http://etetoolkit.org/). They can be most easily installed using `conda`.

3. To install PhyloHerb, simply download it use `git` or decompressing the zip file:

```
git clone https://github.com/lmcai/PhyloHerb.git
```

*IMPORTANT*: PhyloHerb is currently only compatible with **Python 3**.

To update your local version for any future updates, `cd` into the `PhyloHerb` directory then type
```
git fetch --prune origin
git reset --hard origin
git clean -f -d
```

### Optional programs
1. Assembler: [GetOrganelle v1.7.0+](https://github.com/Kinggerm/GetOrganelle). Most easily installed via `conda`. All of the dependancies will be installed automatically.
2. Assembly viewer: [Bandage](https://rrwick.github.io/Bandage/)
3. Optional assembly viewer: [Geneious](https://www.geneious.com/) (the complementary version is sufficient)
4. Aligner: 

	[Pasta](https://github.com/smirarab/pasta) for highly variable regions such as the ITS sequences
	
	[MAFFT](https://mafft.cbrc.jp/alignment/software/) for less variable regions or long alignments (>5 kb) that pasta may not be able to handle when the number of species is high (>500 sp)
5. Manual assembly examination: [Geneious](https://www.geneious.com/) (the licensed version are required)

6. Phylogeny: [IQ-TREE](http://www.iqtree.org/) or [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) or [ExaML](https://cme.h-its.org/exelixis/web/software/examl/index.html)


## II. Quick start

<img src="/images/PH_pipeline.png" width="674" height="335">

**Core functions of PhyloHerb**

### Usage:

```
phyloherb.py [-h] -m mode [-i dir] [-o dir] [-suffix string] [-sp file]
                    [-g file] [-l integer] [-n integer] [-ref file] [-mito]
                    [-rdna] [-t file] [-missing float 0-1] [-f mode]
                    [-gene_def file] [-b file] [-s file]
```

### Options:
```
	-i dir			input directory
	-o dir			output directory
	-m str			execution mode, choose one of the following: qc, ortho, conc, order, getseq, submission
	-suffix str		[qc, ortho, conc mode] suffix of input files
	-sp	file		[ortho mode] a file containing a list of species
	-g file			[ortho and conc mode] a file containing a list of loci
	-l int			[ortho mode] minimum length of blast hits
	-n int			[ortho mode] number of threads for BLAST
	-ref file		[ortho mode] custom reference sequences
  	-mito			[ortho mode] extract mitochondrial genes using build-in references
  	-rdna           	[ortho mode] extract nuclear ribosomal regions using build-in references
  	-t file				[order mode] newick tree file to order alignments based on phylogeny
  	-missing float 0-1  [order mode] maximum proportion of missing data allowed for each species
  	-f mode				[getseq mode] choose one of the following: gene, genetic_block, intergenic
  	-gene_def file      [getseq mode] a gene delimitation file that defines genetic blocks
  	-b file             [submission mode] path to the bash file
  	-s file             [submission mode] path to the taxon sampling sheet
```

**1. Ortholog gene extraction using built-in database**

**A.** If you have your assemblies and want to extract genes using our curated database:

*Input:* Place all fasta formated assemblies in one folder. Make sure they have consistent suffix. To extract genes, use the following command:

```
#For plastid 
python phyloherb.py -m ortho -i <input directory> -o <output directory> -suffix <suffix>
#For rRNA
python phyloherb.py -m ortho -i <input directory> -o <output directory> -suffix <suffix> -rdna
#For mitochondrion
python phyloherb.py -m ortho -i <input directory> -o <output directory> -suffix <suffix> -mito
```
List of genes and species included in our built-in database can be found [here](https://github.com/lmcai/PhyloHerb/tree/main/database).

*Output:* In the output folder, you can find fasta sequences named after genes. The header within each fasta is consistent with the species names (file names of the input assemblies).

**B.** If you need to assemble organelle genomes and rRNA regions, see Section ## below.

**2. Alignment and concatenation**

For small dataset in conserved regions, you can use `MAFFT` for align. For large dataset across distantly related species, we recommend `PASTA`. An example bash file to run MAFFT and PASTA can be found in [mafft_pasta.sh](phyloherbLib/mafft_pasta.sh). To concatenate sequences, place all fasta alignments in one folder and use the following command:

```
python phyloherb.py -m conc -i < input directory> -o <output prefix> -suffix <alignment suffix>
```
*Output:* `*.conc.fas` and `*.conc.nex` are concatenated sequences in fasta and nexus formats, respectively. `*.partition` is the gene partition file that can be used for [PartitionFinder](https://www.robertlanfear.com/partitionfinder/).

**3. Reorder alignments based on phylogeny to assist manual curation**

Reorder sequences based on phylogeny can help distinguish analytical errors versus shared mutations. 

*Input:* One species tree in newick format; multiple FASTA alignments in one folder.

```
python phyloherb.py -m order -t <tree> -i <input dir> -o <output dir> -suffix <alignment suffix>
```
*Output:* `*.ordered.fas` is the reordered alignment for each gene. `*.pasta_ref.tre` is the pruned species tree that can be used as the reference tree in the second PASTA alignment.

## III. General guidelines for genome skimming data collection

**Overview**

If interested in phylogeny alone, up to 384 samples (4 plates * 96 samples/plate) can be multiplexed on a single Illumina HiSeq 2500 lane for most flowering plants. Using the NovaSeq plastform can generate more complete genomes thanks to its larger output, but currently we cannot put more than 384 multiplexed samples due to the barcode limitation. If circularized plastid genomes are needed, >2 Gb data per species can usually get you there, which translates to ~60 samples per lane.

*IMPORTANT*: If your species have fewer-than-usual plastids per cell or exceptionally large genomes, you need to reduce the number of multiplexed species per sequencing lane. Use the following equation to calculate the expected base coverage of plastid genome:

<img src="/images/plastid_perc.png" width="600" height="100">

Minimally, you want the plastid coverage to be larger than 10X.

**FAQ**

1. DNA extraction from herbarium specimens? How?!

We have successfully extracted DNAs from 200-year-old specimens. Age matters less than the preservation methods (see [this paper](https://www.frontiersin.org/articles/10.3389/fevo.2019.00439/full)). Standard commercial DNA extraction kits are frequently used to obtain DNA (e.g, Tiangen DNAsecure Plant Kit, Qiagen DNeasy Plant Mini Kit). We used a [Promega Maxwell](https://www.promega.com/products/lab-automation/maxwell-instruments/maxwell-rsc-instrument/?catNum=AS4500) instrument that can extract DNA from 16 samples simultaneously within an hour. This automatic approach is certainly more labour efficient, but manual extractions have more guaranteed yields for delicate samples.

2. Where can I find the genome sizes of my species?

In addition to searching through the literature or conducting your own flow cytometry experiments, you could also check the [Plant DNA C-value database](https://cvalues.science.kew.org/) organized by Kew.

3. NGS library preparation and multiplexing

We used the [KAPA HyperPlus Kit](https://sequencing.roche.com/en/products-solutions/products/sample-preparation/dna-reagents/library-preparation/kapa-hyperplus.html) for NGS library. Many institutes provided services for NGS library preparation with robots. We have used quarter reaction (1/4 of all reagents) for our NGS libraries, and it works just fine.

4. Where are the limits?

About 0.5-6% of the reads from genome skimming come from plastomes. The base coverage is roughly half for mitochondria and 2X for nuclear ribosomal regions compared to plastids. Theoretically the base coverage vary with the size of the nuclear genome and the abundance of these genetic regions within a cell, but we found it to be relatively consistent across flowering plant species despite the dramatic difference in their genome sizes (200 Mb to 3Gb). Below is a **very rough** estimation of what you may expect for plastome from certain amount of input data.

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

If you are working with a high performance computing cluster with **slurm workload manager**, you can submit individual assembly tasks to the cluster and run them simultaneously. An example bash file is provided in `/phyloherbLib/getorg.sh`. You can submit your job by typing

```
sbatch getorg.sh <forward.fq> <backward.fq> <output prefix>
```

For more details using this bash file, see the [botany2021 tutorial](/botany2021_tutorial/README.md).

*IMPORTANT*: Make sure you load the correct environment and provide absolute path to the input data if they are not in the current directory by modifying relavant variables in `getorg.sh`. Instructions for single-end data can also be found in `getorg.sh`.

### 4. Output

The key output files from Getorganelle include
```
*.path_sequence.fasta, the assembly in fasta format 

extended_K*.assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg, a simplified assembly graph

extended_K*.assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv, a tab-format contig label file for bandage visualization

get_org.log.txt, the log file
```

### 5. Assembly visualization with Bandage

[Bandage](https://rrwick.github.io/Bandage/) is a program for visualising de novo assembly graphs. The assembly graph files *.fastg and *.gfa generated from GetOrganelle be visualized in Bandage and exported into sequences. 

The *.path_sequence.fasta files do not always navigate the right paths for organelle genomes, especially the ones with complicated structures. The authors of GetOrganelle put together a nice [video](https://www.youtube.com/watch?v=cXUV7k-F26w)  introducing how to generate complete (if possible) and accurate sequences from Bandage with different examples.


### 6. Assembly QC 

After the assemblies are completed, you can summarize the results using the `qc` function of phyloherb. For each species, it will extract the following information: the number of total input reads, the number of reads mapped to the target region, average base coverage, the total length of the assembly, GC%, and whether the assembly is circularized. 

```
python phyloherb.py -m qc -s sample_sheet.tsv -i <directory with Getorganelle output folders> -o <output directory>
```
This command will copy all of the assemblies under the input directory to a new directory and rename the files based on their species prefixes. In the output directory, you will also find a summary spreadsheet `assembly_sum.tsv` with the following header:
```
sp_prefix	Total_reads	Reads_in_target_region	Average_base_coverage	Length	GC%	Circularized
```

## IV. Annotation and organellar structure variations

Annotation is **not necessary** if you are interested in phylogeny alone, but if you want to submit your circularized assemblies to GenBank or extract intergenic regions from your spcecies, it is a must.

The most convenient tool I have used is the web-based tool [GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html). I have concatenated 100 plastomes into a single fasta file and annotated them all at once on GeSeq. But if you are annotating hundreds of plastomes, the command-line based tool [PGA](https://github.com/quxiaojian/PGA) might be a better option.

## V. Ortholog identification and alignment

Phyloherb will identify the best-matching region of each gene/intergenic region in the assemblly using BLAST. We provide a build-in database of plastid genes from 355 seed plant families. This database is sufficient for getting genes from species that are not too distantly related to our reference species. You can also supply your own reference in a fasta file following the instructions below. The list of reference species is [here](/phyloherbLib/reference_sp.list). The list of the genes in the database is [here](/phyloherbLib/gene.list). 

1. Generate reference gene sequences **(optional)**

You can obtain your gene references from GenBank or the alignment files generated by GeSeq above. All references can be included in a single fasta file `gene_ref.fas`. It is recommended that for diverse groups, you use **more than one** reference sequences for each gene to catch the genetic variance. 

In `gene_ref.fas`, the headers should begin by the gene name then followed by an underscore and a species name (which will be ignored):

```
>gene1_sp1
atcg...
>gene1_sp2
atcg...
>gene2_sp1
atcg...
>gene2_sp2
```

*IMPORTANT*: If you have annotated plastid assemblies (e.g., from GeqSeq) in genbank format, you can use the `getseq` function of PhyloHerb to obtain gene and intergenic regions (see section V.4 below).

2. Extract orthologous gene or intergenic regions from the assembly

We can  extract the target gene regions using the `ortho` function of phyloherb. This function will conduct BLAST search in the assembly, extract the best matching regions, and output them to a directory.

In the output directory, orthologous genes will be written to separate fasta files and the headers will be the species prefixes.

```
python phyloherb.py -m ortho -i <directory containing assemblies> -o <output directory>
```
You can choose to extract a subset of genes from a subset of the species by supplying a `-g gene_subset.txt` and `-sp species_subset.txt`. Example files can be found in [gene_subset.txt](/example/gene_subset.txt) and [species_subset.txt](/example/species_subset.txt). You can also set a minimum length limit for gene region extraction via `-l <lower limit>`. Blast hits shorter than this will not be use.
```
python phyloherb.py -m ortho -i <directory containing assemblies> -o <output directory> -g <gene list> -sp <species list> -l <length limit>
```

3. Alignment

I like to use the `--adjustdirection` function from `MAFFT` to correct reverse complimentary sequences. Then I will use `pasta` to more accurately align high variable sequences such as the intergenic regions and the ITS regions. `pasta` first generates a guidance tree, then align among closely-related species, finally merge the alignments to produce the output.

This is a potentially time consuming step so I recommend running it on the cluster using the example batch file [mafft_pasta.sh](phyloherbLib/mafft_pasta.sh).

Copy `mafft_pasta.sh` to the same directory where the gene sequences are located. Modify the file to include appropriate environmental parameters. Then the batch job can be submitted to the cluster by typing
```
sbatch mafft_pasta.sh <gene_1>
sbatch mafft_pasta.sh <gene_2>
```

4. Intergenic regions

Intergenic regions are mostly useful for phylogenetic research among closely related species. In addition, some plastid genes are very short, we can combine several adjacent gene and intergenic region to produce a longer genetic block for downstream analyses. Combining short genetic blocks into long ones has the benefit of higher blast accuracy, but you should make sure the target regions **do not have structural variations** among different species. `PhyloHerb` offers three modes for genetic region extraction **from Genbank annotations**. The results can be used as reference sequences for ortholog extraction using the `ortho` function of PhyloHerb.

Let's assume that there are seven genes G1-7 on a scaffold.

<img src="/images/seven_genes_black.png" width="400" height="80">

The `-f gene` mode will extract all annotated gene features from all genbank files in the input directory
```
python phyloherb.py -m getseq -f gene -i <input directory> -suffix <genbank suffix> -o <output directory>
```
<img src="/images/gene.png" width="400" height="90">

Both the `-f genetic_block` and `-f intergenic` mode will take a genetic block definition file supplied by `-gene_def`, and extract corresponding regions. This tab-delimited file have three columns: genetic region name, start gene, and end gene. An example can be found [here](/example/gene_def.txt).
```
name	start	end
INT1	G1	G2
INT2	G3	G7
...
```

The `-f genetic_block` mode will extract the coding regions of the start and end genes as well as everything in between. It is good for combining multiple short genes.
```
python phyloherb.py -m getseq -f genetic_block -i <input directory> -suffix <genbank suffix> -o <output directory> -gene_def <gene definition file>
```

<img src="/images/genetic_block.png" width="400" height="90">

Finally, the `-f intergenic` mode generates similar outcomes, but does not include the genes on both ends. It is good for extracting longer intergenic regions.
```
python phyloherb.py -m getseq -f intergenic -i <input directory> -suffix <genbank suffix> -o <output directory> -gene_def <gene definition file>
```

<img src="/images/intergenic.png" width="400" height="90">


5. Nuclear ribosomal regions

PhyloHerb can extract the coding regions of the rDNA repeat (18S+ITS1+5.8S+ITS2+28S). The highly variable NTS and ETS will be discarded. It will first identify rDNAs (18S, 5.8S, and 28S) in the assembly and then extract the ITS regions in between. The output contains five fasta files `18S.fas, ITS1.fas, 5.8S.fas, ITS2.fas, and 28S.fas` in the output directory. To get rDNA sequences, add the `-rdna` flag under the `ortho` mode. 

```
python phyloherb.py -m ortho -i <directory containing assemblies> -o <output directory> -rdna
```

Below is an illustration of the structure and sequence conservation of the nuclear ribosomal region.

<img src="/images/ITS.png" width="400" height="400">


6. Mitochondrial regions

For most plant groups, mitochondria are not phylogenetically informative because the genes evolve too slowly, but the intergenic regions are highly variable. Moreover, the qualities of mitochondrial genomes are usually not as good as plastomes. So I recommend using mitochondrial genes only. If you want to use our build-in mitochondrial gene sequence database, invoke the `-mito` flag under the `ortho` mode. 

A list of the reference mitochondrial genes can be found [here](/database/mito_gene.list). The reference sequences themselves can be found [here](/database/mito_reference.fas).

```
python phyloherb.py -m ortho -i <directory containing assemblies> -o <output directory> -mito 
```

7. Manual curation in Geneious

At this point, it is recommended that you take a initial look at your alignments. **Initial** means be prepared to complete the "alignment-manual check-phylogeny" cycle for at least two rounds to get publication quality data.

The purpose of the initial check is to remove obvious low-quality sequences. Do not conduct any site-based filtering yet! For example, the two sequences highlighted in red below contain too many SNPs (marked in black). They should be removed.

<img src="/images/Geneious.png" width="600" height="400">

Geneious provides a nice interface to work with alignments. You can view alignment statistics, edit sequences, and concatenate alignments. Alternative automatic tools include [trimAL](http://trimal.cgenomics.org/getting_started_with_trimal_v1.2) and [Seaview](http://doua.prabi.fr/software/seaview).


## VI. Phylogeny reconstruction

1. Concatenation

Many tools are available for concatenating alignments. I recommend the `conc` function of phyloherb or Geneious. I have applied both tools to dataset with 1000 sp x 100 genes. The `conc` function of phyloherb will also output a gene delineation file required by `PartitionFinder`.

To use the `conc` function of phyloherb, use the following commands
```
python phyloherb.py -m conc -i <directory containing alignments> -o <output prefix> -suffix <suffix>
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
python phyloherb.py -m order -t <reference.tre> -i <directory containing alignments> -o <output directory> -suffix <suffix> [optional] -missing <float number 0 to 1>
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
