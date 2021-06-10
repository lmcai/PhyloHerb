from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import fnmatch, os, argparse, sys
from Bio import SeqIO
import phyloherbLib
from ete3 import Tree

parser = argparse.ArgumentParser(description='PhyloHerb is a bioinfomatic utility wrappepr to process genome skimming data for phylogenomics studies.')
parser.add_argument('-a', help='working mode, options include[submision, qc]', required=True)
parser.add_argument('--suffix',  help='suffix of alignment files', required=True)
parser.add_argument('--output',  help='output directory', required=True)
parser.add_argument('--loci_order',  help='(optional) a file containing the order of the loci in the concatenation')

args = parser.parse_args()


def order_aln(sptree,input_dir,suffix,output_dir,max_missing):
	t=Tree(sptree)
	total_taxa=[]
	for leaf in t:
		total_taxa.append(leaf.name)
	genes=os.listdir(input_dir)
	genes=[i for i in genes if i.endswith(suffix)]
	
	for g in genes:
		sp2preserve=[]
		t=Tree(sptree)
		y=SeqIO.parse(input_dir+'/'+g,'fasta')
		out=open(output_dir+'/'+g+'.ordered.fas','a')
		for rec in y:
			missing=float(rec.seq.count('-'))/len(rec.seq)
        	if rec.id in total_taxa and missing<max_missing:
                SeqIO.write(rec,out,'fasta')
                sp2preserve.append(rec.id)
        out.close()
    	t.prune(list(set(total_taxa) & set(sp2preserve))) 
    	t.write(format=1, outfile=output_dir+'/'+g+".pasta_ref.tre")
