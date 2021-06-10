from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import os, argparse, sys, gzip, shutil
from Bio import SeqIO
import phyloherbLib
from ete3 import Tree

parser = argparse.ArgumentParser(description='PhyloHerb is a bioinfomatic utility wrappepr to process genome skimming data for phylogenomics studies.')
parser.add_argument('-a', help='execution mode, options include[submision, qc, ortho, ]', required=True)
parser.add_argument('--suffix',  help='suffix of alignment files', required=True)
parser.add_argument('--output',  help='output directory', required=True)
parser.add_argument('--loci_order',  help='(optional) a file containing the order of the loci in the concatenation')

args = parser.parse_args()



def submiter_gen(bash_file,sample_sheet,output):
	sp_sheet=open(sample_sheet).readlines()
	out=open(output,'a')
	for l in sp_sheet[1:]:
		out.write('sbatch '+bash_file+' '+l.split()[1]+' '+l.split()[2]+' '+l.split()[0]+'\n')
	out.close()


def qc(sample_sheet,input_dir,output_dir):
	sp_sheet=open(sample_sheet).readlines()
	sp_sheet=[i.split()[0] for i in sp_sheet[1:]]
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	out=open(output_dir+'/assembly_sum.tsv','a')
	out.write('\t'.join(['sp_prefix','Total_reads','Reads_in_target_region','Average_base_coverage','Length','GC%','Circularized'])+'\n')
	for sp in sp_sheet:
		try:
			circ='No'
			#extract some information from log file
			log=open(input_dir+'/'+sp+'/get_org.log.txt').readlines()
			for l in log:
				if 'Reads used =' in l:
					two_reads=l.split()[-1]
					total_reads=sum([int(i) for i in two_reads.split('+')])
				elif 'base-coverage =' in l:
					base_cov=l.split()[-1]
				elif 'Result status' in l:
					if l.split(': ')[-1]=='circular genome':circ='Yes'
			#get number of reads in target reagion
			target_reads=open(input_dir+'/'+sp+'/seed/embplant_pt.initial.fq').readlines()
			target_reads=len(target_reads)/4
			#get info from the assembly
			assem=os.listdir(input_dir+'/'+sp)
			assem=[i for i in assem if i.endswith('path_sequence.fasta')]
			shutil.copy(input_dir+'/'+sp+'/'+assem[0], output_dir+'/'+sp+'.assembly.fas')
			assem_seq=open(input_dir+'/'+sp+'/'+assem[0]).readlines()
			assem_len=0
			GC=0
			for l in assem_seq:
				if not l.startswith('>'):
					assem_len=assem_len+len(l)
					GC=GC+l.count('G')+l.count('C')+l.count('g')+l.count('c')
			out.write('\t'.join([sp,str(total_reads),str(target_reads),base_cov,str(assem_len),str(float(GC)/assem_len),circ])+'\n')
		except IOError:
			try:
				log=open(input_dir+'/'+sp+'/get_org.log.txt').readlines()
				for l in log:
					if 'Reads used =' in l:
						two_reads=l.split()[-1]
						total_reads=sum([int(i) for i in two_reads.split('+')])
				target_reads=open(input_dir+'/'+sp+'/seed/embplant_pt.initial.fq').readlines()
				target_reads=len(target_reads)/4
				out.write('\t'.join([sp,str(total_reads),str(target_reads),'NA','NA','NA','NA'])+'\n')
			except IOError:
				out.write('\t'.join([sp,'NA','NA','NA','NA','NA','NA'])+'\n')
	

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

