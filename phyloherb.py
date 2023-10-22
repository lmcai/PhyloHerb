import sys, textwrap
print(textwrap.dedent("""\
 __                 __        ___  __   __  
|__) |__| \ / |    /  \ |__| |__  |__) |__) 
|    |  |  |  |___ \__/ |  | |___ |  \ |__) 
                                            
"""))
print('############################################################\n\
PhyloHerb v1.1.3\n\
A bioinformatic pipeline for herbariomics based biodiversity research\n')

if sys.version_info.major==2:
	print('You are using Python 2. Please upgrade to Python 3. PhyloHerb quit now...')
	quit()

from Bio import AlignIO
import os, argparse, shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description='PhyloHerb is a bioinfomatic utility wrappepr to process genome skimming data for phylogenomics studies.')
parser.add_argument('-m', metavar='mode', help='execution mode, choose one of the following: qc, ortho, conc, order, getseq, assemb, submission', required=True)
parser.add_argument('-i', metavar='dir', help='input directory')
parser.add_argument('-o', metavar='dir', help='output directory')
parser.add_argument('-suffix', metavar='string', help='[qc, ortho, conc mode] suffix of input files')
parser.add_argument('-sp',  metavar='file', help='[ortho mode] a file containing a list of species')
parser.add_argument('-g',  metavar='file', help='[ortho and conc mode] a file containing a list of loci')
parser.add_argument('-l',  metavar='integer', help='[ortho mode] minimum length of blast hits')
parser.add_argument('-evalue',  metavar='float', help='[ortho mode] evalue threshold for BLAST')
parser.add_argument('-n',  metavar='integer', help='[ortho, assemb mode] number of threads')
parser.add_argument('-ref',  metavar='file', help='[ortho mode] custom reference sequences')
parser.add_argument('-mito',  help='[ortho mode] extract mitochondrial genes using build-in references',action='store_true')
parser.add_argument('-rdna',  help='[ortho mode] extract nuclear ribosomal regions using build-in references',action='store_true')
parser.add_argument('-nuc',  help='[ortho mode] extract low-copy nuclear loci',action='store_true')
parser.add_argument('-r1',  metavar='reads file', help='[assemb mode] forward reads')
parser.add_argument('-r2',  metavar='reads file', help='[assemb mode] reverse reads')
parser.add_argument('-rs',  metavar='reads file', help='[assemb mode] single-end or unpaired reads, separate mulitple files by \',\'')
parser.add_argument('-prefix',  metavar='string', help='[assemb mode] assembly output prefix')
parser.add_argument('-t', metavar='file', help='[order mode] newick tree file to order alignments based on phylogeny')
parser.add_argument('-missing', metavar='float 0-1', help='[order mode] maximum proportion of missing data allowed for each species')
parser.add_argument('-f',  metavar='mode', help='[getseq mode] how to extract loci, choose one of the following: gene, genetic_block, intergenic')
parser.add_argument('-gene_def',  metavar='file', help='[getseq mode] a gene delimitation file that defines genetic blocks')
parser.add_argument('-b', metavar='file', help='[submission mode] path to the bash file')
parser.add_argument('-s',  metavar='file', help='[submission mode] path to the taxon sampling sheet')



args = parser.parse_args()

def submiter_gen(bash_file,sample_sheet,output):
	sp_sheet=open(sample_sheet).readlines()
	out=open(output,'a')
	for l in sp_sheet[1:]:
		d=out.write('sbatch '+bash_file+' '+l.split()[1]+' '+l.split()[2]+' '+l.split()[0]+'\n')
	out.close()

def qc_assemblies_only(sp_list,input_dir,output_dir):
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	out=open(output_dir+'/assembly_sum.tsv','a')
	out.write('\t'.join(['sp_prefix','Length','GC%','Number_scaffolds'])+'\n')
	for sp in sp_list:
		try:
			assem_seq=open(input_dir+'/'+sp).readlines()
			assem_len=0
			GC=0
			scaffolds=0
			for l in assem_seq:
				if not l.startswith('>'):
					assem_len=assem_len+len(l)
					GC=GC+l.count('G')+l.count('C')+l.count('g')+l.count('c')
				else:
					scaffolds=scaffolds+1
			out.write('\t'.join([sp,str(assem_len),str(float(GC)/assem_len),str(scaffolds)])+'\n')
		except:
			out.write('\t'.join([sp,'NA','NA','NA'])+'\n')

def qc(sp_sheet,input_dir,output_dir):
	failed=[]
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	out=open(output_dir+'/assembly_sum.tsv','a')
	out.write('\t'.join(['sp_prefix','Total_reads','Reads_in_target_region','Average_base_coverage','Length','GC%','Number_scaffolds','Circularized'])+'\n')
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
					if l.split(': ')[-1]=='circular genome\n':circ='Yes'
			#get number of reads in target reagion
			target_reads_file=os.listdir(input_dir+'/'+sp+'/seed/')
			target_reads_file=[j for j in target_reads_file if j.endswith('.fq')]
			target_reads=open(input_dir+'/'+sp+'/seed/'+target_reads_file[0]).readlines()
			target_reads=len(target_reads)/4
			#get info from the assembly
			assem=os.listdir(input_dir+'/'+sp)
			assem=[i for i in assem if i.endswith('path_sequence.fasta')]
			shutil.copy(input_dir+'/'+sp+'/'+assem[0], output_dir+'/'+sp+'.assembly.fas')
			assem_seq=open(input_dir+'/'+sp+'/'+assem[0]).readlines()
			assem_len=0
			GC=0
			scaffolds=0
			for l in assem_seq:
				if not l.startswith('>'):
					assem_len=assem_len+len(l)
					GC=GC+l.count('G')+l.count('C')+l.count('g')+l.count('c')
				else:
					scaffolds=scaffolds+1
			out.write('\t'.join([sp,str(total_reads),str(target_reads),base_cov,str(assem_len),str(float(GC)/assem_len),str(scaffolds),circ])+'\n')
		except (IOError, IndexError) as e:
			try:
				log=open(input_dir+'/'+sp+'/get_org.log.txt').readlines()
				for l in log:
					if 'Reads used =' in l:
						two_reads=l.split()[-1]
						total_reads=sum([int(i) for i in two_reads.split('+')])
				target_reads_file=os.listdir(input_dir+'/'+sp+'/seed/')
				target_reads_file=[j for j in target_reads_file if j.endswith('.fq')]
				target_reads=open(input_dir+'/'+sp+'/seed/'+target_reads_file[0]).readlines()
				target_reads=len(target_reads)/4
				out.write('\t'.join([sp,str(total_reads),str(target_reads),'NA','NA','NA','NA'])+'\n')
			except IOError:
					out.write('\t'.join([sp,'NA','NA','NA','NA','NA','NA'])+'\n')
					failed.append(sp)
	if len(failed)>0:
		if len(failed)==len(sp_sheet):
			print('No GetOrganelle outputs found. If you are using assembly fasta files only, be sure to add the -suffix argument.')
		else:
			print('Cannot find GetOrganelle outputs in the directory '+input_dir+' for the following species: '+', '.join(failed))
	
	
def ortho_extraction(sp,reference_seq,input_dir,output_dir,genes,min_len,threads,evalue):
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	print('processing species '+sp)
	lib_ID=sp
	S= 'makeblastdb -in ' +input_dir+'/'+ lib_ID +' -dbtype nucl -out '+lib_ID+' &>/dev/null'
	os.system(S)
	S = 'blastn -task dc-megablast -db '+lib_ID+' -query ' + reference_seq + ' -num_threads '+threads+' -outfmt 6 -evalue '+evalue+' -out '+ lib_ID +'.blast.out &>/dev/null'
	os.system(S)
	x=open(lib_ID+'.blast.out').readlines()
	y=SeqIO.index(input_dir+'/'+lib_ID,'fasta')
	a={}
	for g in genes:
		#print g
		best=0
		a[g]=[r for r in x if g+'_' in r or g+'\t' in r]
		min_evalue=1
		length=1
		for rec in a[g]:
			#print rec
			if float(rec.split('\t')[10])<=min_evalue and float(rec.split('\t')[3])>length:
				min_evalue=float(rec.split('\t')[10])
				length=float(rec.split('\t')[3])
				best=rec
		#best is the best match
		try:
			hit=best.split('\t')[1]
			start=min(int(best.split('\t')[8]),int(best.split('\t')[9]))
			end=max(int(best.split('\t')[8]),int(best.split('\t')[9]))
			if end-start>min_len:
				seq=y[hit].seq[(start-1):(end-1)]
				SeqIO.write(SeqRecord(seq,lib_ID.split('.')[0], '', ''),open(output_dir+'/'+g+'.fas','a'),'fasta')
		except (NameError,AttributeError):continue
	os.remove(lib_ID+'.nhr')
	os.remove(lib_ID+'.nin')
	os.remove(lib_ID+'.nsq')
	
def get_ITS(sp,blast_file,input_dir,output_dir,min_len):
	ITS1=[]
	ITS2=[]
	blast_txt=open(blast_file).readlines()
	y=SeqIO.index(input_dir+'/'+sp,'fasta')
	a={}
	for g in ['18S','5.8S','28S']:
		#print g
		best=0
		a[g]=[r for r in blast_txt if g+'_' in r]
		min_evalue=1
		length=1
		for rec in a[g]:
			#print rec
			if float(rec.split('\t')[10])<=min_evalue and float(rec.split('\t')[3])>length:
				min_evalue=float(rec.split('\t')[10])
				length=float(rec.split('\t')[3])
				best=rec
		a[g]=best
	#extract ITS1
	try:
		if a['18S'].split('\t')[1] == a['5.8S'].split('\t')[1]:
			ITS1=[int(a['18S'].split('\t')[8]),int(a['18S'].split('\t')[9]),int(a['5.8S'].split('\t')[8]),int(a['5.8S'].split('\t')[9])]
			ITS1.sort()
			start=ITS1[1]
			end=ITS1[2]
			seq=y[a['18S'].split('\t')[1]].seq[start:(end-1)]
			output_handle=open(output_dir+'/ITS1.fas','a')
			d=output_handle.write(">%s\n%s\n" % (sp.split('.')[0],seq))
			output_handle.close()
	except:pass
	try:
		if a['5.8S'].split('\t')[1] == a['28S'].split('\t')[1]:
			ITS2=[int(a['5.8S'].split('\t')[8]),int(a['5.8S'].split('\t')[9]),int(a['28S'].split('\t')[8]),int(a['28S'].split('\t')[9])]
			ITS2.sort()
			start=ITS2[1]
			end=ITS2[2]
			seq=y[a['5.8S'].split('\t')[1]].seq[start:(end-1)]
			output_handle=open(output_dir+'/ITS2.fas','a')
			d=output_handle.write(">%s\n%s\n" % (sp.split('.')[0],seq))
			output_handle.close()
	except:pass


def order_aln(sptree,input_dir,suffix,output_dir,max_missing):
	t=Tree(sptree)
	total_taxa=[]
	for leaf in t:
		total_taxa.append(leaf.name)
	genes=os.listdir(input_dir)
	genes=[i for i in genes if i.endswith(suffix)]
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	for g in genes:
		sp2preserve=[]
		t=Tree(sptree)
		y=SeqIO.index(input_dir+'/'+g,'fasta')
		out=open(output_dir+'/'+'.'.join(g.split('.')[:-1])+'.ordered.fas','a')
		for sp in total_taxa:
			try:
				missing=float(y[sp].seq.count('-')+y[sp].seq.count('N'))/len(y[sp].seq)
				if missing<max_missing:
					SeqIO.write(y[sp],out,'fasta')
					sp2preserve.append(sp)
			except KeyError:pass
		out.close()
		t.prune(list(set(total_taxa) & set(sp2preserve))) 
		t.write(format=1, outfile=output_dir+'/'+'.'.join(g.split('.')[:-1])+".pasta_ref.tre")

def concatenation(input_dir,files,output):
    total_sp=[]
    for fn in files:
        seq_name=open(input_dir+'/'+fn).readlines()
        seq_name=[j[1:].strip() for j in seq_name if j.startswith('>')]
        total_sp=total_sp+seq_name
    total_sp=list(set(total_sp))
    seq={}
    for i in total_sp:
        seq[i]=''
    b=1
    out1=open(output+'.conc.fas','a')
    out2=open(output+'.partition','a')
    for fn in files:
        recs=SeqIO.index(input_dir+'/'+fn,'fasta')
        for rec in recs:
            seq_len=len(recs[rec].seq)
            break
        for j in total_sp:
            try:seq[j]=seq[j]+str(recs[j].seq)
            except KeyError:seq[j]=seq[j]+'-'*seq_len
        out2.write(fn + '='+str(b)+'-'+str(b+seq_len-1)+';\n')
        b=b+seq_len
    for sp in total_sp:
        out1.write('>'+sp+'\n'+seq[sp]+'\n')


def gene_extra(input_dir,suffix,output_dir):
	filenames=os.listdir(input_dir)
	filenames=[j for j in filenames if j.endswith(suffix)]
	if len(filenames)==0:
		print('No file ends with '+suffix+' found in the input directory. PhyloHerb quit...')
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	for f in filenames:
		gb_recs=SeqIO.read(input_dir+'/'+f,'genbank')
		for feature in gb_recs.features:
			if feature.type=='gene':
				try:
					output_handle=open(output_dir+'/'+feature.qualifiers['gene'][0]+'.gene.fas','a')
					d=output_handle.write(">%s\n%s\n" % (feature.qualifiers['gene'][0]+'_'+f.split('.')[0],feature.extract(gb_recs).seq))
					output_handle.close()
				except KeyError:pass

def geneblock_extra(input_dir,suffix,output_dir,gene_def):
	filenames=os.listdir(input_dir)
	filenames=[j for j in filenames if j.endswith(suffix)]
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	genes=open(gene_def).readlines()
	for f in filenames:
		try:
			gb_recs=SeqIO.read(input_dir+'/'+f,'genbank')
		except ValueError:
			print('############################################################\n\
#ERROR:At least one of the files with suffix '+suffix+' is not in GenBank format!\n\
Make sure you put all your genbank reference in the input directory with correct suffix.\n')
			exit()
		#get all gene positions. for genes in the IR region, this will be the position in the second IR.
		gene_start={}
		gene_end={}
		for feature in gb_recs.features:
			if feature.type=='gene':
				try:
					gene_name=feature.qualifiers['gene'][0]
					gene_start[gene_name]=int(feature.location.start)
					gene_end[gene_name]=int(feature.location.end)
				except KeyError:pass
		for l in genes[1:]:
			loci=l.split()[0]
			start_g=l.split()[1]
			end_g=l.split()[2]
			two_gene_pos=[]
			try:
				two_gene_pos=[gene_start[start_g],gene_end[start_g],gene_start[end_g],gene_end[end_g]]
				two_gene_pos.sort()
				start=min(two_gene_pos)
				end=max(two_gene_pos)
				seq=gb_recs.seq[(start-1):end]
				output_handle=open(output_dir+'/'+loci+'.geneblock.fas','a')
				d=output_handle.write(">%s\n%s\n" % (loci+'_'+f.split('.')[0],seq))
				output_handle.close()
			except KeyError:
				print('Cannot find the following gene: '+l.strip()+' in the GenBank file: '+f)

def intergenic_extra(input_dir,suffix,output_dir,gene_def):
	filenames=os.listdir(input_dir)
	filenames=[j for j in filenames if j.endswith(suffix)]
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	genes=open(gene_def).readlines()
	for f in filenames:
		gb_recs=SeqIO.read(input_dir+'/'+f,'genbank')
		#get all gene positions. for genes in the IR region, this will be the position in the second IR.
		gene_start={}
		gene_end={}
		for feature in gb_recs.features:
			if feature.type=='gene':
				try:
					gene_name=feature.qualifiers['gene'][0]
					gene_start[gene_name]=int(feature.location.start)
					gene_end[gene_name]=int(feature.location.end)
				except KeyError:pass
		#get intergenic regions and write to file
		for l in genes[1:]:
			loci=l.split()[0]
			start_g=l.split()[1]
			end_g=l.split()[2]
			two_gene_pos=[]
			try:
				two_gene_pos=[gene_start[start_g],gene_end[start_g],gene_start[end_g],gene_end[end_g]]
				two_gene_pos.sort()
				seq=gb_recs.seq[(two_gene_pos[1]-1):two_gene_pos[2]]
				output_handle=open(output_dir+'/'+loci+'.intergenic.fas','a')
				d=output_handle.write(">%s\n%s\n" % (loci+'_'+f.split('.')[0],seq))
				output_handle.close()
			except KeyError:
				print('Cannot find the following gene: '+l.strip()+' in the GenBank file: '+f)

def assembly(r1,r2,rs,threads,prefix,reference):
	reference_ID=reference.split('/')[-1]
	reference_ID=reference_ID.split('.')[0]
	#bowtie
	#if os.path.exists(reference_ID+'.rev.1.bt2'):
	#	print('Bowtie database found.')
	if True:
		print('Building Bowtie database using '+ reference)
		bowtieCommand='bowtie2-build '+reference+' '+reference_ID+' > /dev/null'
		print('Bowtie command: '+bowtieCommand)
		os.system(bowtieCommand)
		print('Bowtie database completed.')
	#reads mapping
	print('Mapping reads using '+threads+' threads. This may take a while...')
	if r1!='NA' and rs=='NA':
		#PE only
		bowtieCommand='bowtie2 --very-sensitive-local -p '+threads+' -x '+reference_ID+' -1 '+r1+' -2 '+r2+' -S '+prefix+'.'+reference_ID+'.sam'
	elif r1!='NA' and rs!='NA':
		#PE+SE
		bowtieCommand='bowtie2 --very-sensitive-local -p '+threads+' -x '+reference_ID+' -1 '+r1+' -2 '+r2+' -U '+rs+' -S '+prefix+'.'+reference_ID+'.sam'
	elif r1=='NA':
		#SE
		bowtieCommand='bowtie2 --very-sensitive-local -p '+threads+' -x '+reference_ID+' -U '+rs+' -S '+prefix+'.'+reference_ID+'.sam'
	print('Bowtie command: '+bowtieCommand)	
	os.system(bowtieCommand)
	print('Bowtie reads mapping completed.')
	#get mapped reads
	samCommand='samtools bam2fq -F 4 -1 '+prefix+'.filtered.R1.fq -2 '+prefix+'.filtered.R2.fq -s '+prefix+'.filtered.unpaired.fq '+prefix+'.'+reference_ID+'.sam > /dev/null'
	print('Extracting mapped reads.\nSamtools command: '+samCommand)
	os.system(samCommand)
	print('Reads extraction completed.\nStart spades assembly.')
	spadesCommand='spades.py -t '+threads+' --phred-offset 33 -1 '+prefix+'.filtered.R1.fq -2 '+prefix+'.filtered.R2.fq -s '+prefix+'.filtered.unpaired.fq -k 21,55,85,115 -o '+prefix+'_spades'
	print('Assembling reads using '+threads+' threads in Spades. This might take a while.\nSpades command: '+spadesCommand)
	os.system(spadesCommand+' >/dev/null')
	print('Assembly completed. Output can be found in the folder: '+prefix+'_spades')


def filter_rec_per_query(blast_results):
	pergene_rec={}
	filtered_hits=[]
	#sort blast hits by gene
	for l in blast_results:
		try:pergene_rec[l.split()[0]].append(l)
		except KeyError:pergene_rec[l.split()[0]]=[l]
	#get one best rec per gene
	for g in pergene_rec.keys():
		blast_lst=[l.strip().split() for l in pergene_rec[g]]
		blast_table= pd.DataFrame(blast_lst, columns =['V'+str(i) for i in list(range(12))]) 
		#sort dataframe
		#blast_table["V2"] = pd.to_numeric(blast_table["V2"])
		#blast_table=blast_table.sort_values(by='V2', ascending=False)
		blast_table['ref_start'] = [0] *len(blast_table)
		blast_table['ref_end'] = [0] *len(blast_table)
		for i in range(0,len(blast_table)):
			blast_table.iloc[i,12]=min(int(blast_table.iloc[i,8]),int(blast_table.iloc[i,9]))
			blast_table.iloc[i,13]=max(int(blast_table.iloc[i,8]),int(blast_table.iloc[i,9]))
		#blast_table=blast_table.sort_values(by='ref_start')
		blast_table["V6"] = pd.to_numeric(blast_table["V6"])
		blast_table["V7"] = pd.to_numeric(blast_table["V7"])
		blast_table=blast_table.sort_values(by='V6')
		cur_rec=blast_table.iloc[0]
		#consolidate ranges and preserve multiple hits on one sequences
		for i in range(0,len(blast_table)):
			if blast_table.iloc[i,6]<=cur_rec.V7:
				#extend range
				if blast_table.iloc[i,7]>cur_rec.V7:
					cur_rec.V7=blast_table.iloc[i,7]
					if blast_table.iloc[i,13]>cur_rec.ref_end:
						cur_rec.ref_end=blast_table.iloc[i,13]
			else:
				#different region
				cur_rec.V6=str(cur_rec.V6); cur_rec.V7=str(cur_rec.V7); cur_rec.ref_start=str(cur_rec.ref_start); cur_rec.ref_end=str(cur_rec.ref_end)
				filtered_hits.append('\t'.join(cur_rec))
				cur_rec=blast_table.iloc[i]
		#append the last record
		cur_rec.V6=str(cur_rec.V6); cur_rec.V7=str(cur_rec.V7); cur_rec.ref_start=str(cur_rec.ref_start); cur_rec.ref_end=str(cur_rec.ref_end)
		filtered_hits.append('\t'.join(cur_rec))
	#return list of strings as best hits
	return(filtered_hits)



def gene_assem(blast_results,assemb,species_name,outdir):
	#initiate values
	cur_gene=blast_results[0].split('\t')[1].split('_')[0]
	#if '_old' in blast_results[0]:cur_gene=cur_gene+'_old'
	blast_syntax={}
	blast_syntax[cur_gene]=[blast_results[0]]
	for l in blast_results:
		gene=l.split('\t')[1].split('_')[0]
		#if '_old' in l:
		#	gene=gene+'_old'
		if gene==cur_gene:
			blast_syntax[gene].append(l)
		else:
			#the blast syntax for last gene is complete, process them to extract sequences
			hits=filter_rec_per_query(blast_syntax[cur_gene])
			hits=[l.split() for l in hits]
			hit_table=pd.DataFrame(hits, columns =['V'+str(i) for i in list(range(14))]) 
			#hit_table["V2"] = pd.to_numeric(hit_table["V2"])
			hit_table["V10"] = pd.to_numeric(hit_table["V10"])
			hit_table['V12'] = pd.to_numeric(hit_table["V12"])
			hit_table['V13'] = pd.to_numeric(hit_table["V13"])
			hit_table=hit_table.sort_values(by='V12')
			if len(hit_table)==1:
			#if only one rec, output it
				if int(hit_table.iloc[0,8])>int(hit_table.iloc[0,9]):
					seq_str = str(assemb[hit_table.iloc[0,0]].seq[(int(hit_table.iloc[0,6])-1):int(hit_table.iloc[0,7])].reverse_complement())
				else:
					seq_str = str(assemb[hit_table.iloc[0,0]].seq[(int(hit_table.iloc[0,6])-1):int(hit_table.iloc[0,7])])
			#print(hit_table)
			#examine if the boundaries of exon hits are overlapping
			else:
			#multiple regions
				cur_region = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
				cur_region = cur_region.append(hit_table.iloc[0])
				right_end = cur_region.iloc[0,13]
				output_recs = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
				for i in range(0,len(hit_table)):
					if hit_table.iloc[i,12]-right_end>-10:
						#move to the next block, get optimum hit for the current region
						cur_region = cur_region.sort_values(by='V10')
						#print(cur_region)
						output_recs = output_recs.append(cur_region.iloc[0])
						cur_region = pd.DataFrame(columns=['V'+str(i) for i in list(range(14))])
						cur_region = cur_region.append(hit_table.iloc[i])
						right_end = cur_region.iloc[0,13]
					else:
						#left end overlap
						cur_region = cur_region.append(hit_table.iloc[i])
						if hit_table.iloc[i,13]>right_end:right_end=hit_table.iloc[i,13]
				#add the last rec
				cur_region = cur_region.sort_values(by='V10')
				output_recs = output_recs.append(cur_region.iloc[0])
				seq_str=''
				for i in range(0,len(output_recs)):
					#reverse compliment the sequence if needed
					if int(output_recs.iloc[i,8])>int(output_recs.iloc[i,9]):
						seq_str=seq_str+'NNNNN'+str(assemb[output_recs.iloc[i,0]].seq[(int(output_recs.iloc[i,6])-1):int(output_recs.iloc[i,7])].reverse_complement())
					else:
						seq_str=seq_str+'NNNNN'+str(assemb[output_recs.iloc[i,0]].seq[(int(output_recs.iloc[i,6])-1):int(output_recs.iloc[i,7])])
			#write sequence to output file
			out=open(outdir+'/'+cur_gene+'.fas','a')
			out.write('>'+species_name+'\n'+seq_str+'\n')
			out.close()
			#update the name of the working gene
			cur_gene=gene
			blast_syntax[gene]=[l]


mode=args.m
if mode =='submission':
	try:
		print('Generating submission commands for '+str(len(open(args.s).readlines())-1)+' species...')
		submiter_gen(args.b,args.s,args.o)
		print('Done.\nTo submit the jobs to the cluster, type:\nsh <output file name>')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m submision -b <bash file> -s <sample sheet> -o <output file name>')
	except IOError as e:print(e.errno)
elif mode =='qc':
	try:
		if args.suffix:
			sp_list=[i for i in os.listdir(args.i) if i.endswith(args.suffix)]
			print('processing '+str(len(sp_list))+' species for QC analysis...')
			qc_assemblies_only(sp_list,args.i,args.o)
		else:
			if args.s:
				sp_list=open(args.s).readlines()
				sp_list=[i.split()[0] for i in sp_list[1:]]
			else:
				sp_list=[i for i in os.listdir(args.i) if os.path.isdir(args.i+'/'+i)]
			if len(sp_list)==0:
				print('############################################################\n\
#ERROR:No GetOrganelle output folders found. If you are using assembly fasta files only, be sure to add the -suffix argument.\n\
Usage:\n\
With assemblies only:\n\
python phyloherb.py -m qc -i <input dir> -o <output dir> -suffix <suffix>\n\
With Getorganelle output folders:\n\
python phyloherb.py -m qc -i <input dir> -o <output dir> [optional] -s <sample sheet>')
				exit()
			print('processing '+str(len(sp_list))+' species for QC analysis...')
			qc(sp_list,args.i,args.o)
		print('output assembly_sum.tsv and assembly fasta files to '+args.o)
		print('Done.')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
With Getorganelle output folders:\n\
python phyloherb.py -m qc -i <input dir> -o <output dir> [optional] -s <sample sheet>\n\
With assemblies only:\n\
python phyloherb.py -m qc -i <input dir> -o <output dir> -suffix <suffix>')
elif mode =='ortho' and not args.nuc:
	try:
		PH_path=os.path.dirname(os.path.abspath(__file__))
		#print(PH_path)
		#get species list
		if args.sp:
			species=open(args.sp).readlines()
			species=[j.strip() for j in species]
			print('Using custom species set in '+args.sp+': '+', '.join(species))
		else:
			species=os.listdir(args.i)
			if args.suffix:
				species=[j for j in species if j.endswith(args.suffix)]
			else:
				species=[j for j in species if j.endswith('.assembly.fas')]
			print('Using all '+str(len(species))+' species found in '+args.i+': '+', '.join([j.split('.')[0] for j in species]))
		if len(species)==0:
			print('############################################################\n\
#ERROR:Zero species found! It looks like your assemblies does not end with \'.assembly.fas\'. Please use the -suffix flag!\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue> -g <gene list> -sp <species list> -l <minimum length for blast> -ref <custom ref seq> -mito <mito mode> -rdna <rDNA mode> -nuc <nuclear mode>')
			exit()
		#get minimum length for blast hit
		if args.l:
			min_len=int(args.l)
		else:min_len=60
		if args.n:
			threads=args.n
		else:threads='1'
		if args.evalue:
			evalue=args.evalue
		else:
			evalue='1e-20'
		print('Using length cutoff ' + str(min_len)+' bp for BLAST result filtering')
		#genes=["ycf2","ycf1","rpoC2","rpoB","rpoC1","rrn23","ndhF","ndhB","psaB","ndhA","clpP","ycf3","psbB","atpA","matK","rpl2","ndhD","atpB","rrn16","accD","rbcL","psbC","atpF","psaA","rps16","ndhH","psbA","psbD","rpoA","trnE-UUC","ccsA","petA","trnS-CGA","atpI","ndhK","rps2","cemA","rps3","petB","rps4","ycf4","ndhG","petD","ndhI","ndhJ","rps7","rps11","rpl22","rps8","atpE","rpl14","ndhC","rpl16","rpl20","rps18","ndhE","rps14","rps19","rpl23","rps15","psbE","atpH","psaC","psbH","rpl33","ycf15","psbZ","psbK","psaJ","pbf1","psbJ","rrn5","psbF","psbL","rpl32","psaI","petG","rpl36","psbI","psbT","psbM","petL","petN"]
		if args.rdna:
			#rDNA mode
			min_len=10
			print('Using build-in ribosomal gene set')
			genes=['18S','5.8S','28S']
			reference=PH_path+'/database/rDNA_reference.fas'
			for sp in species:
				ortho_extraction(sp,reference,args.i,args.o,genes,min_len,threads,evalue)
				#get ITSs
				get_ITS(sp,sp+'.blast.out',args.i,args.o,min_len)
		else:
			#get gene list
			if args.g:
				genes=open(args.g).readlines()
				genes=[j.strip() for j in genes]
				print('Using custom gene set')
			elif args.ref:
				genes=open(args.ref).readlines()
				genes=[j[1:].split('_')[0].strip() for j in genes if j.startswith('>')]
				genes=list(set(genes))
				print('Using custom gene set from the reference sequence: '+args.ref)
			elif args.mito:
				genes=open(PH_path+'/database/mito_gene.list').readlines()
				genes=[j.strip() for j in genes]
				print('Using build-in mitochondrial gene set')
			else:
				genes=open(PH_path+'/database/plastid_gene.list').readlines()
				genes=[j.strip() for j in genes]
				print('Using build-in plastid gene set')
			#get reference sequences
			if args.ref:
				reference=args.ref
			elif args.mito:
				reference=PH_path+'/database/mito_reference.fas'
			else:
				reference=PH_path+'/database/plastid_reference.fas'
			#extract blast hits:				
			for sp in species:
				ortho_extraction(sp,reference,args.i,args.o,genes,min_len,threads,evalue)	
		print('Completed gene extraction for '+str(len(species))+' species. Orthologous gene sequences are written to the output directory: '+args.o+'/.')		
	except TypeError:
			print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue> -g <gene list> -sp <species list> -l <min length for blast> -ref <custom ref seq> -mito <mito mode> -rdna <rDNA mode> -nuc <nuclear mode>')
	except FileNotFoundError:
		print('############################################################\n\
#ERROR:Input files not found or not in FASTA format!\n')
	except IOError as e:print(e.errno)
elif mode =='assemb':
	#check dependencies
	exit_now=0
	if not shutil.which('bowtie2'):
		print('############################################################\n\
#ERROR: Bowtie2 not in the environment. Please add its path or reinstall!')
		exit_now=1
	if not shutil.which('spades.py'):
		print('############################################################\n\
#ERROR: Spades not in the environment. Please add its path or reinstall!')
		exit_now=1
	if not shutil.which('samtools'):
		print('############################################################\n\
#ERROR: samtools not in the environment. Please add its path or reinstall!')
		exit_now=1
	if exit_now==1:quit()
	try:
		threads='1'
		if args.n:threads=args.n
		prefix=args.prefix
		reference=args.ref
		rs='NA'
		r1='NA'
		r2='NA'
		if args.rs:rs=args.rs
		else:
			r1=args.r1
			r2=args.r2
			if r1==r2 and r1!='NA':
				print('############################################################\n\
#ERROR: The forward and reverse reads must be different!\n\
Usage:\n\
python phyloherb.py -m assemb -ref <reference fasta> -prefix <prefix> [optional] -r1 <R1.fq> -r2 <R2.fq> -rs <Single1.fq,Single2.fq> -n <threads>')
				quit()
		if not (rs=='NA' and r1=='NA' and r2=='NA'):
			assembly(r1,r2,rs,threads,prefix,reference)
		else:
			print('############################################################\n\
#ERROR:Supply at least the pair-end reads or single-end reads\n\
Usage:\n\
python phyloherb.py -m assemb -ref <reference fasta> -prefix <prefix> [optional] -r1 <R1.fq> -r2 <R2.fq> -rs <Single1.fq,Single2.fq> -n <threads>')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m assemb -ref <reference fasta> -prefix <prefix> [optional] -r1 <R1.fq> -r2 <R2.fq> -rs <Single1.fq,Single2.fq> -n <threads>')
elif mode =='ortho' and args.nuc:
	try:
		import pandas as pd
	except ModuleNotFoundError:
		print('Module pandas not found. Please install.')
		quit()
	try:
		reference_ID=args.ref
		reference_ID=reference_ID.split('/')[-1]
		reference_ID=reference_ID.split('.')[0]
		threads='1'
		if args.n:threads=args.n
		evalue='1e-40'
		if args.evalue:evalue=args.evalue
		#if os.path.exists(reference_ID+'.nhr'):
		#	print('Use exsisting BLAST database.\nBLAST evalue threashold: '+evalue)				
		while True:
			print('Building BLAST database.\nBLAST evalue threashold: '+evalue)
			blastCommand='makeblastdb -in '+args.ref+' -dbtype nucl -out '+reference_ID+' >/dev/null'
			os.system(blastCommand)
		if not os.path.isdir(args.o):os.mkdir(args.o)
		#get assemblies
		if args.suffix:
			assemb_list=os.listdir(args.i)
			assemb_list=[j for j in assemb_list if j.endswith(args.suffix)]
			if len(assemb_list)>0:
				print('Using '+str(len(assemb_list)) +' species with suffix '+args.suffix+' in the directory: '+args.i)
				for sp in assemb_list:
					spID=sp.split(args.suffix)[0]
					print('Processing '+spID)
					#blast
					blastCommand='blastn -task dc-megablast -db '+reference_ID+' -num_threads '+threads+' -query '+args.i+'/'+sp+' \
-outfmt 6 -evalue '+evalue+' | sort -k2,2 -k12,12gr -k11,11g -k3,3gr - >'+spID+'.blast'
					os.system(blastCommand)
					#gene assembly
					blast_results=open(spID+'.blast').readlines()
					#add a dummy ending line
					blast_results=blast_results+['end\tend\tdummy\tline\t37\t2\t60\t383\t1315\t1000\t4.70e-106\t379\n']
					assemb=SeqIO.index(args.i+'/'+sp,'fasta')
					gene_assem(blast_results,assemb,spID,args.o)
			else:
				print('############################################################\n\
#ERROR:Zero species found! It looks like your assemblies does not end with '+args.suffix+'. Please Check!\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> -ref <custom ref seq> -nuc [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue>')
		else:
			print('No assembly suffix supplied. PhyloHerb is searching Spades outputs in the input directory: '+args.i)
			spades_list=[i for i in os.listdir(args.i) if os.path.isdir(args.i+'/'+i)]
			spades_list=[i.split('_spades')[0] for i in spades_list if os.path.exists(args.i+'/'+i+'/scaffolds.fasta')]
			if len(spades_list)>0:
				print('Using '+str(len(spades_list)) +' spades assemblies found in the input directory: '+args.i)
				for sp in spades_list:
					print('Processing '+sp)
					blastCommand='blastn -task dc-megablast -db '+reference_ID+' -num_threads '+threads+' -query '+args.i+'/'+sp+'_spades/scaffolds.fasta \
-outfmt 6 -evalue '+evalue+' | sort -k2,2 -k12,12gr -k11,11g -k3,3gr - >'+sp+'.blast'
					os.system(blastCommand)
					blast_results=open(sp+'.blast').readlines()
					blast_results=blast_results+['end\tend\t86.111\t324\t37\t2\t60\t383\t1315\t1000\t4.70e-106\t379\n']
					assemb=SeqIO.index(args.i+'/'+sp+'_spades/scaffolds.fasta','fasta')
					gene_assem(blast_results,assemb,sp,args.o)
			else:
				print('############################################################\n\
#ERROR:Zero species found! Please make sure this is the parent directory of SPADES outputs. If you have assemblies in FASTA, please use the -suffix flag.\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> -ref <custom ref seq> -nuc [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue>')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> -ref <custom ref seq> -nuc [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue>')
	except AttributeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m ortho -i <input dir> -o <output dir> -ref <custom ref seq> -nuc [optional] -suffix <assembly suffix> -n <number of threads> -evalue <evalue>')
elif mode =='conc':
	try:
		if args.g:
			files=open(args.g).readlines()
			files=[j.strip() for j in files]
			print('Concatenate '+str(len(files))+' alignments based on the custom gene set')
		else:
			files=os.listdir(args.i)
			files=[j for j in files if j.endswith(args.suffix)]
			print('Concatenate '+str(len(files))+' alignments in the directory '+args.i + ' with suffix '+args.suffix + ':'+', '.join(files))
		concatenation(args.i,files,args.o)
		print('Done.')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m conc -i <input dir> -o <output prefix> -suffix <alignment suffix> [optional] -g <gene list file>')
elif mode =='order':
	try:
		from ete3 import Tree
		if args.missing:
			missing=float(args.missing)
		else:
			missing=1
		files=os.listdir(args.i)
		files=[j for j in files if j.endswith(args.suffix)]
		print('Reorder ' +str(len(files))+' alignments based on the species tree '+args.t)
		order_aln(args.t,args.i,args.suffix,args.o,missing)
		print('Done.')
	except ModuleNotFoundError as e:
		print('############################################################\n\
#ERROR:Module not correctly loaded:')
		print(e)
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m order -t <species tree> -i <input dir> -o <output dir> -suffix <alignment suffix> [optional] -missing <missing proportion>')
elif mode =='getseq':
	try:
		print('Extracting gene and/or intergenic regions under the '+args.f+' mode from the Genbank files in '+args.i)
		if args.f=='gene':
			gene_extra(args.i,args.suffix,args.o)
			print('Done!')
		elif args.f=='genetic_block':
			geneblock_extra(args.i,args.suffix,args.o,args.gene_def)
			print('Done!')
		elif args.f=='intergenic':
			intergenic_extra(args.i,args.suffix,args.o,args.gene_def)
			print('Done!')
		else:
			print('############################################################\n\
#ERROR:Argument error!\n\
Usage:\n\
python phyloherb.py -m getseq -f <gene|genetic_block|intergenic> -i <input dir> -o <output dir> -suffix <genbank file suffix> [optional] -gene_def <gene definition file>')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python phyloherb.py -m getseq -f <gene|genetic_block|intergenic> -i <input directory> -o <output directory> -suffix <genbank file suffix> [optional] -gene_def <gene definition file>')
	except IOError as e:print(e)
else:
	print('############################################################\n\
#ERROR: Please choose one of the following execution mode using -m: qc, ortho, conc, order, getseq, assemb, submission\n\
')




