from Bio import AlignIO
import os, argparse, sys, shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus

parser = argparse.ArgumentParser(description='PhyloHerb is a bioinfomatic utility wrappepr to process genome skimming data for phylogenomics studies.')
parser.add_argument('-m', metavar='mode', help='execution mode, choose one of the following [submission, qc, ortho, conc, order, getseq]', required=True)
parser.add_argument('-i', metavar='dir', help='input directory')
parser.add_argument('-o', metavar='dir', help='output directory')
parser.add_argument('-b', metavar='file', help='[submission mode] path to the bash file')
parser.add_argument('-s',  metavar='file', help='[submission mode] path to the taxon sampling sheet')
parser.add_argument('-sp',  metavar='file', help='[ortho mode] a file containing a list of species')
parser.add_argument('-g',  metavar='file', help='[ortho and conc mode] a file containing a list of loci')
parser.add_argument('-l',  metavar='integer', help='[ortho mode] minimum length of blast hits')
parser.add_argument('-ref',  metavar='file', help='[ortho mode] custom reference sequences')
parser.add_argument('-mito',  help='[ortho mode] extract mitochondrial genes using build-in references',action='store_true')
parser.add_argument('-rdna',  help='[ortho mode] extract nuclear ribosomal regions using build-in references',action='store_true')
parser.add_argument('-suffix', metavar='string', help='[conc mode] suffix of alignment files')
parser.add_argument('-t', metavar='file', help='[order mode] newick tree file to order alignments based on phylogeny')
parser.add_argument('-missing', metavar='float 0-1', help='[order mode] maximum proportion of missing data allowed for each species')
parser.add_argument('-f',  metavar='mode', help='[getseq mode] how to extract loci, choose one of the following [gene, genetic_block, intergenic]')
parser.add_argument('-gene_def',  metavar='file', help='[getseq mode] a gene delimitation file that defines genetic blocks')


args = parser.parse_args()

def submiter_gen(bash_file,sample_sheet,output):
	sp_sheet=open(sample_sheet).readlines()
	out=open(output,'a')
	for l in sp_sheet[1:]:
		d=out.write('sbatch '+bash_file+' '+l.split()[1]+' '+l.split()[2]+' '+l.split()[0]+'\n')
	out.close()


def qc(sample_sheet,input_dir,output_dir):
	failed=[]
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
				target_reads_file=os.listdir(input_dir+'/'+sp+'/seed/')
				target_reads_file=[j for j in target_reads_file if j.endswith('.fq')]
				target_reads=open(input_dir+'/'+sp+'/seed/'+target_reads_file[0]).readlines()
				target_reads=len(target_reads)/4
				out.write('\t'.join([sp,str(total_reads),str(target_reads),'NA','NA','NA','NA'])+'\n')
			except IOError:
				out.write('\t'.join([sp,'NA','NA','NA','NA','NA','NA'])+'\n')
				failed.append(sp)
	if len(failed)>0:
		print('Cannot find GetOrganelle outputs in the directory '+input_dir+' for the following species: '+', '.join(failed))
	
def ortho_extraction(sp,reference_seq,input_dir,output_dir,genes,min_len):
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	print('processing species '+sp)
	lib_ID=sp
	S= 'makeblastdb -in ' +input_dir+'/'+ lib_ID +'.assembly.fas -dbtype nucl -out '+lib_ID
	os.system(S)
	S = 'blastn -task dc-megablast -db '+lib_ID+' -query ' + reference_seq + ' -outfmt 6 -evalue 1e-20 -out '+ lib_ID +'.blast.out'
	os.system(S)
	x=open(lib_ID+'.blast.out').readlines()
	y=SeqIO.index(input_dir+'/'+lib_ID+'.assembly.fas','fasta')
	a={}
	for g in genes:
		#print g
		best=0
		a[g]=[r for r in x if g+'_' in r]
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
				SeqIO.write(SeqRecord(seq,lib_ID, '', ''),open(output_dir+'/'+g+'.fas','a'),'fasta')
		except (NameError,AttributeError):continue
	os.remove(lib_ID+'.nhr')
	os.remove(lib_ID+'.nin')
	os.remove(lib_ID+'.nsq')
	
def get_ITS(sp,blast_file,input_dir,output_dir,min_len):
	lib_ID=sp
	blast_txt=open(blast_file).readlines()
	y=SeqIO.index(input_dir+'/'+lib_ID+'.assembly.fas','fasta')
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
	if a['18S'].split('\t')[1] == a['5.8S'].split('\t')[1]:
		ITS1.append([int(a['18S'].split('\t')[8]),int(a['18S'].split('\t')[9]),int(a['5.8S'].split('\t')[8]),int(a['5.8S'].split('\t')[9])])
		ITS1.sort()
		start=ITS1[1]
		end=ITS1[2]
		seq=y.seq[(start-1):end]
		output_handle=open(output_dir+'/ITS1.fas','a')
		d=output_handle.write(">%s\n%s\n" % (loci+'_'+f.split('.')[0],seq))
		output_handle.close()
	if a['5.8S'].split('\t')[1] == a['28S'].split('\t')[1]:
		ITS2.append([int(a['5.8S'].split('\t')[8]),int(a['5.8S'].split('\t')[9]),int(a['28S'].split('\t')[8]),int(a['28S'].split('\t')[9])])
		ITS2.sort()
		start=ITS2[1]
		end=ITS2[2]
		seq=y.seq[(start-1):end]
		output_handle=open(output_dir+'/ITS2.fas','a')
		d=output_handle.write(">%s\n%s\n" % (loci+'_'+f.split('.')[0],seq))
		output_handle.close()



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
		y=SeqIO.parse(input_dir+'/'+g,'fasta')
		out=open(output_dir+'/'+'.'.join(g.split('.')[:-1])+'.ordered.fas','a')
		for rec in y:
			missing=float(rec.seq.count('-')+rec.seq.count('N'))/len(rec.seq)
			if rec.id in total_taxa and missing<max_missing:
				SeqIO.write(rec,out,'fasta')
				sp2preserve.append(rec.id)
		out.close()
		t.prune(list(set(total_taxa) & set(sp2preserve))) 
		t.write(format=1, outfile=output_dir+'/'+'.'.join(g.split('.')[:-1])+".pasta_ref.tre")

def concatenation(input_dir,files,output):
	nexus_filenames=[]
	if not os.path.isdir(output+'_tem'):os.mkdir(output+'_tem')
	for fn in files:
		#the alphabet argument will not be used, but not including it will trigger an error
		#x=AlignIO.read(input_dir+'/'+fn,'fasta',alphabet=Gapped(IUPAC.protein))
		new_filename=output+'_tem'+'/'+'.'.join(fn.split('.')[:-1])+'.nex'
		nexus_filenames.append(new_filename)
		#g = open(new_filename, "w")
		#d=g.write(x.format("nexus"))
		#g.close()
		AlignIO.convert(input_dir+'/'+fn, 'fasta', new_filename, 'nexus', molecule_type='DNA')
	nexi =  [(fname, Nexus.Nexus(fname)) for fname in nexus_filenames]
	combined = Nexus.combine(nexi)
	out=open(output+'.conc.nex', 'w')
	combined.write_nexus_data(out)
	out.close()
	d=AlignIO.convert(output+'.conc.nex', 'nexus', output+'.conc.fas', 'fasta', molecule_type='DNA')
	#os.rmdir(output+'_tem')
	shutil.rmtree(output+'_tem', ignore_errors=True)
	out=open(output+'.partition','a')
	x=open(output+'.conc.nex').readlines()
	begin_write=0
	for l in x:
		if l.startswith('charset'):
			begin_write=1
		elif l.startswith('charpartition'):
			begin_write=0
		if begin_write==1:out.write(l)
	out.close()

def gene_extra(input_dir,suffix,output_dir):
	filenames=os.listdir(input_dir)
	filenames=[j for j in filenames if j.endswith(suffix)]
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
				print('Cannot find the following genes in the Genbank annotation: '+l)

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
				print('Cannot find the following genes in the Genbank annotation: '+l)


mode=args.m
print('############################################################\n\
PhyloHerb v1.0\n\
A bioinformatic pipeline for herbariomics based biodiversity research\n')
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
		print('processing '+str(len(open(args.s).readlines())-1)+' species for QC analysis...')
		qc(args.s,args.i,args.o)
		print('Done.')
	except TypeError:
		print('############################################################\n\
		#ERROR:Insufficient arguments!\n\
		Usage:\n\
		python phyloherb.py -m qc -s <sample sheet> -i <input directory containing Getorganelle output> -o <output directory>')
elif mode =='ortho':
	try:
		PH_path=os.path.dirname(__file__)
		#print(PH_path)
		#get species list
		if args.sp is not None:
			species=open(args.sp).readlines()
			species=[j.strip() for j in species]
			print('Using custom species set in '+args.sp+': '+', '.join(species))
		else:
			species=os.listdir(args.i)
			species=[j.split('.')[0] for j in species if j.endswith('.assembly.fas')]
			print('Using all species found in '+args.i+': '+', '.join(species))
		if len(species)==0:
			print('############################################################\n\
			#ERROR:Zero species found! Check your input!\n\
			Usage:\n\
			python phyloherb.py -m ortho -i <input directory> -o <output directory> [optional] -g <gene list file> -sp <species list> -l <minimum length for blast hit> -ref <fasta file of custom reference> -mito <use build-in mitochondrial reference genes>')
		#get minimum length for blast hit
		if args.l is not None:
			min_len=int(args.l)
		else:min_len=60
		print('Using length cutoff ' + str(min_len)+' bp for BLAST result filtering')
		#genes=["ycf2","ycf1","rpoC2","rpoB","rpoC1","rrn23","ndhF","ndhB","psaB","ndhA","clpP","ycf3","psbB","atpA","matK","rpl2","ndhD","atpB","rrn16","accD","rbcL","psbC","atpF","psaA","rps16","ndhH","psbA","psbD","rpoA","trnE-UUC","ccsA","petA","trnS-CGA","atpI","ndhK","rps2","cemA","rps3","petB","rps4","ycf4","ndhG","petD","ndhI","ndhJ","rps7","rps11","rpl22","rps8","atpE","rpl14","ndhC","rpl16","rpl20","rps18","ndhE","rps14","rps19","rpl23","rps15","psbE","atpH","psaC","psbH","rpl33","ycf15","psbZ","psbK","psaJ","pbf1","psbJ","rrn5","psbF","psbL","rpl32","psaI","petG","rpl36","psbI","psbT","psbM","petL","petN"]
		if args.rdna is not None:
			#rDNA mode
			print('Using build-in ribosomal gene set')
			genes=['18S','5.8S','28S']
			reference=PH_path+'/database/rDNA_reference.fas'
			for sp in species:
				ortho_extraction(sp,reference,args.i,args.o,genes,min_len)
			#get ITSs
			
		else:
			#get gene list
			if args.g is not None:
				genes=open(args.g).readlines()
				genes=[j.strip() for j in genes]
				print('Using custom gene set')
			elif args.mito is not None:
				genes=open(PH_path+'/database/mito_gene.list').readlines()
				genes=[j.strip() for j in genes]
				print('Using build-in mitochondrial gene set')
			else:
				genes=open(PH_path+'/database/plastid_gene.list').readlines()
				genes=[j.strip() for j in genes]
				print('Using build-in plastid gene set')
			#get reference sequences
			if args.ref is not None:
				reference=args.ref
			elif args.mito is not None:
				reference=PH_path+'/database/mito_reference.fas'
			else:
				reference=PH_path+'/database/plastid_reference.fas'
			#extract blast hits:				
			for sp in species:
				ortho_extraction(sp,reference,args.i,args.o,genes,min_len)	
		print('Completed gene extraction for '+str(len(species))+' species.')		
	except TypeError:
			print('############################################################\n\
		#ERROR:Insufficient arguments!\n\
		Usage:\n\
		python phyloherb.py -m ortho -i <input directory> -o <output directory> [optional] -g <gene list file> -sp <species list> -l <minimum length for blast hit> -ref <fasta file of custom reference> -mito <use build-in mitochondrial reference> -rdna <use build-in rDNA reference>')
	except IOError as e:print(e.errno)
elif mode =='conc':
	try:
		if args.g is not None:
			files=open(args.g).readlines()
			files=[j.strip() for j in files]
			print('Concatenate '+str(len(files))+' alignments based on the custom gene set')
		else:
			files=os.listdir(args.i)
			files=[j for j in files if j.endswith(args.suffix)]
			print('Concatenate '+str(len(files))+' alignments in the directory '+args.i + ' with suffix '+args.suffix)
		concatenation(args.i,files,args.o)
		print('Done.')
	except TypeError:
		print('############################################################\n\
		#ERROR:Insufficient arguments!\n\
		Usage:\n\
		python phyloherb.py -m conc -i <input directory> -o <output prefix> -suffix <alignment suffix> [optional] -g <gene list file>')
elif mode =='order':
	try:
		from ete3 import Tree
		if args.missing is not None:
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
		python phyloherb.py -m order -t <species tree> -i <input directory> -o <output directory> -suffix <alignment suffix> [optional] -missing <missing proportion>')
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
		python phyloherb.py -m getseq -f <gene|genetic_block|intergenic> -i <input directory> -o <output directory> -suffix <alignment suffix> [optional] -gene_def <gene definition file>')
	except TypeError:
		print('############################################################\n\
		#ERROR:Insufficient arguments!\n\
		Usage:\n\
		python phyloherb.py -m getseq -f <gene|genetic_block|intergenic> -i <input directory> -o <output directory> -suffix <alignment suffix> [optional] -gene_def <gene definition file>')
	except IOError as e:print(e)
else:
	print('############################################################\n\
	#ERROR: Please choose one of the following execution mode using -m: submission, qc, ortho, conc, order, getseq\n\
	')




