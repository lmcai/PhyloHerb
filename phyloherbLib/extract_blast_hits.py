import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

lib_ID=sys.argv[1]
x=open(lib_ID+'.gene.blast.out').readlines()
y=SeqIO.index('/n/davis_lab/Users/lmcai/MalpChl/assemblies/chl/'+lib_ID+'.fas','fasta')

gene=["ycf2","ycf1","rpoC2","rpoB","rpoC1","rrn23","trnK-UUU","ndhF","ndhB","psaB","ndhA","clpP","ycf3","psbB","atpA","matK","rpl2","ndhD","atpB","rrn16","accD","rbcL","psbC","atpF","psaA","rps16","ndhH","psbA","psbD","rpoA","trnE-UUC","ccsA","petA","trnS-CGA","atpI","ndhK","rps2","cemA","trnV-UAC","rps3","petB","trnL-UAA","rps4","ycf4","ndhG","petD","ndhI","ndhJ","rps7","rps11","rpl22","rps8","atpE","rpl14","ndhC","rpl16","rpl20","rps18","ndhE","rps14","rps19","rpl23","rps15","psbE","atpH","psaC","psbH","rpl33","ycf15","psbZ","psbK","psaJ","pbf1","psbJ","rrn5","psbF","psbL","rpl32","psaI","petG","rpl36","psbI","psbT","psbM","trnA-UGC","petL","petN","INT22","INT32","INT34","INT36","INT40","INT44","INT47","INT49","INT61","INT67",'INT68','INT69','INT70','INT71','INT72']

a={}
for g in gene:
	#print g
	best=0
	a[g]=(r for r in x if g+'_' in r)
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
		if end-start>60:
			seq=y[hit].seq[(start-1):(end-1)]
			SeqIO.write(SeqRecord(seq,lib_ID, '', ''),open('/n/davis_lab/Users/lmcai/MalpChl/alignment_phylogeny/chl/gene_aln/'+g+'.fas','a'),'fasta')
	except (NameError,AttributeError):continue

xx=open(lib_ID+'.intergenic.blast.out').readlines()
intergene=["INT1","INT2","INT3","INT4","INT5","INT6","INT7","INT8","INT9","INT10","INT11","INT12","INT13","INT14","INT15","INT16","INT17","INT18","INT19","INT20","INT21","INT23","INT24","INT25","INT26","INT27","INT28","INT29","INT30","INT31","INT33","INT35","INT37","INT38","INT39","INT41","INT42","INT43","INT45","INT46","INT48","INT50","INT51","INT52","INT53","INT54","INT55","INT56","INT57","INT58","INT59","INT60","INT62","INT63","INT64","INT65"]

a={}
for g in intergene:
	#print g
	best=0
	a[g]=(r for r in xx if g+'_' in r)
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
		if end-start>60:
			seq=y[hit].seq[(start-1):(end-1)]
			SeqIO.write(SeqRecord(seq,lib_ID, '', ''),open('/n/davis_lab/Users/lmcai/MalpChl/alignment_phylogeny/chl/intergenic_aln/'+g+'.fas','a'),'fasta')
	except (NameError,AttributeError):continue

