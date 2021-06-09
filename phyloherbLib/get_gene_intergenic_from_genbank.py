import fnmatch
import os
from Bio import SeqIO

filenames=list()
for file in os.listdir('.'):
	 if fnmatch.fnmatch(file, '*.gb'):filenames.append(file)

#gene=["ycf2","ycf1","rpoC2","rpoB","rpoC1","rrn23","trnK-UUU","ndhF","ndhB","psaB","ndhA","clpP","ycf3","psbB","atpA","matK","rpl2","ndhD","atpB","rrn16","accD","rbcL","psbC","atpF","psaA","rps16","ndhH","psbA","psbD","rpoA","trnE-UUC","ccsA","petA","trnS-CGA","atpI","ndhK","rps2","cemA","trnV-UAC","rps3","petB","trnL-UAA","rps4","ycf4","ndhG","petD","ndhI","ndhJ","rps7","rps11","rpl22","rps8","atpE","rpl14","ndhC","rpl16","rpl20","rps18","ndhE","rps14","rps19","rpl23","rps15","psbE","atpH","psaC","psbH","rpl33","ycf15","psbZ","psbK","psaJ","pbf1","psbJ","rrn5","psbF","psbL","rpl32","psaI","petG","rpl36","psbI","psbT","psbM","trnA-UGC","rrn4.5","petL","petN"]

########################################
#get several additional gene sequences
gene2write=['trnE-UUC', 'trnS-CGA', 'pbf1', 'psbT', 'psbM', 'rrn4.5', 'petL', 'petN']


for f in filenames:
	gb_recs=SeqIO.read(f,'genbank')
	for feature in gb_recs.features:
		if feature.type=='gene':
			if feature.qualifiers['gene'][0] in gene2write:
				output_handle=open(feature.qualifiers['gene'][0]+'.fas','a')
				d=output_handle.write(">%s\n%s\n" % (feature.qualifiers['gene'][0],feature.extract(gb_recs).seq))
				output_handle.close()

########################################
#get some combined gene blocks
#region_ID	start	end
#INT22	ndhJ	ndhC
#INT32	psbJ	psbE
#INT34	petL	psaJ
#INT36	rpl33	rpl20
#INT40	psbT	psbH
#INT44	rps11	rpl36
#INT47	rpl14	rpl16
#INT49	rps3	rpl22
#INT61	trnA-UGC	trnE-UUC
#INT67	rpl23	rps19
#INT68	petN	psbM	coding
#INT69	psaI	ycf4	coding
#INT70	rpl32	ccsA	coding
#INT71	atpA	atpH	coding
#INT72	ndhD	ndhG	coding
#INT73	ndhH	rps15	coding

i=0
for f in filenames:
	gb_recs=SeqIO.read(f,'genbank')
	i=i+1
	combined_block={}
	for feature in gb_recs.features:
		if feature.type=='gene':
			if feature.qualifiers['gene'][0]=='ndhJ' or feature.qualifiers['gene'][0]=='ndhC':
				try:
					combined_block['INT22']=combined_block['INT22']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT22']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='psbJ' or feature.qualifiers['gene'][0]=='psbE':
				try:
					combined_block['INT32']=combined_block['INT32']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT32']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='petL' or feature.qualifiers['gene'][0]=='psaJ':
				try:
					combined_block['INT34']=combined_block['INT34']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT34']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rpl33' or feature.qualifiers['gene'][0]=='rpl20':
				try:
					combined_block['INT36']=combined_block['INT36']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT36']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='psbT' or feature.qualifiers['gene'][0]=='psbH':
				try:
					combined_block['INT40']=combined_block['INT40']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT40']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rps11' or feature.qualifiers['gene'][0]=='rps36':
				try:
					combined_block['INT44']=combined_block['INT44']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT44']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rpl14' or feature.qualifiers['gene'][0]=='rpl16':
				try:
					combined_block['INT47']=combined_block['INT47']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT47']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rps3' or feature.qualifiers['gene'][0]=='rpl22':
				try:
					combined_block['INT49']=combined_block['INT49']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT49']=[int(feature.location.start),int(feature.location.end)]
			#the following two are located in IR regions, need to be filtered
			if feature.qualifiers['gene'][0]=='trnA-UGC' or feature.qualifiers['gene'][0]=='trnE-UUC':
				if int(feature.location.start)>105000 and int(feature.location.start)<108000:
					try:
						combined_block['INT61']=combined_block['INT61']+[int(feature.location.start),int(feature.location.end)]
					except KeyError:
						combined_block['INT61']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rpl23' or feature.qualifiers['gene'][0]=='rps19':
				if int(feature.location.start)>86000 and int(feature.location.start)<92000:
					try:
						combined_block['INT67']=combined_block['INT67']+[int(feature.location.start),int(feature.location.end)]
					except KeyError:
						combined_block['INT67']=[int(feature.location.start),int(feature.location.end)]
			
			#few more after first round of blast and phylogeny
			if feature.qualifiers['gene'][0]=='petN' or feature.qualifiers['gene'][0]=='psbM':
				try:
					combined_block['INT68']=combined_block['INT68']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT68']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='psaI' or feature.qualifiers['gene'][0]=='ycf4':
				try:
					combined_block['INT69']=combined_block['INT69']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT69']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='rpl32' or feature.qualifiers['gene'][0]=='ccsA':
				try:
					combined_block['INT70']=combined_block['INT70']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT70']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='atpA' or feature.qualifiers['gene'][0]=='atpH':
				try:
					combined_block['INT71']=combined_block['INT71']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT71']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='ndhD' or feature.qualifiers['gene'][0]=='ndhG':
				try:
					combined_block['INT72']=combined_block['INT72']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT72']=[int(feature.location.start),int(feature.location.end)]
			if feature.qualifiers['gene'][0]=='ndhH' or feature.qualifiers['gene'][0]=='rps15':
				try:
					combined_block['INT73']=combined_block['INT73']+[int(feature.location.start),int(feature.location.end)]
				except KeyError:
					combined_block['INT73']=[int(feature.location.start),int(feature.location.end)]

	for k in combined_block.keys():
		start=min(combined_block[k])
		end=max(combined_block[k])
		seq=gb_recs.seq[(start-1):(end-1)]
		output_handle=open(k+'.fas','a')
		d=output_handle.write(">%s\n%s\n" % (k+'_'+str(i),seq))
		output_handle.close()


####################
#get intergenic region
x=open('/Users/limingcai/Downloads/intergenic_and_multigene_index.csv').readlines()
int_dict={}
for l in x[12:]:
	int_dict[l.split(',')[0]]=l.split(',')[1:3]
	
i=0
for f in filenames:
	gb_recs=SeqIO.read(f,'genbank')
	i=i+1
	#get all gene positions. for genes in the IR region, this will be the position in the second IR.
	gene_start={}
	gene_end={}
	for feature in gb_recs.features:
		if feature.type=='gene':
			try:
				gene_name=feature.qualifiers['gene'][0]
				gene_start[gene_name]=int(feature.location.start)
				gene_end[gene_name]=int(feature.location.end)
			except KeyError:
				print(feature.qualifiers,f)
	
	#get intergenic regions and write to file
	for k in int_dict.keys():
		two_gene_pos=[]
		try:
			two_gene_pos=[gene_start[int_dict[k][0]],gene_end[int_dict[k][0]],gene_start[int_dict[k][1]],gene_end[int_dict[k][1]]]
			two_gene_pos.sort()
			seq=gb_recs.seq[(two_gene_pos[1]-1):(two_gene_pos[2]-1)]
			output_handle=open(k+'.fas','a')
			d=output_handle.write(">%s\n%s\n" % (int_dict[k][0]+'_'+int_dict[k][1]+'_'+str(i)+'_'+gb_recs.id,seq))
			output_handle.close()
		except KeyError:
			print(k,f)
	
	#write the first interval 
	seq=gb_recs.seq[0:(gene_start['psbA']-1)]
	output_handle=open('INT1.fas','a')
	d=output_handle.write(">%s\n%s\n" % ('start_psbA_'+str(i)+'_'+gb_recs.id,seq))
	output_handle.close()
	
