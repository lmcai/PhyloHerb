echo $1

nohup makeblastdb -in $1.assembly.fas -dbtype nucl -out tem >/dev/null 2>&1
nohup blastn -task dc-megablast -db tem -query malp.chl.gene.fas -outfmt 6 -evalue 1e-40 -out $1.gene.blast.out >/dev/null 2>&1 

python extract_blast_hits.py $1

rm tem.*
rm $1.*.blast.out
