#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 5-24:00
#SBATCH -p <partition>
#SBATCH --mem=30000
#SBATCH -o mafft_pasta.out
#SBATCH -e mafft_pasta.err

mafft --adjustdirection $1\.fas | sed 's/_R_//g' >$1\.mafft.aln.fas
mkdir $1\.pasta
cd $1\.pasta
run_pasta.py -i ../$1\.mafft.aln.fas
rm *.txt
rm *.gz
mv pastajob.tre ../$1\.pasta.tree
mv pastajob.*$1\.mafft.aln.fas.aln ../$1\.pasta.aln.fas
cd ..
rm -r $1\.pasta