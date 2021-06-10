#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 5-24:00
#SBATCH -p <partition>
#SBATCH --mem=30000
#SBATCH -o mafft_pasta.out
#SBATCH -e mafft_pasta.err

mkdir $1.pasta
cd $1.pasta

run_pasta.py -i ../$1\.pasta.aln.fas -a -t ../$1\.pasta_ref.tre -o .
mv pastajob.*$1\.*.missing_filtered.aln ../$1\.pasta_rnd2.aln.fas
mv pastajob.tre ../$1\.pasta.tree
cd ..
rm -r $1\.pasta
