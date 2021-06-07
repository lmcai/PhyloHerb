#!/bin/bash
#
#SBATCH -n 8                 # Number of cores
#SBATCH -N 1                 # Number of nodes for the cores
#SBATCH -t 0-10:00           # Runtime in D-HH:MM format
#SBATCH -p serial_requeue    # Partition to submit to
#SBATCH --mem=8000            # Memory pool for all CPUs
#SBATCH -o getorg.out      # File to which standard out will be written
#SBATCH -e getorg.err      # File to which standard err will be written


#Activate conda environment, assuming the name of the conda environment is 'getorg'
#module load Anaconda
#source activate getorg

#export DATA_DIR=[absolute path to input data]

#To assemble plant plastome
get_organelle_from_reads.py -1 $1 -2 $2 -o chl/$3 -R 15 -k 21,45,65,85,95,105 -F embplant_pt

#To assemble plant nuclear ribosomal RNA
get_organelle_from_reads.py -1 $1 -2 $2 -o ITS/$3 -R 10 -k 35,85,105,115 -F embplant_nr

#To assemble plant mitochondria:
get_organelle_from_reads.py -1 $1 -2 $2 -o mito/$3 -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt


#To submit the job to the cluster
#sbatch getorg.sh [forward reads] [backward reads] [output prefix]

###########################
#if your reads are single-ended, use the following commands instead
#To assemble plant plastome
#get_organelle_from_reads.py -1 $1 -o chl/$2 -R 15 -k 21,45,65,85,95,105 -F embplant_pt
#To assemble plant nuclear ribosomal RNA
#get_organelle_from_reads.py -1 $1 -o ITS/$2 -R 10 -k 35,85,105,115 -F embplant_nr
#To assemble plant mitochondria:
#get_organelle_from_reads.py -1 $1 -o mito/$2 -R 50 -k 21,45,65,85,105 -P 1000000 -F embplant_mt

#To submit the job to the cluster
#sbatch getorg.sh [forward reads] [backward reads] [output prefix]