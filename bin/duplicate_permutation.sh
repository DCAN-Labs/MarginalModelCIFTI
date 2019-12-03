#!/bin/bash

#duplicate_permutation.sh will take an R script for performing a series of permutations and create copies of the script. 
#usage: duplicate_permutation.sh R_permutation_script.R 10 /mnt/rose/shared/projects/ABCD/avg_pconn_maker/subset_analyses/analyses/scratch/perms_group1/permfile.txt sbatch

Rfile=$1; shift
nduplicates=$1; shift
permfile_path=$1; shift
run_sbatch=$1; shift

echo "$@"

for iter in `seq 1 $nduplicates`; do 
	sed "s|$permfile_path|${permfile_path}${iter}|"<${Rfile} >${Rfile}${iter}.R
	if [ $run_sbatch = sbatch ]; then
		sbatch "$@" run_permutation.sh ${Rfile}${iter}.R
		sleep 10
	fi
done
