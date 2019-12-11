#!/bin/bash

#duplicate_permutation.sh will take an R script for performing a series of permutations and create copies of the script. 
#usage without sbatching: duplicate_permutation.sh R_permutation_script.R 10
#usage with sbatching: duplicate_permutation.sh R_permutation_script.R 10 sbatch -t 04:00:00 --mem 24G

Rfile=$1; shift
nduplicates=$1; shift
run_sbatch=$1; shift

#pull the path for the permutation file output so that we we can change it iteratively
permfile_path=`grep 'output_permfile' ${Rfile} | grep -m 1 '='`
permfile_path_suffix=`echo ${permfile_path} | awk -F"'" '{print $2}'`

#the permfile may be encapsulated by double quotes -- so if its empty use double quotes instead to grab it
if [[ -z ${permfile_path_suffix} ]]; then permfile_path_suffix=`echo ${permfile_path} | awk -F'"' '{print $2}'`; fi

#pull the output directory path so that we can change it iteratively
output_path=`grep 'output_directory' ${Rfile} | grep -m 1 '='`
output_path_suffix=`echo ${output_path} | awk -F'"' '{print $2}'`
output_path_flag=2

#the output directory path may be encapsulated by single quotes -- so if its empty use single quotes instead to grab it
if [[ -z ${output_path_suffix} ]]; then output_path_suffix=`echo ${output_path} | awk -F"'" '{print $2}'`; output_path_flag=1 ;fi

for iter in `seq 1 $nduplicates`; do 
	sed "s|"${permfile_path_suffix}"|"${permfile_path_suffix}1"|"<${Rfile} >${Rfile}${iter}_temp.R
	if [[ $output_path_flag = 2 ]]; then
		sed "s|"${output_path_suffix}"|"${output_path_suffix}/perm${iter}/"|"<${Rfile}${iter}_temp.R >${Rfile}${iter}.R
	fi
	if [[ $output_path_flag = 1 ]]; then
		sed "s|${output_path_suffix}'|${output_path_suffix}/perm${iter}'|"<${Rfile}${iter}_temp.R >${Rfile}${iter}.R
	fi
	rm ${Rfile}${iter}_temp.R
	mkdir ${output_path_suffix}/perm${iter}
	if [[ ! -z ${run_sbatch} ]]; then
		if [ $run_sbatch = sbatch ]; then
			sbatch "$@" run_permutation.sh ${Rfile}${iter}.R
			sleep 10
		fi
	fi
done
