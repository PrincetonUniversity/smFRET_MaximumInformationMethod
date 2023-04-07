#!/bin/bash


if [ ! -d $boot_dir ]
then
	mkdir $boot_dir
fi


RANDOM=$PPID

for A in $alphas
do
	cd $work_dir/$A
	export num_files=$(ls -l *.fr|wc -l)
	export data_location=$(pwd)
	break
done

iter=1

while [ $iter -le $boot_num ]
do
	echo "working on Run $iter"
	export run_dir=$boot_dir/run$iter
	if [ ! -d $run_dir ]
	then
		mkdir $run_dir
	fi
	
	
	if [ ! -e $run_dir/${output_file}.mat ]
	then
		#first step - create resampled data sets
		$scripts_dir/bootstrap_choose_files.sh 
		
		#second step - submit jobs for boots to scheduler
		$scripts_dir/bootstrap_boots.sh
		
		#third step - calculate time resolutions
		$scripts_dir/bootstrap_gettr.sh
		
		#forth step - fit to 2-state model
		$scripts_dir/bootstrap_GevaSkiner.sh
	fi
	
	let iter++
done 



cp $scripts_dir/matlab/hist_error_gen.m $final_dir/
cp $scripts_dir/matlab/bootstrap_check_runs.m $final_dir/