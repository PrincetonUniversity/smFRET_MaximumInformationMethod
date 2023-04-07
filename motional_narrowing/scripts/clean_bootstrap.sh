#!/bin/bash


boots_count=1
while [ $boots_count -le $boot_num ]
do
	export run_dir=$boot_dir/run$boots_count
	echo "working on $run_dir"
	cd $run_dir
	if [ -e $run_dir/${output_file}.mat ]
	then
		rm *.mat
	fi

	for A in $alphas
	do
		cd $run_dir/$A
		if [  -e hist$A.out ]
		then
			rm hist$A.out
			echo "cleaning on $run_dir/$A"
		fi
	done
	let boots_count++
done


