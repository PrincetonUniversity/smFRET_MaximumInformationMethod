#!/bin/bash

cd $data_location
file_names_list=( $(echo $(ls *.fr|grep fr) ) ) 

cd $run_dir

if [ ! -d files ]
then
	mkdir files
fi

for A in $alphas
do
	if [ ! -d $A ]
	then 
		mkdir $A
	fi
done


cd files
file_test=$(ls -l|grep .fr)
if [ $? -ge 1 ]
then
	boot_count=1
else 
	boot_count=$(ls -l *.fr|wc -l)
fi




while [ $boot_count -le $num_files ]
do
	number=$RANDOM
	let "number %= $num_files"
	choosen_file=${file_names_list[${number}]}
	if [ ! -e $run_dir/files/$choosen_file ]
	then
		ln -s $data_location/$choosen_file $run_dir/files/
		fr_dest=$choosen_file
		xc_target=${choosen_file%%.fr}.xc
		xc_dest=${choosen_file%%.fr}.xc
	else
		inc=1
		while [ -e $run_dir/files/${choosen_file%%.fr}.$inc.fr ]
		do	
			let inc++
		done
		ln -s $data_location/$choosen_file $run_dir/files/${choosen_file%%.fr}.$inc.fr
		fr_dest=${choosen_file%%.fr}.$inc.fr
		xc_target=${choosen_file%%.fr}.xc
		xc_dest=${choosen_file%%.fr}.$inc.xc
	fi
	for A in $alphas
	do
		ln -s $work_dir/$A/$choosen_file $run_dir/$A/$fr_dest
		ln -s $work_dir/$A/$xc_target $run_dir/$A/$xc_dest
	done	
	
	
	
	let boot_count++
done