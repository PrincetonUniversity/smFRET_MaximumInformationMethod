#!/bin/bash


for A in $alphas
do
	if [ -d $work_dir/$A ]
	then
		cd $work_dir/$A
		if [ -e hist$A.out ]
		then
			hists="$A $hists"
		fi
	fi	
done


if [ -n "$hists" ]
then	
	echo The following analysis files already exist from a previous run:
	for H in $hists
	do
		echo $work_dir/$H/hist$H.out
	done
	echo 
	timeout=30
	echo "Would you like to recalculate these distributions (y/n)?"
	read -t $timeout keep_dist
	if [ -z $keep_dist ]
	then
		echo timed out ... keeping preexisting distributions
	else
		if [ "$keep_dist" = y ]
		then
			echo deleting old files
			for H in $hists
			do
				rm $work_dir/$H/hist$H.out
			done
			for A in $alphas
			do
				if [ -e $final_dir/$sample$A.meanX ]
				then
        				rm $final_dir/*.meanX
					echo deleting meanX files 
				fi
			done
		else
			echo Old files have been kept	
		fi
	fi	
fi



data_dirs_file=$final_dir/data_dirs.txt

err_msg="You must put a list of all directories \n
	containing raw data in the file: data_dits.txt \n
	Also, this script is only compatible with \n
	auto-trajectory selection so before alpha \n
	analysis can be performed run: \n\n
	
	(in bash)\t ./do_analysis frsort\n\n
	
	(in matlab)\t  select_traj"

if [ ! -e $data_dirs_file ]
then

	echo -e $err_msg		
	exit 1
fi	
	

if [ -n "$skip_data_selection" ]
then
	exit 0
fi	

for M in $work_dir $work_dir/m{1..4}
do
	if [ ! -d $M ]
	then
		mkdir $M
	fi	
done
	
use_archive=0	

while read DIR
do
	#get date
	date=$(echo $DIR | tr / ' '| awk '{print $3}')
		
	if [ $date == $USER ]
	then
		use_archive=1
		continue
		
	fi
	pass_file=$log_dir/${date}_pass_manual.txt

	if [ ! -e $pass_file ]
	then
		echo -e $err_msg
		echo 
		echo "File not found:"
		echo
		echo $pass_file
	fi
	
	num_pass_corr=$(cat ${DIR}good/pass_corr.txt 2>/dev/null| wc -l) 
	num_pass_man=$(cat $pass_file | wc -l)
	num_fail_man=$(cat $log_dir/${date}_fail_manual.txt | wc -l)
	
	let "num_total_man=$num_pass_man + $num_fail_man"
	
	if [ ! $num_total_man -eq $num_pass_corr ]
	then
		echo "sort_traj.m did not execute successfully in the dir:"
		echo 
		echo ${DIR}good #removed slash before good
		echo
		echo "Please re-run the program"
	fi
	

	if [ -h $work_dir/$date ]
	then
		rm $work_dir/$date
	fi	
	
	if [ ! -d $work_dir/$date ]
	then
		mkdir $work_dir/$date
	else
		rm $work_dir/$date/* 2>/dev/null
	fi		
	
	cd $work_dir/$date
	
	while read F
	do
		ln -s $F.* .
	done < $pass_file	
	
	
done < $data_dirs_file

for M in m{1..4}
do
	cd $work_dir/$M
	rm -r prune 2>/dev/null
	rm * 2>/dev/null
	
	if [ $use_archive -eq 1 ]
	then
		ls $data_archive/good/manual_select/*$M*.fr 1>/dev/null 2>/dev/null
		archived_files=${?}
	else
		archived_files=1
	fi	
	
	ls ../20*/*$M*fr 1>/dev/null 2>/dev/null
	data_files=${?}
	
	if [[ $data_files -eq 0 || $archived_files -eq 0 ]]
	then
		if [ $data_files -eq 0 ]
		then
			ln -s ../20*/*$M*.fr .
			ln -s ../20*/*$M*.log .
			ln -s ../20*/*$M*.xc .
		fi
		if [ $archived_files -eq 0 ]
		then
			ln -s $data_archive/good/manual_select/*$M*.fr .
			ln -s $data_archive/good/manual_select/*$M*.xc .
			ln -s $data_archive/good/manual_select/*$M*.log .
		fi		
		mkdir prune
	else
		cd $work_dir
		rm -r $M 
	fi		
done
	
	
cd $work_dir

cp $scripts_dir/matlab/sort_prune{'',.ctf} . 
cp $scripts_dir/matlab/run_sort_prune.sh .
./run_sort_prune.sh /usr/local/MATLAB/R2012a/ 2>/dev/null 1>/dev/null 
rm -rf sort_prune* 2>/dev/null
	
	
for M in m{1..4}
do
	if [ -d $work_dir/$M ]
	then
		if [ -e $work_dir/$M/outliers.txt ]
		then
			cd $work_dir/$M
			while read out
			do
				mv $out* prune/
			done < outliers.txt
		fi	
	fi
done
	
