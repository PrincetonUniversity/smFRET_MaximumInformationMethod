#!/bin/bash


fr_dir=$(pwd)

cp $scripts_dir/matlab/sort_reg_I{'',.ctf} . 2>/dev/null


source $scripts_dir/get_job_id_list.sh

QSUB_OPTS="-d $fr_dir -o regI_output.out -e regI_run_error.out  -m n -M nomail -l walltime=1800 $depend_switch"

NAME=" -N reg_I"
echo "
	./sort_reg_I > r1_out
	rm -r sort_reg_I*
	" | qsub $QSUB_OPTS $NAME>> $work_dir/job_id_list
	
	
NAME=" -N batchhist"

source $scripts_dir/sort_batchhist.sh


date=$(echo $fr_dir | tr / ' '| awk '{print $3}')

if [ -e pass_manual.txt ]
then
	if [ ! -e $log_dir/${date}_pass_manual.txt ]
	then
		cp pass_manual.txt $log_dir/${date}_pass_manual.txt
	fi
fi

if [ -e fail_manual.txt ]
then
	if [ ! -e $log_dir/${date}_fail_manual.txt ]
	then
		cp fail_manual.txt $log_dir/${date}_fail_manual.txt
	fi
fi

	 
