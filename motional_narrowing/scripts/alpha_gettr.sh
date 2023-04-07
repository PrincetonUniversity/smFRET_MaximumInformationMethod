#!/bin/bash

cd $final_dir


cp $scripts_dir/matlab/alpha_gettr{'',.ctf} .
cp $scripts_dir/matlab/run_alpha_gettr.sh .


source $scripts_dir/get_job_id_list_master.sh


echo "
	./run_alpha_gettr.sh /usr/local/MATLAB/R2012a/
	rm -rf alpha_gettr*
" | qsub -d $final_dir -kn -M nomail -m n -N alpha_gettr $depend_switch -o alpha_gettr_out -e alpha_gettr_error >>  $work_dir/job_id_list_master

echo "alpha_gettr.m job submitted to scheduler"
echo
