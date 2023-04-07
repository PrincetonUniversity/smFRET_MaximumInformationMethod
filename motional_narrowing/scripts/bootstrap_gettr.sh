#!/bin/bash


cd $run_dir

for CF in .alphas .sample .r0
do
	cp $final_dir/$CF .
done


echo $run_dir > .final_dir
echo $run_dir > .work_dir


cp $scripts_dir/matlab/alpha_gettr{'',.ctf} .
cp $scripts_dir/matlab/run_alpha_gettr.sh .

source $scripts_dir/get_job_id_list.sh


echo "
	./run_alpha_gettr.sh /usr/local/MATLAB/R2012a/
        rm -rf alpha_gettr*
" | qsub -d $run_dir -k n -M nomail -m n -N alpha_gettr $depend_switch >> $work_dir/job_id_list

