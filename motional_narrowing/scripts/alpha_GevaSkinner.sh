#!/bin/bash

#copy files into directory
cp $scripts_dir/matlab/process_global_fit{'',.ctf} .
cp $scripts_dir/matlab/process_global_fit.m .
cp $scripts_dir/matlab/GevaSkinner.m .
cp $scripts_dir/matlab/fit_GevaSkinner_model_global.m .
cp $scripts_dir/matlab/alpha_plot_results.m .
cp $scripts_dir/matlab/run_process_global_fit.sh .

source $scripts_dir/get_job_id_list_master.sh

echo "
	echo 	\$HOSTNAME > host_info

        ./run_process_global_fit.sh /usr/local/MATLAB/R2012a/ > globalfit.out
	rm -rf process_global* fit_GevaSkinner_model_global.m
	" | qsub -d $final_dir -o global_fit_out -e global_fit_error -k n -N global_fit -l walltime=90000 -M nomail -m n $depend_switch >> /dev/null
	
	

echo "2-state fitting job submitted to scheduler"
echo 	
	
	
