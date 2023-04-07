#!/bin/bash


cd $run_dir


cp $final_dir/.output_file .

cp $scripts_dir/matlab/process_global_fit{'',.ctf} .
cp $scripts_dir/matlab/process_global_fit.m .
cp $scripts_dir/matlab/GevaSkinner.m .
cp $scripts_dir/matlab/fit_GevaSkinner_model_global.m .
cp $scripts_dir/matlab/alpha_plot_results.m .
cp $scripts_dir/matlab/run_process_global_fit.sh .

source $scripts_dir/get_job_id_list.sh


echo "
	cp $final_dir/.final_dir .
	echo 	\$HOSTNAME > host_info

        ./run_process_global_fit.sh /usr/local/MATLAB/R2012a/
	rm -rf process_global* fit_GevaSkinner_model_global.m
	" | qsub -d $run_dir -k n -N global_fit -l walltime=12000 -M nomail -m n $depend_switch >> /dev/null 
	
	
	
	
