#!/bin/bash


if [ -e $final_dir/.crosstalk ]
then
	XT_SWITCH="-x$(cat $final_dir/.crosstalk)"

else
	XT_SWITCH=''
fi


boots_count=0
for A in $alphas
do
	cd $run_dir/$A
	if [ ! -e hist$A.out ]
	then
		jobname=hist$A
		
		echo "
			boots *.fr -a.$A -m -b5 -ohist$A.out $XT_SWITCH
 		" | qsub -d$run_dir/$A -kn -M nomail -m n -N $jobname >> $work_dir/job_id_list
	
		let boots_count++
	fi
done


cd $run_dir
rm *.out -f > /dev/null
for A in $alphas
do
	ln -s $run_dir/$A/hist$A.out ${sample}$A.out
done