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
	cd $work_dir/$A
	if [ ! -e hist$A.out ]
	then
		jobname=hist$A
		
		echo "
			boots *.fr -a.$A -m -b25 -ohist$A.out $XT_SWITCH
		" | qsub -d$work_dir/$A -kn -M nomail -m n -N $jobname -o boots_out -e boots_err >> $work_dir/job_id_list
	
		let boots_count++
	fi
done
echo
echo "$boots_count boots jobs have been submitted"



cd $final_dir
rm *.out -f > /dev/null
for A in $alphas
do
	ln -s $work_dir/$A/hist$A.out ${sample}$A.out
done


frag_cnt=1

num_jobs=20 
submit_job_cnt=1
submit_string=''

for A in $alphas
do
	cd $work_dir/$A
	total_cnt=$(ls *.fr|wc -l)
	traj=1
	for F in *.fr
	do
		if [ ! -e ${F%%.fr}.xc ]
		then
			
			jobname=fr${A}_$traj
			

			if [[ $total_cnt -le $traj ]]
			then
				id_list=$work_dir/job_id_list_master
				source $scripts_dir/get_job_id_list.sh
				frag_cnt=1
				submit_job_cnt=$num_jobs
			else
				id_list=$work_dir/job_id_list
				depend_switch=''
			fi
			
			submit_string="$submit_string $F"
			let submit_job_cnt++
			
			if [ $submit_job_cnt -ge $num_jobs ]
			then
				echo "
					for FRAG in $submit_string
					do
						frag \$FRAG -a.$A $XT_SWITCH
					done	
				" | qsub -d $work_dir/$A -kn -M nomail -m n -N $jobname $depend_switch -o /dev/null -e /dev/null >> $id_list
			
				submit_job_cnt=1
				submit_string=''	
				
			fi
			
			let frag_count++	
		fi
		let traj++
	done	
done


echo "$frag_count frag jobs have been submitted"

exit 0
