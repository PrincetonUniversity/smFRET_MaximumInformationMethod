

if [ -e $work_dir/job_id_list ]
then

	count=1
	
	while read id
	do
		if [ $count -eq 1 ]
		then
			job_ids=$id
		else
			job_ids=$job_ids:$id
		fi
		
		let count++
	done < $work_dir/job_id_list
	
	rm $work_dir/job_id_list
	
	depend_switch="-W depend=afterany:$job_ids"

else
	depend_switch=''
fi