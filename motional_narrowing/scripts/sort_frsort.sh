# This is the main script which controls automatic selection of 
# fret (.fr) files
# 
# The tasks are is as follows
# 
# 1 - analyze each trajectory with frag
# 2 - sort by signal to noise
# 3 - analyze each trajectory with a correlation function
# 4 - analyze region I of the trajectory

	

if [ ! -e data_dirs.txt ]
then
	echo
	echo "Can't find the file: data_dirs.txt"
	echo "This file should contain the path to"
	echo "each directory containing the primary"
	echo "fret data to be analyzed, exiting ..."
	
	exit 1
fi


if [ -n "$skip_data_selection" ]
then
	echo "This script does not need to be run"
	echo "trajectory selection has been skipped in config.alpha"
	echo
	
	exit 1
fi	
	
while read fr_dir
do	
	
	if [ -e $work_dir/job_id_list ]
	then
		rm $work_dir/job_id_list
	fi
	
	cp $scripts_dir/matlab/select_traj.m . 2>/dev/null
	

	date=$(echo $fr_dir | tr / ' '| awk '{print $3}')
		
	echo ----------------------
	echo "Working on $date"
	
	$scripts_dir/sort_frag.sh $fr_dir

	if [ $? -eq 0 ]
	then
		good_dir=$fr_dir/good
		
		cd $good_dir
		
		$scripts_dir/sort_sn.sh
 		
 		$scripts_dir/sort_correlation.sh
 		
 		$scripts_dir/sort_reg_I.sh
 		
 		
		
		
	fi
done < 	data_dirs.txt



echo "
Job submission is complete.  However, for unknown reasons,
the job dependencies sometimes get scrambled, causing
execution to hang and preventing all analysis from completing.

If there are still jobs being held in the queue, please run the
script a second or even third time to ensure that all analysis
has run.

Re-running of this script is still recommended to ensure the
completion of boots.

To complete the manual selection portion of the trajectory selection
process run the matlab script:
select_traj.m


Usefull commands:
qstat -n1 	(check the progress of jobs in the scheduler)
qdel \$JOBID 	(cancel a job)
"
