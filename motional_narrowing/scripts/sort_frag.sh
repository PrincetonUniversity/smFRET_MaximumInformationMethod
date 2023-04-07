
fr_dir=$1
date=$(echo $fr_dir | tr / ' '| awk '{print $3}')

function fr_count () {

	fr_cnt=0
	
	cnt_dir=$1
	
	for F in $cnt_dir/*.fr
	do
		let fr_cnt++
	done
	echo $fr_cnt
}
function fr_count_m () {

	fr_cnt=0
	
	cnt_dir=$1
	
	for F in $cnt_dir/*.m?.*fr
	do
		let fr_cnt++
	done
	echo $fr_cnt
}
total_fr=$(fr_count $fr_dir)

total_fr_m=$(fr_count_m $fr_dir)

if [ $total_fr -eq 0 ]
then
	echo 'No fret files were found in the directory:'
	echo $fr_dir
	echo Skipping...
	exit 0
else
	echo $total_fr fret files found
fi

if [ ! $total_fr -eq $total_fr_m ]
then
	
	let bad=total_fr-total_fr_m

	echo
	echo "$bad of the files in this directory are misnamed"
	echo "this will prevent them from being included in subsequent analysis"
	echo
	echo "each filename must include the microscope on which the data was collected"
	echo "somesample.m1.fr"
	echo 
fi
	

if [ -d $fr_dir/good ]
then
	good_fr=$(fr_count $fr_dir/good )
	bad_fr=$(fr_count $fr_dir/bad )
	
	let check_fr=good_fr+bad_fr
	
 	if [ $check_fr -eq $total_fr ]
 	then
 		exit 0
 	else
 		let total_fr=total_fr-good_fr
 		let total_fr=total_fr-bad_fr
 	fi
fi

alpha=10
QSUB_OPTS="-d $fr_dir -e frsort_error -o frsort_out -m n -M nomail"

if [ ! -d $fr_dir/good ]
then
	if [ ! -w $fr_dir ]
	then
		echo "You do not have permission to write to the directory:"
		echo $fr_dir
		echo "Have the owner or sysadmin change the permissions..."
		echo "chmod 777 $fr_dir"
		
		exit 1
		
	fi

	mkdir $fr_dir/good
	chmod 777 $fr_dir/good
else
	if [ ! -w $fr_dir/good ]
	then
		echo "You don't have permissions to write to the directory"
		echo $fr_dir/good
		echo "Have the owner or sysadmin change the permissions..."
		echo "chmod 777 $fr_dir/good $fr_dir/bad"
		
		exit 1
	fi	
fi


if [ ! -e $fr_dir/bad ]
then
	mkdir $fr_dir/bad
	chmod 777 $fr_dir/bad
fi


if [ -e $final_dir/.crosstalk ]
then
	XT_SWITCH="-x$(cat $final_dir/.crosstalk)"

else
	XT_SWITCH=''
fi


cnt=1
total_cnt=1
cd $fr_dir

num_jobs=20 # number of frag jobs to submit to a single pbs job
submit_job_cnt=1
submit_string=''

for F in *.fr
do

	JOBNAME=fr-$date-$total_cnt

	if [[  $total_cnt -eq $total_fr ]]
	then
		id_list=$work_dir/job_id_list_master
		source $scripts_dir/get_job_id_list.sh
		cnt=1
		submit_job_cnt=$num_jobs
	else
		id_list=$work_dir/job_id_list
		depend_switch=''
	fi
		
	QSUB_OPTS_F="$QSUB_OPTS -N$JOBNAME $depend_switch"	
		

	if [[ ! -e $fr_dir/good/$F && ! -e $fr_dir/bad/$F ]]
	then
	
		submit_string="$submit_string $F"
		let submit_job_cnt++
	fi	
	
	
	if [ $submit_job_cnt -ge $num_jobs ]
	then
	
		echo "
			for FRAG in $submit_string
			do
				frag \${FRAG} -a.$alpha $XT_SWITCH
				if [ \${?} -eq 1 ]
				then
					ln \${FRAG%%.fr}.* bad
				else
					ln \${FRAG%%.fr}.* good
				fi
			done	
			sleep 15 
		" | qsub $QSUB_OPTS_F  >> $id_list
		let cnt++
		submit_job_cnt=1
		submit_string=''
			
	fi
	
	let total_cnt++
done


mv $work_dir/job_id_list_master $work_dir/job_id_list 2>/dev/null

exit 0
