

fr_dir=$(pwd)


date=$(echo $fr_dir | tr / ' '| awk '{print $3}')


cp $scripts_dir/matlab/sort_corsort{'',.ctf} $fr_dir/ 2>/dev/null
cp $scripts_dir/matlab/run_sort_corsort.sh $fr_dir

if [ ! -d $log_dir ]
then
	mkdir $log_dir
fi	


if [ -e $fr_dir/.nt ]
then	
	delete=0
	files=".nt .sigma .corbinII .corbinIII"
	for F in $files
	do	
		diff $final_dir/$F $fr_dir/$F 1>/dev/null 2>/dev/null
		if [ ${?} -gt 0  ]
		then
			delete=1
		fi
	done	

	if [ $delete -eq 1 ]
	then
		rm $fr_dir/pass_corr.txt $fr_dir/fail_corr.txt 2>/dev/null
	else
		cp $log_dir/${date}_pass_corr.txt pass_corr.txt 2>/dev/null
		cp $log_dir/${date}_fail_corr.txt fail_corr.txt 2>/dev/null
	fi
	
else
	rm $fr_dir/pass_corr.txt $fr_dir/fail_corr.txt 2>/dev/null	
	
fi



cp $final_dir/.nt $final_dir/.sigma $final_dir/.corbinII $final_dir/.corbinIII $fr_dir/ 1>/dev/null 2>/dev/null
chmod 770 $fr_dir/.nt $fr_dir/.sigma $fr_dir/.corbinII $fr_dir/.corbinIII $fr_dir/ 1>/dev/null 2>/dev/null

source $scripts_dir/get_job_id_list.sh


QSUB_OPTS="-d $fr_dir  -o corr.out -e corr.error -m n -M nomail $depend_switch"


job_name="-N sortcorr"
echo "
	./run_sort_corsort.sh /usr/local/MATLAB/R2012a/ > corsort_output
	rm -r sort_corsort*
	chmod 770 $fr_dir/pass_corr.txt $fr_dir/fail_corr.txt
	cp pass_corr.txt $log_dir/${date}_pass_corr.txt
	cp fail_corr.txt $log_dir/${date}_fail_corr.txt
	" | qsub $QSUB_OPTS $job_name >> $work_dir/job_id_list
	
