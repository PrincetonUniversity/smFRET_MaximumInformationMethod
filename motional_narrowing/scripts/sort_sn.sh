

source $scripts_dir/get_job_id_list.sh

echo $sn_threshold > .sn_threshold
chmod 770 .sn_threshold

QSUB_OPTS="-d $(pwd) -o sn_output.out -e sn_error.out -m n -M nomail -N sort_sn $depend_switch"

echo "
	pass_file=pass_sn.txt
	sn_threshold=\$(cat .sn_threshold)

	if [ -e \$pass_file ]
	then
		rm \$pass_file
	fi
		
	touch \$pass_file
	
	for F in *.fr
	do
		log=\${F%%.fr}.log
		dsn=\$(cat \$log | sed -n 4p | tr . ' ' |awk '{print \$1}')
		asn=\$(cat \$log | sed -n 5p | tr . ' ' |awk '{print \$1}')
		
		if [[ \$dsn -ge \$sn_threshold && \$asn -ge \$sn_threshold ]]
		then
			echo \$F >> \$pass_file
		fi
	done
	chmod 770 \$pass_file
	" | qsub $QSUB_OPTS >> $work_dir/job_id_list
