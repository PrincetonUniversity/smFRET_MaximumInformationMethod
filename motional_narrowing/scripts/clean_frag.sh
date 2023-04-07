#!/bin/bash


if [ ! -e data_dirs.txt ]
then
	echo
	echo "Can't find the file: data_dirs.txt"
	echo "This file should contain the path to"
	echo "each directory containing the primary"
	echo "fret data to be analyzed, exiting ..."
	
	exit 1
fi


echo "Cleaning in data directory"

while read fr_dir
do	
	if [ -e $work_dir/job_id_list ]
	then
		rm $work_dir/job_id_list
	fi
	
	date=$(echo $fr_dir | tr / ' '| awk '{print $3}')
		
	echo ----------------------
	echo "Cleaning on $date"
	
	cd $fr_dir
        rm -rf good
	rm -rf bad
	rm *.log
	rm *.sam
	rm *.xc
	rm gmon.out
	rm frsort_out
	rm frsort_error	
done < data_dirs.txt

echo ----------------------
echo "Cleaning frag files done."
