#!/bin/bash



cd $work_dir


for A in $alphas
do
	if [ ! -d $A ]
	then
		mkdir $A
	fi
done


if [ -n "$skip_data_selection" ]
then
	test_dir_num=1
	while read D
	do
		if [ $test_dir_num -gt 1 ]
		then
			# print a warning if data_dirs.txt has more than one line
			echo
			echo "trajectory selection has been skipped"
			echo "only a single line of the data_dirs.txt file will be used"
			echo "please update your file to reflect this..." 
			echo "ie put all data in a single directory"
			echo
			break
		fi
			
		data_dirs="$D"
		let test_dir_num++			
	done < 	$final_dir/data_dirs.txt
fi


for D in $data_dirs
do
	#test for existance of directories
	if [ ! -d $D ]
	then
		echo 
		echo No directories found matching the criteria:
		echo $D
		echo "Check the config.alpha file"
		echo
		
		exit 1
	fi
	
	for A in $alphas
	do
		cd $work_dir/$A
		for F in $D/*.fr
		do
			F_BASE=$(echo $F | sed s%$D/%%)
			if [ ! -e $F_BASE ]
			then
				ln -s $F .
			fi
		done
	done
done


for A in $alphas
do
	cd $work_dir/$A
	for F in *.fr
	do
		if [ ! -e $F ]
		then
			for EXT in fr xc log sam hist
			do
				rm -f ${F%%.fr}.$EXT >> /dev/null
			done
		fi	
	done
done

num=$(ls -l *.fr|wc -l)
echo Found $num fret trajectories

exit 0
