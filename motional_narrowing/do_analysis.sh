#!/bin/bash


export VERSION=5.1

echo Alpha analysis - $VERSION
echo Jeffery Hanson - 20090327


global_vars="sample alphas work_dir boot_dir boot_num final_dir data_dirs scripts_dir output_file r0 data_archive nt sigma corbinIII corbinII log_dir reduce_pc alpha_dir sn_threshold"
./config.alpha
for G in $global_vars
do
	export $G="$(cat .$G)"
	
	#test to make sure variables are set
	if [ -z "$(cat .$G)" ]
	then
		echo
		echo "The variable $G is unset"
		echo "Make sure the file config.alpha"
		echo "is present and all values are set"
		echo "perhaps your config.alpha file is old?"
		echo "rerun install_analysis.sh to update it."
		
		exit 1
	fi
done


if [ -e .skip_data_selection ]
then
	export skip_data_selection=1
fi	

$scripts_dir/version.sh
if [ $? -eq 1 ]
then
	exit 0
fi	

if [ -z ${1} ]
then
	input=none
else
	input=${1}
fi

case "$input" in 
	alpha )
 		$scripts_dir/alpha_analysis.sh
	;;

	bootstrap )
 		$scripts_dir/bootstrap_analysis.sh	
	;;
	
	frsort )
		$scripts_dir/sort_frsort.sh
	;;

	clean_frag )
		$scripts_dir/clean_frag.sh
	;;
	
	clean_local)
		$scripts_dir/clean_local.sh
	;;	

	clean_bootstrap)
		$scripts_dir/clean_bootstrap.sh
	;;

	clean_bootsall)
		$scripts_dir/clean_bootsall.sh
	;;

	version )
		export print_version=1
		$scripts_dir/version.sh
	;;
	
	save )
		$scripts_dir/save_results.sh
	;;

				
	
	* )
		$scripts_dir/usage.sh
	;;
esac



