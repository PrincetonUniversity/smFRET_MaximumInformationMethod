#!/bin/bash

LATEST_VERSION='5.1'

if [ "$LATEST_VERSION" != "$VERSION" ]
then
	echo
	echo Your analysis files may be out of date...
	echo 
	echo They have been automatically updated.
	echo
	
	export alpha_dir="$(cat .alpha_dir)"

	CONFIG_DIFF=$(diff config.alpha $alpha_dir/config.alpha)
	
	if [ $? -eq 1 ]
	then
		echo "Your config.alpha file has been overwritten"
		echo "Here is the out put of the diff command"
		echo "\"diff old_config new_config\":"
		echo
		echo $CONFIG_DIFF
		echo
		echo "Please update your new config.alpha file accordingly"
	fi
	echo	
	echo "Rerun the \"./do analysis\" command to continue"
	install_analysis.sh 1>/dev/null
	
	exit 1

fi


if [ -z $print_version ]
then
	exit 0
fi	


