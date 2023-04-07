#!/bin/bash


if [ -d $boot_dir ]
then	
	echo "$boot_dir is being removed."
	rm -r $boot_dir
else
	echo -e "no $boot_dir found. Exit.\n"
fi

