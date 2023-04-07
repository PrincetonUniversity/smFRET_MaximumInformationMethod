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


if [ -n "$skip_data_selection" ]
then
	echo "This script does not need to be run"
	echo "trajectory selection has been skipped in config.alpha"
	echo
	
	exit 1
fi	
	

echo "Cleaning Locally"
rm *.jpg
rm *.out
rm *.m
rm *.mat
rm host_info
rm svn-commit.tmp
rm -rf data/
rm *.txt
rm config.alpha

