#!/bin/bash
#

echo
echo Usage:
echo 'do_analysis.sh frsort		- autosort trajectories'
echo 'do_analysis.sh alpha 		- run alpha analysis'
echo 'do_analysis.sh bootstrap 	- do bootstraping'
echo 'do_analysis.sh version		- print version information'
echo 'do_analysis.sh save		- save current results to a directory'
echo 'do_analysis.sh clean_frag	- clean up all the analysis files in raw data folder'
echo 'do_analysis.sh clean_local	- clean up all the analysis files in home local folder'
echo 'do_analysis.sh clean_bootstrap	- clean up analysis results in runs but retain selected data'
echo 'do_analysis.sh clean_bootsall	- remove all boot_dir folder'
echo
echo
echo "all parameters for the analysis are set in the file:"
echo config.alpha
echo
