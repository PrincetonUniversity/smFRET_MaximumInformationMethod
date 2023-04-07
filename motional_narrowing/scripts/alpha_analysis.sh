
#!/bin/bash


##########################################################
#Phase 1 - create directories and link data   ############
##########################################################
echo  =====================================================
echo Creating directories and gathering data...
echo

$scripts_dir/alpha_get_data.sh

if [ ${?} -eq 1 ]
then
	exit 1
fi


$scripts_dir/alpha_create_dir_structure.sh
if [ ${?} -eq 1 ]
then
	exit 1
fi

##########################################################
#Phase 2 - analyze data with frag and boots   ############
##########################################################

echo =====================================================
echo "Submitting batch boots and frag jobs to scheduler..."
echo

$scripts_dir/alpha_frag_boots.sh


if [ ${?} -eq 1 ]
then
	exit 1
fi

##########################################################
#Phase 3 - Using Matlab, calculate the time    ###########
# resolution of each alpha dataset             ###########
##########################################################

echo =====================================================
echo "Calculating mean time resolution for alpha sets"
echo 

$scripts_dir/alpha_gettr.sh

##########################################################
#Phase 4 - Using Matlab, fit data to the 2-state   #######
# motional narrowing model                          #######
##########################################################

echo =====================================================
echo "Fit to 2-state motional narrowing model"
echo

$scripts_dir/alpha_GevaSkinner.sh


