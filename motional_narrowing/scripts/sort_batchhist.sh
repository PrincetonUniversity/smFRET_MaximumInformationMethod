#!/bin/bash


alpha=10 	
boot_single=3	
boot_all=5	

wd=$(pwd)

QSUB_OPTS_ALL="-d $wd -o boots_out_all -e boots_error_all -l walltime=7200 -m n -V -M nomail  $depend_switch"
QSUB_OPTS_SINGLE="-d $wd -o boots_out_single -e boots_error_single -l walltime=7200 -m n -V -M nomail  $depend_switch"


if [ -e $final_dir/.crosstalk ]
then
	XT_SWITCH="-x$(cat $final_dir/.crosstalk)"

else
	XT_SWITCH=''
fi

OUT_F=hist$alpha.out
OPTS="-a0.$alpha -m -b$boot_all -o$OUT_F $XT_SWITCH"

NAME=" -N boots_all"

echo "
		 
	 	
	if [ ! -e $OUT_F ]
	then

			boots $OPTS \$(cat pass_corr.txt) > boots_all

	fi
" | qsub $QSUB_OPTS_ALL $NAME >> /dev/null



OPTS="-a0.$alpha -m -b$boot_single $XT_SWITCH"

NAME=" -N boots"
echo "

	
	for F in \$(cat pass_corr.txt)
	do
		if [ ! -e \${F%%.fr}.hist ]
		then

				boots $OPTS \$F  > boots_single

		fi
	
	done
" | qsub $QSUB_OPTS_SINGLE $NAME>> /dev/null
