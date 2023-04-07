#!/bin/bash


cd $final_dir


if [ ! -e ${output_file}.mat ]
then
	echo Cant find file ${output_file}.mat
	echo Finish alpha analysis before saving data
	exit 1
fi

for A in $alphas
do
	if [ ! -e ${sample}${A}.out ]
	then
		echo Cant find file ${sample}${A}.out
		echo Finish alpha analysis before saving data
		exit 1
	fi
done


date=`date +%Y%m%d`
timeout=30
echo
echo "Input Directory name for saved data, without spaces (default = $date):"
read -t $timeout save_dir

if [ -z $save_dir ]
then
	save_dir=save_$date
else
	save_dir=save_$save_dir
fi


if [ -d $save_dir ]
then
	echo the directory:
	echo
	echo $(pwd)/$save_dir 
	echo
	echo already exists ... choose a different name
	exit 1
else
	echo saving current results to: 
	echo
	echo $(pwd)/$save_dir
	echo
fi

mkdir $save_dir


cp *.out $save_dir
cp *.mat $save_dir
cp *.m $save_dir
cp data_dirs.txt $save_dir

if [ ! -z $log_dir ]
then
	cp -r $log_dir $save_dir/
fi	


cp config.alpha $save_dir/
cd $save_dir
./config.alpha

echo $(pwd) > .final_dir
echo "no_saving_bootstrap" > .boot_dir
cp ../.sample .
cp ../.output_file .