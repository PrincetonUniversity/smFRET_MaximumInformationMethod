#!/bin/bash


which bc 2>/dev/null 1>/dev/null

if [ ! $? -eq 0 ]
then
	echo "This script requires the program \"bc\""
	echo "Please ask Jeff or Dan to install it on this machine."
	exit
fi	

echo
echo "Statistics for daily data collection:"
echo

function fr_count () {

	fr_cnt=0
	
	cnt_dir=$1
	
	for F in "$cnt_dir/"*.fr
	do
		if [ -e "$F" ]
		then
			let fr_cnt++
		fi	
	done
	
	echo $fr_cnt
}


function fr_count_m () {

	fr_cnt=0
	
	cnt_dir=$1
	
	M=$2
	
	for F in $cnt_dir/*$M*.fr
	do
		if [ -e "$F" ]
		then
			let fr_cnt++
		fi	
	done
	
	echo $fr_cnt
}

printf "%8s | %6s | %6s | %6s | %6s | %6s | %6s\n" '' 'total' 'good' 's/n' 'corr' 'manual' 'Eff'

log_dir=$(cat .log_dir)

while read fr_dir
do
	date=$(echo $fr_dir | tr / ' '| awk '{print $3}')
	
	if [ $date == $USER ]
	then
		date=archive
	fi	
			
	
	printf "%8s |-%6s-|-%6s-|-%6s-|-%6s-|-%6s-|-%6s-\n" $date '------' '------' '------' '------' '------' '------'
		
		
	sample='Total'	
		
	total_n=$(fr_count $fr_dir)
	good_n=$(fr_count $fr_dir/good/)
	
	sn_n=$(cat $fr_dir/good/pass_sn.txt 2>/dev/null| wc -l)
	 
	cor_n=$(cat $fr_dir/good/pass_corr.txt 2>/dev/null| wc -l)
	
	cor_fn=$(cat $fr_dir/good/fail_corr.txt 2>/dev/null| wc -l)
	
	let "cor_total=cor_n + cor_fn"
	
	man_n=$(cat $log_dir/${date}_pass_manual.txt 2>/dev/null| wc -l)
	
	if [ ! $total_n -eq 0 ]
	then
		eff=$( echo "scale=9; $man_n/$total_n*100" | bc )
	else
		eff=0
	fi
	
	if [ $cor_total -eq $sn_n ]
	then
		printf "%8s | %6s | %6s | %6s | %6s | %6s | %4.1f %%\n" $sample $total_n $good_n $sn_n $cor_n $man_n $eff
	else
		printf "%8s | %6s | %6s | %6s | %3s%3s | %6s | %4.1f %%\n" $sample $total_n $good_n $sn_n 'XXX' $cor_n $man_n $eff
	fi	
	for M in m{1..4}
	do
		total_n=$(fr_count_m $fr_dir $M)
		good_n=$(fr_count_m $fr_dir/good/ $M)
		
		sn_n=$(cat $fr_dir/good/pass_sn.txt 2>/dev/null | grep $M | wc -l)
	 
		cor_n=$(cat $fr_dir/good/pass_corr.txt 2>/dev/null | grep $M | wc -l)
		
		man_n=$(cat $log_dir/${date}_pass_manual.txt 2>/dev/null | grep $M | wc -l)
		
		if [ ! $total_n -eq 0 ]
		then		
			eff=$( echo "scale=9; $man_n/$total_n*100" | bc )
		else
			eff=0
		fi
		
		printf "%8s | %6s | %6s | %6s | %6s | %6s | %4.1f %%\n" $M $total_n $good_n $sn_n $cor_n $man_n $eff
	done
	

done < data_dirs.txt


ls data/m?/prune 1>/dev/null 2>/dev/null

if [ ! $? -eq 0 ]
then
	exit
fi
	
echo
echo "Statistics for ioa/iod ratio test:"
echo
printf "%8s | %6s| %6s\n" '' 'Initial' 'Final'
printf "%8s-|-%6s-|-%6s-\n" '--------' '------' '------'
for M in m{1..4}
do
	good_n=$(fr_count data/$M/)
	bad_n=$(ls -l data/$M/prune/*.fr 2>/dev/null | wc -l)
	
	let "total_n= good_n+ bad_n"
	
	let all_good+=good_n
	let all_total+=total_n
	
	if [ $total_n -gt 0 ]
	then
		printf "%8s | %6s | %6s\n" $M $total_n $good_n
	fi
done
printf "%8s-|-%6s-|-%6s-\n" '--------' '------' '------'
printf "%8s | %6s | %6s\n" 'Total' $all_total $all_good
