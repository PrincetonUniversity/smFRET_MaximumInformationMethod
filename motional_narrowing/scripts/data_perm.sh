while read D
do
        chmod a+rw $D
	chmod a+rw $D/good $D/bad
	chmod a+rw $D/good/.nt $D/good/.sigma $D/good/.corbinII $D/good/.corbinIII 2>/dev/null
done < data_dirs.txt
