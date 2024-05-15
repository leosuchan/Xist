#!/bin/bash

START=$(date +%s.%N)

cd ~/Chaco-2.2/exec
if [ -d ~/Chaco-2.2/exec/ChacoOutput ]; then rm -rf ChacoOutput; fi
mkdir -p ChacoOutput
for ((i = 0; i<$2; i += 1));
do
	printf $1'\ncout'$i'.txt\n1\n500\n1\nn' > chacoinput.txt
	./chaco < chacoinput.txt
	mv "cout${i}.txt" ./ChacoOutput
	printf 'PROMPT=FALSE\nOUTPUT_TIME=1\nOUTPUT_ASSIGN=TRUE\nPRINT_HEADERS=FALSE\nKL_IMBALANCE='$(bc <<< "scale=2;(($i/100))") > User_Params
done

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
printf '\nTotal execution time: '$DIFF' seconds\n'
