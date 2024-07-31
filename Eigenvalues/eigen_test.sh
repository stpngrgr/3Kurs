#!/bin/bash

if [ -n "$1" ]
then
	fl=1
	echo "-------------------------flag---------------------------"
else
	fl=0
fi

for (( s = 1; s < 5; s++ ))
do \
	for (( n = 1; n <= 100; n++ ))
	do \
		echo "================ n = ${n} s = ${s} =================="
		result=`./a.out ${n} 10 1.e-14 ${s} | tr -d '\0'`\

		if [ $fl -eq 1 ]
		then
			error=`echo "${result}" | grep -a 'ERROR'`
			residual1=`echo "${result}" | grep -a 'Residual1'`
			residual2=`echo "${result}" | grep -a 'Residual2'`
			iterations=`echo "${result}" | grep -a 'Iterations'`
			iterations1=`echo "${result}" | grep -a 'Iterations1'`
			time1=`echo "${result}" | grep -a 'Elapsed1'`
			time2=`echo "${result}" | grep -a 'Elapsed2'`
			echo "${error}"
			echo "${residual1}"
			#echo "${residual2}"
			#echo "${iterations}"
			#echo "${iterations1}"
			#echo "${time1}"
			#echo "${time2}"
		else
			echo "${result}"
		fi
		echo "-----------------------------------------------------"
	done
done


