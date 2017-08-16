#! /bin/bash

if [ $# -ne 1 ]
then
    numprocs=1
else
	numprocs=$1
fi

echo "compiling with" $numprocs "processor(s)" | tee -a "compile_all.log"
while read suite; do
	echo $suite | tee -a "compile_all.log"
	make -j $numprocs $suite | tee -a "compile_all.log"
done <test/available_tests.txt
