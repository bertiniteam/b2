#! /bin/bash

if [ $# -ne 1 ]
then
    numprocs=1
else
	numprocs=$1
fi

printf "\n\n\ncompiling with %d processor(s)\n\n\n" $numprocs | tee -a "compile_all.log"
while read suite; do
	printf "\n\n%s\n\n" $suite | tee -a "compile_all.log"
	make -j $numprocs $suite | tee -a "compile_all.log"
done <test/available_tests.txt
