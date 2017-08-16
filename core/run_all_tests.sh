#! /bin/bash


printf "running all tests\n\n\n" | tee -a "run_all.log"

while read suite; do
	printf "%s\n\n" $suite | tee -a "run_all.log"
	./$suite | tee -a "run_all.log"
done <test/available_tests.txt
