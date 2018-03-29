#!/bin/bash

for number in {0..359}
do
./main $number 0
done

for number in {360..719}
do
./main $number 1
done
exit 0
