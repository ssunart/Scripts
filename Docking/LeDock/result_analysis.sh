#!/bin/bash

while read line
do
    echo $line
    cp ledock_anal* $line/
    cd /home/soowon/Project/Docking_benchmark_2023/LeDock/$line
    ./ledock_anal.csh
    cd ..
done < dir_list.txt
