#!/bin/bash
###Check Set NUMBER###
LIST_FILE_NAME=("set499_result")

echo $LIST_FILE_NAME

cd $LIST_FILE_NAME
ls *.mol2 > ligands
ledock_linux_x86 dock.in &

cp ../ledock_anal* .
./ledock_anal.csh
