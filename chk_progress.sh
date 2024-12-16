#!/bin/bash
clear

date
cd ~/project/cache_challenge/suyeon/vina_adt

echo -e "\n----- Suyeon Cache_challenge Vina progress ----- \n"

echo `ls -v | wc -l`"/452414"


cd ~/project/cache_challenge/suyeon/results

echo -e "\n----- Suyeon Cache_challenge AD_GPU progress --- \n"
ls -v > fin_list 

echo `grep -r 'dlg' fin_list | wc -l`"/452414"



cd ~/project/ConsensusDocking_ML/part1/PDBbind/refined_set

echo -e "\n----- PDBbind AD_GPU bigbox progress ----------- \n"

echo `ls -d -v *_big | wc -l`"/5316"
echo -e "\n------------------------------------------------ \n"

