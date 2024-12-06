###glg 파일을 읽어 단백질 전체를 grid box로 지정하는 global docking 수행

#!/bin/bash

###Big size preprep###



###part마다 기본값 AD-GPU 도킹 수행하기###
###PART 수정 필요###
cd part1


while read line
do
    echo "$line"
    cd "$line"_big
    mkdir gpuR
##Grid prepare for Big##
    cp ../"$line"/ligands/"$line"_ligand.pdbqt ligands/ 
    rm midpoint
    rm sum.txt

    cat "$line"_prep.glg | grep "Midpoint" | cut -d "," -f1 | head -1 | tr -d " |  |___|__| Midpoint = (" >> midpoint
    cat "$line"_prep.glg | grep "Midpoint" | cut -d "," -f2  | head -1 | tr -d " " >> midpoint

    cat "$line"_prep.glg | grep "Midpoint" | cut -d "," -f3  | head -1 | tr -d " )" >> midpoint

    sed -e "1 s/^/center_x = /" midpoint | sed -e "2s/^/center_y = /" | sed -e "3s/^/center_z = /" > config.txt
 
    cat "$line"_prep.glg | grep "Maximum coordinates" | cut -d "," -f1 | head -1 | tr -d "Maximum coordinates :           (" | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt
    cat "$line"_prep.glg | grep "Maximum coordinates" | cut -d "," -f2 | head -1 | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt
    cat "$line"_prep.glg | grep "Maximum coordinates" | cut -d "," -f3 | head -1 | tr -d ")" | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt

    cat "$line"_prep.glg | grep "Minimum coordinates" | cut -d "," -f1 | head -1 | tr -d "Maximum coordinates :           (" | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt
    cat "$line"_prep.glg | grep "Minimum coordinates" | cut -d "," -f2 | head -1 | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt
    cat "$line"_prep.glg | grep "Minimum coordinates" | cut -d "," -f3 | head -1 | tr -d ")" | awk '{ sub(/^[ \t]+/, ""); print }' >> sum.txt

    cat sum.txt | tr -d '-' > sum_val.txt
    
    IFS=$'\n' values=(`cat sum_val.txt`)
#for VALUE in "${values[@]}"; do echo "<---- $VALUE ---->"; done
#echo ${values[0]}
#echo ${values[3]}
    result_a=$(echo "scale=4; ${values[0]} + ${values[3]} +2" | bc)
    echo size_x = "$result_a" >> config.txt
    result_b=$(echo "scale=4; ${values[1]} + ${values[4]} +2" | bc)
    echo size_y = "$result_b" >> config.txt
    result_c=$(echo "scale=4; ${values[2]} + ${values[5]} +2" | bc)
    echo size_z = "$result_c" >> config.txt

    NPTS=$(echo $result_a,$result_b,$result_c)

    IFS=$'\n' values_2=(`cat midpoint`)

    echo ${values[0]}
    echo ${values_2[0]}
    prepare_gpf4.py -l ligands/"$line"_ligand.pdbqt -r "$line"_prep.pdbqt -p  ligand_types=HD,Br,Cl,A,C,N,P,NA,OA,F,S,SA,I  -p gridcenter=${values_2[0]},${values_2[1]},${values_2[2]} -p spacing=1.0
    sed -i "1s/.*/npts 50 50 50/g" "$line"_prep.gpf
###prepare_gpf4.py –l ligands/"$line"_ligand.pdbqt –r "$line"_prep.pdbqt –p ligand_types='HD,Br,Cl,A,C,N,P,NA,OA,F,S,SA,I' –p npts='"$result_a","$resul    t_b","$result_c"' –p spacing='1' -p gridcenter='${values_2[0]},${values_2[1],${values_2[2]}}' 2>&1
    ./autogrid4 -p "$line"_prep.gpf -l "$line"_prep.glg

##Run AD-GPU##
    autodock_gpu_128wi --lfile ligands/"$line"_ligand.pdbqt --ffile "$line"_prep.maps.fld --resnam gpuR/"$line"_gpuR --nrun 100
    cd gpuR/
    grep '^DOCKED' "$line"_gpuR.dlg | cut -c9- > "$line"_gpuR.pdbqt
    cd ../../
done < dir_list   
