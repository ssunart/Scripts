source /home/soowon/anaconda3/envs/openmm/amber.sh
ante-MMPBSA.py  -p /home/soowon/project/practice/pract_openmm/test_2/SYS_gaff2_nw.prmtop -c com.prmtop -r rec.prmtop -l ligand.prmtop -s :WAT:Na+:Cl-:Mg+:K+ -n :LIG --radii mbondi2
MMPBSA.py -O -i mmpbsa.in -o /home/soowon/project/practice/pract_openmm/test_2/FINAL_RESULTS_MMPBSA.dat -sp /home/soowon/project/practice/pract_openmm/test_2/SYS_gaff2_nw.prmtop -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y /home/soowon/project/practice/pract_openmm/test_2/prot_lig_prod1-2_nw.pdb
mkdir /home/soowon/project/practice/pract_openmm/test_2/MMPBSA_igb_2
mv _MMPBSA* com.prmtop rec.prmtop ligand.prmtop reference.frc mmpbsa.in /home/soowon/project/practice/pract_openmm/test_2/MMPBSA_igb_2
