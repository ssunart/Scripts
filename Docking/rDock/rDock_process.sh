#!/bin/bash

###prep with Chimera 'Dock prep' & obabel###
###prepare preprocessed receptor.mol2 & reference ligand.sdf###

rbcavity -was -d -r <PRMFILE>

###transfer lig with OBABEL###
obabel -i sdf 1onp_Lig_0.sdf -o sdf -O 1onp_Lig_0.sdf


rbdock -i <INPUT>.sdf -o <OUTPUT> -r <PRMFILE> -p dock.prm -n 50
