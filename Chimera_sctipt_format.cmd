open 1onp
select
~select ligand
delete sel
write format mol2 #0 1onp_Nlig.mol2
select :.A
select invert
delete sel
addh
addcharge std 
addcharge nonstd #0 0 method gas
write format mol2 #0 1onp_Nlig_prep.mol2
