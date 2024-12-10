###prep_vina config file
import os
import re
from glob import glob

# Step 1: Check if the vinaR directory exists, if not, create it
if not os.path.exists('vinaR'):
    os.makedirs('vinaR')

# Step 2: Find the PDBID from the filename in the current directory
validation_files = glob("*Validation_Affinities*.sdf")

for file in validation_files:
    # Extract PDBID from the filename (assuming the pattern is 'setXXX_PDBID...Validation_Affinities...')
    pdbid_match = re.search(r'_([A-Z0-9]{4})-[\w]+_Validation_Affinities', file)
    
    if pdbid_match:
        pdbid = pdbid_match.group(1)
    else:
        raise ValueError(f"PDBID not found in filename: {file}")
    
    # Step 3: Locate the corresponding gpf file and extract grid center
    gpf_file = f"{pdbid}_prep.gpf"
    
    if os.path.exists(gpf_file):
        with open(gpf_file, 'r') as gpf:
            content = gpf.read()
            gridcenter_match = re.search(r'gridcenter\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)', content)
            if gridcenter_match:
                center_x = gridcenter_match.group(1)
                center_y = gridcenter_match.group(2)
                center_z = gridcenter_match.group(3)
            else:
                raise ValueError(f"Grid center not found in {gpf_file}")
    
    # Step 4: Create config file for each ligand file in ligands directory
    ligand_files = glob('ligands/*.pdbqt')
    
    for ligand_file in ligand_files:
        ligand_base = os.path.basename(ligand_file)
        
        # Prepare config file content
        config_content = f"""
receptor = {pdbid}_prep.pdbqt
ligand = ligands/{ligand_base}
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 32
num_modes = 100
min_rmsd = 5
energy_range = 5
out = vinaR/{ligand_base.replace('.pdbqt', '_vinaR.pdbqt')}
"""
        
        # Write the config file to the vinaR directory
        config_filename = f"vinaR/config_{ligand_base.replace('.pdbqt', '')}.txt"
        with open(config_filename, 'w') as config_file:
            config_file.write(config_content.strip())

print("Config files have been generated.")
