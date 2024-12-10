###VinaGPU 작업 자동화
import os
import re
import subprocess
from glob import glob

# Define base directory for your project and the Vina GPU executable path
base_dir = "/home/soowon/project/01_ConsensusDocking/DockingBenchmark/BindingDB/part1"
vina_gpu_dir = "/home/soowon/apps/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1"

# Step 1: Check if the VGPUR directory exists, if not, create it
vgpur_dir = os.path.join(base_dir, 'VGPUR')
if not os.path.exists(vgpur_dir):
    os.makedirs(vgpur_dir)

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
    gpf_file = os.path.join(base_dir, f"{pdbid}_prep.gpf")
    
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
    ligands_dir = os.path.join(base_dir, 'ligands/pdbqt')
    ligand_files = glob(os.path.join(ligands_dir, '*.pdbqt'))
    
    for ligand_file in ligand_files:
        ligand_base = os.path.basename(ligand_file)
        
        # Prepare config file content with absolute paths
        config_content = f"""
receptor = {os.path.join(base_dir, f'{pdbid}_prep.pdbqt')}
ligand = {ligands_file}
opencl_binary_path = /home/soowon/apps/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}
size_x = 20.0
size_y = 20.0
size_z = 20.0
thread = 8000
num_modes = 100
energy_range = 5
output_directory = {vgpur_dir}
"""
        
        # Write the config file to the VGPUR directory with an absolute path
        config_filename = os.path.join(vgpur_dir, f"config_{ligand_base.replace('.pdbqt', '')}.txt")
        with open(config_filename, 'w') as config_file:
            config_file.write(config_content.strip())
        
        # Step 5: Run the VGPU command from the specified directory
        VGPU_command = f"./AutoDock-Vina-GPU-2-1 --config {config_filename}"
        print(f"Running: {VGPU_command} from {vina_gpu_dir}")
        
        # Change the working directory to the Vina-GPU directory
        original_dir = os.getcwd()  # Save the current working directory
        os.chdir(vina_gpu_dir)  # Change to the Vina-GPU directory
        
        try:
            # Execute the command
            result = subprocess.run(VGPU_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(result.stdout.decode())
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running {VGPU_command}:")
            print(e.stderr.decode())
        
        # Change back to the original directory
        os.chdir(original_dir)

