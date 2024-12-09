###두 분자의 좌표를 비교하여 RMSD (Root Mean Square Deviation)를 계산
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.spatial import distance_matrix
from scipy.optimize import linear_sum_assignment
import glob
import csv

def read_pdb(file):
    coords = []
    with open(file, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

def align_coordinates(coords1, coords2):
    distance_mat = distance_matrix(coords1, coords2)
    row_ind, col_ind = linear_sum_assignment(distance_mat)
    aligned_coords1 = coords1[row_ind]
    aligned_coords2 = coords2[col_ind]
    return aligned_coords1, aligned_coords2

def rmsd(aligned_coords1, aligned_coords2):
    return np.sqrt(np.mean(np.sum((aligned_coords1 - aligned_coords2) ** 2, axis=1)))

file1 = '2pvu_ligand.pdb'
file_pattern = '*gpu_out*.pdb'
output_csv_file = 'rmsd_results.csv'

coords1 = read_pdb(file1)

opt_dock_files = glob.glob(file_pattern)

rmsd_results = []

for file2 in opt_dock_files:
    coords2 = read_pdb(file2)
    aligned_coords1, aligned_coords2 = align_coordinates(coords1, coords2)
    rmsd_value = rmsd(aligned_coords1, aligned_coords2)
    print(f'RMSD between {file1} and {file2}: {rmsd_value:.4f} Å')

    rmsd_results.append((file1, file2, rmsd_value))

with open(output_csv_file, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(['File 1', 'File 2', 'RMSD (Å)'])
    csv_writer.writerows(rmsd_results)


# In[ ]:
