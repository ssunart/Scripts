###mol2 파일로부터 분자의 질량 중심(중심 위치)을 계산
import numpy as np

def get_mass_center(mol2_file):
    atom_masses = {
        'H': 1.008,
        'C': 12.01,
        'N': 14.007,
        'O': 15.999,
        'F': 18.998
    }
    total_mass = 0
    mass_center = np.array([0.0, 0.0, 0.0])
    with open(mol2_file, 'r') as file:
        for line in file:
            if line.startswith('@<TRIPOS>ATOM'):
                break
        for line in file:
            if line.startswith('@<TRIPOS>BOND'):
                break
            parts = line.split()
            atom_type = parts[5][0]  # Assume the atom type is the first character of the 6th field
            if atom_type not in atom_masses:
                continue
            mass = atom_masses[atom_type]
            position = np.array([float(parts[2]), float(parts[3]), float(parts[4])])
            mass_center += mass * position
            total_mass += mass
    mass_center /= total_mass
    return mass_center

print(get_mass_center('z3n_ligand.mol2'))
