# Open Babel 예시
from openbabel import openbabel, pybel

smiles = 'COC1=CC=C(C=C1)C1=CC(=O)C2=C(O)C(OC)=C(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)C=C2O1'
mol = pybel.readstring('smi', smiles)
mol.make3D()
mol.write('sdf', 'Pectolinarin_obabel.sdf', overwrite=True)
mol.write('mol2', 'Pectolinarin_obabel.mol2', overwrite=True)

# RDKit 예시
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)
writer = Chem.SDWriter('Pectolinarin_rdkit.sdf')
writer.write(mol)
writer.close()
