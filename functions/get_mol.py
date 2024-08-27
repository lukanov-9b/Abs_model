import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pymatgen.core import Structure, Lattice


def SmilesToMol(smi: str):
    mHs = Chem.AddHs(Chem.MolFromSmiles(smi))
    # start generating 3D coordinates and optimize the conformation
    embedError = AllChem.EmbedMolecule(mHs, useRandomCoords=True)
    if embedError == 0:
        UffoptError = AllChem.UFFOptimizeMolecule(mHs, 3000)
        if UffoptError == 0:
            return Chem.MolToMolBlock(mHs)

        print("UFF optimization failed, trying MMFF optimization")
        MMFFoptError = AllChem.MMFFOptimizeMolecule(mHs, 3000)
        if MMFFoptError == 0:
            return Chem.MolToMolBlock(mHs)

        print("MMFF optimizaiton has failed")
        return None
    else:
        print("Embedding Failed")
        return None


def get_structure(smiles: str) -> Structure:
    try:
        bias = 5
        lattice = Lattice([[50, 0, 0], [0, 50, 0], [0, 0, 50]])
        mol = SmilesToMol(smiles)
        pattern = '([-\.\d]+)\s+([-\.\d]+)\s+([-\.\d]+)\s+([A-Za-z]+)\s+\d+'
        atoms = re.findall(pattern, mol)
        coords = np.array([[float(j) for j in [x[0], x[1], x[2]]] for x in atoms])
        coords /= 50
        h_x = abs(coords[:, 0].min())
        h_y = abs(coords[:, 1].min())
        h_z = abs(coords[:, 2].min())
        coords = coords + np.array([h_x, h_y, h_z]) + bias
        sites = [x[3] for x in atoms]

        return Structure(lattice, sites, coords)
    except:
        return None
