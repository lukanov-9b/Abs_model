from rdkit import Chem


def get_formula(smiles) -> str:
    atoms = {}
    try:
        for atom in Chem.MolFromSmiles(smiles).GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in atoms:
                atoms[symbol] = 1
            else:
                atoms[symbol] += 1
        return ''.join([k + str(v) for k,v in atoms.items()])
    except:
        return None


def get_element_composition(smiles: str) -> str:
    element_list = set([x.GetSymbol() for x in Chem.MolFromSmiles(smiles).GetAtoms()])
    element_list = list(element_list)
    element_list.sort()

    return ''.join(element_list)