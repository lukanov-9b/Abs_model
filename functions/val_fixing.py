from rdkit import Chem
# исправление валентностей на атомах бора и азота
def ValFixing(smiles:str) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        errors = Chem.DetectChemistryProblems(mol)
        if not errors:
            return Chem.CanonSmiles(smiles)
        else:
            for error in errors:
                if error.GetType()=='AtomValenceException':
                    atom = mol.GetAtomWithIdx(error.GetAtomIdx())
                    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge()==0 and atom.GetExplicitValence()==4:
                        atom.SetFormalCharge(1)
                    elif atom.GetAtomicNum() == 5 and atom.GetFormalCharge()==0 and atom.GetExplicitValence()==4:
                        atom.SetFormalCharge(-1)
            
            Chem.SanitizeMol(mol)
            return Chem.CanonSmiles(Chem.MolToSmiles(mol))
    except:
        return None