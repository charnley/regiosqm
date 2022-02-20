from typing import List, Tuple

from ppqm import chembridge
from rdkit import Chem
from rdkit.Chem import AllChem, Mol


def get_target_atoms(molobj, target):
    """Find target atom indices from SMART"""
    atoms = molobj.GetSubstructMatches(target)
    # convert tuple of tuple to one-dimensional list
    atoms = [element for tupl in atoms for element in tupl]
    return atoms


def generate_protonations(mol: Mol) -> Tuple[List[int], List[Mol]]:
    """
    Protonate all aromatic carbons

    returns
        atom centers
        protonated molobjs per atom center

    """

    mol_prime = chembridge.copy_molobj(mol)

    Chem.Kekulize(mol_prime, clearAromaticFlags=True)

    reaction1 = AllChem.ReactionFromSmarts("[C;R;H1:1]=[C,N;R;H1:2]>>[CH2:1][*H+:2]")
    reaction2 = AllChem.ReactionFromSmarts("[C;R;H1:1]=[C,N;R;H0:2]>>[CH2:1][*+;H0:2]")

    molobjs = []
    target_atoms = []

    smarts_1 = Chem.MolFromSmarts("[C;R;H1:1]=[C,N;R;H1:2]")
    smarts_2 = Chem.MolFromSmarts("[C;R;H1:1]=[C,N;R;H0:2]")
    atoms_1 = get_target_atoms(mol_prime, smarts_1)
    atoms_2 = get_target_atoms(mol_prime, smarts_2)

    i = 0
    products_1 = reaction1.RunReactants((mol_prime,))
    for x in products_1:

        mol_ = x[0]
        smiles = Chem.MolToSmiles(mol_)
        smiles = smiles.replace("NH2+", "N+")
        mol_ = Chem.MolFromSmiles(smiles)

        molobjs.append(mol_)
        target_atoms.append(atoms_1[i])

        i += 1

    isav = i

    products_2 = reaction2.RunReactants((mol_prime,))
    for x in products_2:

        mol_ = x[0]
        smiles = Chem.MolToSmiles(mol_)
        smiles = smiles.replace("NH2+", "N+")
        mol_ = Chem.MolFromSmiles(smiles)

        molobjs.append(mol_)
        target_atoms.append(atoms_2[2 * (i - isav) - 2])

        i += 1

    return target_atoms, molobjs


def generate_brominations(mol: Mol):
    """Bromniate all aromatic carbons"""
    # __rxn1__ = AllChem.ReactionFromSmarts("[C;R;H1:1]=[C,N;R;H1:2]>>[CH:1](Br)[*H+:2]")
    # __rxn2__ = AllChem.ReactionFromSmarts("[C;R;H1:1]=[C,N;R;H0:2]>>[CH:1](Br)[*+;H0:2]")
    raise NotImplementedError
