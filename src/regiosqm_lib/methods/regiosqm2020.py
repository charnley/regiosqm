# Notes
# regiosqm20 opt
# cmd = f'{XTBHOME}/bin/xtb --gfn{method} {start_structure} --opt --gbsa {solvent} --chrg {chrg} --uhf 0'
# - generate conformers (jan defaults)
# - optimize with GFN-ff / Methanol
# - Cluster and find unique conformers
# - optimize with GFN-1 / Methanol
# - return lowest energy


import logging
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import ppqm
from ppqm import chembridge, units
from rdkit import Chem
from rdkit.Chem import Mol

import regiosqm_lib

_logger = logging.getLogger(__name__)

calculation_opt1 = {
    "gfn": "ff",
    "gbsa": "Methanol",
    "opt": None,
    "uhf": 0,
}

calculation_opt2 = {
    "gfn": 1,
    "gbsa": "Methanol",
    "opt": None,
    "uhf": 0,
}


def predict_regioselective_dataframe(mol: Mol, xtb_options={}) -> pd.DataFrame:
    """Expects neutralized molecule"""

    # Generate tautomers
    # Re-order atoms based on base mol
    tautomers, energies = generate_relevant_tautomers(mol)
    tautomers = reorder_atoms(mol, tautomers)

    rows = []

    # Protonate tautomers
    for mol, energy in zip(tautomers, energies):
        atoms, protomers = regiosqm_lib.generate_protonations(mol)
        for atom, protomer in zip(atoms, protomers):

            mol_ = chembridge.copy_molobj(mol)
            mol_.__sssAtoms = [atom]

            row = {"target": mol_, "molobj": protomer, "atom": atom, "tautomer_energy": energy}

            rows.append(row)

    pdf_protomers = pd.DataFrame(rows)
    pdf_protomers["protonation_energy"] = pdf_protomers["molobj"].progress_apply(get_best_energy)

    return pdf_protomers


def predict_regioselective_sites(
    mol: Mol, energy_cut1=1.0, energy_cut2=3.0, xtb_options={}, pdf_protomers=None
) -> Optional[Mol]:
    """
    Return molobj with annotated atoms

    Assumes salt is removed
    """

    # Neutralize molobj
    mol = chembridge.neutralize_molobj(mol)

    energies_series = pdf_protomers["protonation_energy"] - pdf_protomers["tautomer_energy"]
    energies_series -= min(energies_series)

    atoms_cut1 = pdf_protomers[energies_series < energy_cut1]["atom"].tolist()
    atoms_cut2 = pdf_protomers[energies_series < energy_cut2]["atom"].tolist()
    atoms_cut2 = [idx for idx in atoms_cut2 if idx not in atoms_cut1]

    # Return mol with embedded properties
    mol.SetProp("regiosqm2020_cut1", str(atoms_cut1))
    mol.SetProp("regiosqm2020_cut2", str(atoms_cut2))

    return mol


def optimize_and_filter(mol, calculation, xtb_options={}):
    """Return energies in kcal/mol

    only place where xtb is called currently

    """
    # TODO should set global xtb
    xtb = ppqm.XtbCalculator(**xtb_options)

    results = xtb.calculate(mol, calculation)
    converged_indicies = []
    energies = []

    for i, properties in enumerate(results):

        if properties is None:
            # Conformer uncoverged
            continue

        coord = properties["coord"]
        energy = properties["total_energy"] * units.hartree_to_kcalmol
        chembridge.molobj_set_coordinates(mol, coord, confid=i)
        energies.append(energy)
        converged_indicies.append(i)

    # Remove uncoverged conformers
    mol = chembridge.molobj_select_conformers(mol, converged_indicies)

    return mol, energies


def get_best_energy(mol: Mol, conformer_options={}, xtb_options={}):
    """Generate conformers, minimize and keep the best conformer. Return the minimum energy"""

    mol = regiosqm_lib.conformers.generate_conformers(mol, **conformer_options)
    mol, energies = optimize_and_filter(mol, calculation_opt1, xtb_options=xtb_options)

    mol = regiosqm_lib.conformers.find_unique_conformers(mol)
    mol, energies = optimize_and_filter(mol, calculation_opt2, xtb_options=xtb_options)

    energy = min(energies)

    return energy


def get_energies(molobjs: Mol) -> np.ndarray:

    # TODO Make parallel

    energies = list()
    for mol in molobjs:
        energy = get_best_energy(mol)
        energies.append(energy)

    energies = np.array(energies)

    return energies


def generate_relevant_tautomers(mol: Mol, energy_cutoff=15.0) -> Tuple[List[Mol], np.ndarray]:
    """
    Find energy relevant tautoemrs for a molobj
    Generate tautomers and min(3 + 3*rot_bond, max_conf) conformations followed by xTB calculations.
    """

    tautomers = regiosqm_lib.tautomers.generate_tautomers(mol)
    _logger.info("Calculating tautomer energies")
    energies = get_energies(tautomers)

    # Find relevant graphs
    (indices,) = np.where(energies - np.min(energies) < energy_cutoff)

    return [tautomers[i] for i in indices], energies[indices]


def generate_singlebond_mol(mol: Mol) -> Mol:
    """
    Removes all double bonds, but keeps the atoms in order.
    """

    rd_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for atom in mol.GetAtoms():
        rd_atom = Chem.rdchem.Atom(atom.GetSymbol())
        rd_atom.SetAtomMapNum(atom.GetIdx())
        rd_mol.AddAtom(rd_atom)

    for atom in mol.GetAtoms():
        idx1 = atom.GetIdx()
        for atomNeighbor in atom.GetNeighbors():
            idx2 = atomNeighbor.GetIdx()
            if idx1 < idx2:
                rd_mol.AddBond(idx1, idx2, Chem.rdchem.BondType.SINGLE)

    return rd_mol.GetMol()


def reorder_atoms(reference_mol: Mol, target_mols: List[Mol]) -> List[Mol]:
    """
    Takes a list of mol objects and reorder the atoms according to a reference
    mol
    """

    reordered_mols = []
    single_bonded_ref = generate_singlebond_mol(reference_mol)
    for mol in target_mols:
        atoms = generate_singlebond_mol(mol).GetSubstructMatch(single_bonded_ref)
        reordered_mols.append(Chem.rdmolops.RenumberAtoms(mol, atoms))

    return reordered_mols
