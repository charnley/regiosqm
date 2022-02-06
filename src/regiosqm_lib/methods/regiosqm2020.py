# Notes
# regiosqm20 opt
# cmd = f'{XTBHOME}/bin/xtb --gfn{method} {start_structure} --opt --gbsa {solvent} --chrg {chrg} --uhf 0'
# - generate conformers (jan defaults)
# - optimize with GFN-ff / Methanol
# - Cluster and find unique conformers
# - optimize with GFN-1 / Methanol
# - return lowest energy
from typing import List

from rdkit.Chem import Mol

import regiosqm_lib


def get_best_conformer_and_energy(mol: Mol, conformer_options={}):
    """Generate conformers, minimize and keep the best conformer."""

    mol = regiosqm_lib.conformers.generate_conformers(mol, **conformer_options)
    # TODO Optimize conformers

    mol = regiosqm_lib.conformers.find_unique_conformers(mol)
    # TODO Optimize conformers 2

    # TODO Find conformer with lowest energy
    energy = 0.0
    return mol, energy


def find_relevant_tautomers(mol: Mol, energy_cutoff=15.0) -> List[Mol]:
    """
    Find energy relevant tautoemrs for a molobj
    Generate tautomers and min(3 + 3*rot_bond, max_conf) conformations followed by xTB calculations.
    """

    tautomers = regiosqm_lib.tautomers.find_tautomers(mol)

    for tautomer in tautomers:
        conformer, energy = get_best_conformer_and_energy(tautomer)

    return tautomers
