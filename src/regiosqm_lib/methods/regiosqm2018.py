import logging
from typing import Optional

import numpy as np
import pandas as pd
import ppqm
from ppqm import chembridge
from rdkit.Chem import Mol

import regiosqm_lib

_logger = logging.getLogger(__name__)


mopac_calculation = {
    "pm3": None,
    "eps": 4.8,
    "cycles": 200,
}


def predict_regioselective_dataframe(mol: Mol, mopac_options={}) -> pd.DataFrame:
    """Expects neutralized molecule"""

    rows = []

    # Protonate tautomers
    atoms, protomers = regiosqm_lib.generate_protonations(mol)
    for atom, protomer in zip(atoms, protomers):

        mol_ = chembridge.copy_molobj(mol)
        mol_.__sssAtoms = [atom]

        row = {"target": mol_, "molobj": protomer, "atom": atom, "tautomer_energy": 0.0}

        rows.append(row)

    pdf_protomers = pd.DataFrame(rows)
    pdf_protomers["protonation_energy"] = pdf_protomers["molobj"].progress_apply(get_best_energy)

    return pdf_protomers


def predict_regioselective_sites(
    mol: Mol, pdf_protomers: pd.DataFrame, energy_cut1=1.0, energy_cut2=3.0, mopac_options={}
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
    mol.SetProp("regiosqm2018_cut1", str(atoms_cut1))
    mol.SetProp("regiosqm2018_cut2", str(atoms_cut2))

    return mol


def get_best_energy(mol: Mol, conformer_options={}):

    mol = regiosqm_lib.generate_conformers(mol, random_seed=4, **conformer_options)

    assert mol.GetNumConformers(), "No conformers"

    mopac = ppqm.MopacCalculator()

    results = mopac.calculate(mol, mopac_calculation)

    energies = [res["h"] for res in results if "h" in res]

    energies = np.asarray(energies)
    energy = np.min(energies)

    return energy
