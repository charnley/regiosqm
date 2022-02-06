import logging
from typing import List

from ppqm import chembridge
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize

_logger = logging.getLogger(__name__)


default_xtb_options = {
    "solvent": "Methanol",
    "method": "gfn1",
}


def find_tautomers(mol: Mol, sort=True) -> List[Mol]:
    """Find tautomers for a molobj"""
    # Find canonical taut
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(mol)

    # Find all unique tauts
    tauts = enumerator.Enumerate(mol)
    tauts = list(tauts)
    tauts += [canon]
    tauts = chembridge.unique(tauts)

    # Sort according to highest taut score
    if sort:
        scores = [enumerator.ScoreTautomer(x) for x in tauts]
        tauts = [tau for _, tau in sorted(zip(scores, tauts), key=lambda pair: pair[0])]

    return tauts


# def rdkit_taut_gen(rdkit_mol, name='mol', min_conf=1, rot_conf=3, max_conf=20, chrg=0,
#                 method='1', solvent='Methanol', conf_cutoff=3, taut_cutoff=15, calc_dir=os.getcwd(), dirprefix=''):
#     """
#     Generate tautomers and min(3 + 3*rot_bond, max_conf) conformations followed by xTB calculations.
#     """
#
#     # Create tauts
#     tauts = reorderTautomers(rdkit_mol) # Canonicalize list and omit duplicates
#
#     # Run calculations on tauts
#     taut_names_list = list()
#     taut_smiles_list = list()
#     taut_energy_list = list()
#     for taut_index, taut in enumerate(tauts):
#
#         taut_name=f"{name}_t{taut_index}"
#
#         best_conf, best_energy = acc_conf_search(taut, name=taut_name, min_conf=min_conf, rot_conf=rot_conf, max_conf=max_conf, chrg=chrg,
#                                                 method=method, solvent=solvent, conf_cutoff=conf_cutoff, calc_dir=calc_dir, dirprefix=dirprefix)
#
#         taut_names_list.append(taut_name) #best_conf or taut_name
#         taut_smiles_list.append(Chem.MolToSmiles(taut, canonical=False))
#         taut_energy_list.append(best_energy)
#
#     # Find the tautomers below cutoff
#     rel_energies = np.array(taut_energy_list) - np.min(taut_energy_list) #covert to relative energies
#     below_cutoff = (rel_energies <= taut_cutoff).sum() #get number of conf below cutoff
#     conf_tuble = list(zip(taut_names_list,taut_smiles_list,taut_energy_list,rel_energies)) #make a tuble
#     conf_tuble = sorted(conf_tuble, key=itemgetter(3))[:below_cutoff] #get only the best conf below cutoff
#     sorted_taut_names_list = [item[0] for item in conf_tuble]
#     sorted_taut_smiles_list = [item[1] for item in conf_tuble]
#     sorted_taut_energy_list = [item[2] for item in conf_tuble]
#
#     return sorted_taut_names_list, sorted_taut_smiles_list, sorted_taut_energy_list
#
