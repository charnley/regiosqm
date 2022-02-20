import logging

from ppqm import chembridge
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdMolDescriptors
from rdkit.ML.Cluster import Butina

_logger = logging.getLogger(__name__)


def generate_conformers(
    mol: Mol,
    minimum_conformers=1,
    maximum_conformers=20,
    per_rotation_conformers=3,
    random_seed: int = -1,
    inplace=False,
):
    """Generate conformers with rdkit"""

    if not inplace:
        mol = chembridge.copy_molobj(mol)

    mol = Chem.AddHs(mol)

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(mol)
    confs = min(minimum_conformers + per_rotation_conformers * rot_bond, maximum_conformers)

    AllChem.EmbedMultipleConfs(
        mol,
        numConfs=confs,
        clearConfs=True,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        ETversion=2,
        randomSeed=random_seed,
    )

    return mol


def find_unique_conformers(mol: Mol, threshold=0.5):
    """
    Clustering conformers with RDKit's Butina algorithm to find unique
    conformer using either heavy-atom root mean square deviation (RMSD) or
    heavy-atom torsion fingerprint deviation (TFD)
    """

    # calculate difference matrix
    diffmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)

    # Cluster conformers
    num_confs = mol.GetNumConformers()
    groups = Butina.ClusterData(diffmat, num_confs, threshold, isDistData=True, reordering=True)

    idxs = [group[0] for group in groups]

    mol = chembridge.molobj_select_conformers(mol, idxs)

    return mol
