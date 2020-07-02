"""This module contains functions to assess docking performance."""

# External libraries
from rdkit import Chem
from rdkit.Chem import AllChem


def rdkit_mols_from_sdf(file_path):
    """
    Retrieve rdkit molecules from an SDF file.

    Parameters
    ----------
    file_path: str
        Path to SDF file.

    Returns
    -------
    molecules: list of rdkit.Chem.rdchem.Mol
        List of read molecules.
    """

    molecules = list()
    supplier = Chem.SDMolSupplier(file_path)
    for molecule in supplier:
        molecules.append(molecule)

    return molecules


def rmsds_from_sdfs(file_path_1, file_path_2):
    """
    Calculate RMSDs between two sets of molecules stored in two SDF files.

    Parameters
    ----------
    file_path_1: str
        Path to first SDF file.
    file_path_2: str
        Path to second SDF file.

    Returns
    -------
    rmsds: list of float
        List of RMSD values for each pair of molecules stored in provided SDF files.
    """

    # read SDF files
    poses1 = rdkit_mols_from_sdf(file_path_1)
    poses2 = rdkit_mols_from_sdf(file_path_2)

    # calculate RMSDs
    rmsds = list()
    for pose1, pose2 in zip(poses1, poses2):
        rmsds.append(AllChem.GetBestRMS(pose1, pose2))

    return rmsds
