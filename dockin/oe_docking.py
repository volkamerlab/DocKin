"""This module contains functions to dock small molecules into proteins using the OpenEye toolkit."""

# Standard libraries
import pathlib
import tempfile

# External libraries
from openeye import oechem, oespruce, oedocking, oequacpac, oeomega
import requests

# DocKin library
from dockin.utils import protein_resnames


def get_structure_from_pdb(pdb_id):
    """
    Retrieve a structure from PDB and convert it into an OEChem molecule.

    Parameters
    ----------
    pdb_id: str
        PDB ID.

    Returns
    -------
    mol: oechem.OEGraphMol
        An oechem.OEGraphMol object holding the structure retrieved from PDB.
    """

    # get structure
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)

    # store structure in temporary file
    with tempfile.NamedTemporaryFile(suffix='.pdb') as temp_file:
        temp_file.write(response.content)

        # read structure from temporary file
        with oechem.oemolistream() as ifs:
            ifs.SetFlavor(oechem.OEFormat_PDB,
                          oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)

            # Throw error if temporary file is not available
            if not ifs.open(temp_file.name):
                oechem.OEThrow.Fatal(f'Unable to open {temp_file.name} for reading.')

            # Throw error if content of temporary file is not readable
            mol = oechem.OEGraphMol()
            if not oechem.OEReadMolecule(ifs, mol):
                oechem.OEThrow.Fatal(f'Unable to read molecule from {temp_file.name}.')

    return mol


def select_chain(mol, chain_id):
    """
    Select a chain from an OEChem molecule.

    Parameters
    ----------
    mol: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a molecular structure.

    chain_id: str
        Chain identifier.

    Returns
    -------
    selection: oechem.OEGraphMol
        An oechem.OEGraphMol object holding selected components.
    """

    # do not change input mol
    selection = mol.CreateCopy()

    # delete other chains
    for atom in selection.GetAtoms():
        residue = oechem.OEAtomGetResidue(atom)
        if residue.GetChainID() != chain_id:
            selection.DeleteAtom(atom)

    return selection


def select_altloc(mol, altloc_id):
    """
    Select an alternate location from an OEChem molecule.

    Parameters
    ----------
    mol: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a molecular structure.

    altloc_id: str
        Alternate location identifier.

    Returns
    -------
    selection: oechem.OEGraphMol
        An oechem.OEGraphMol object holding selected components.
    """

    # do not change input mol
    selection = mol.CreateCopy()

    allowed_altloc_ids = [' ', altloc_id]

    # delete other alternate location
    for atom in selection.GetAtoms():
        residue = oechem.OEAtomGetResidue(atom)
        if oechem.OEResidue.GetAlternateLocation(residue) not in allowed_altloc_ids:
            selection.DeleteAtom(atom)

    return selection


def select_ligand(mol, ligand_id):
    """
    Select all atoms from an OEChem molecule that are protein, water or are part of the specified ligand.

    Parameters
    ----------
    mol: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a molecular structure.

    ligand_id: str
        Ligand identifier.

    Returns
    -------
    selection: oechem.OEGraphMol
        An oechem.OEGraphMol object holding selected components.
    """

    # do not change input mol
    selection = mol.CreateCopy()

    allowed_resnames = protein_resnames + ['HOH', ligand_id]

    # delete other residues
    for atom in selection.GetAtoms():
        residue = oechem.OEAtomGetResidue(atom)
        if residue.GetName() not in allowed_resnames:
            selection.DeleteAtom(atom)

    return selection


def prepare_complex(protein_ligand_complex, protein_save_path=None, ligand_save_path=None):
    """
    Prepare a OEChem molecule holding a protein ligand complex for docking.

    Parameters
    ----------
    protein_ligand_complex: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a structure with protein and ligand.

    protein_save_path: str
        File path for saving prepared protein. If protein_save_path is not provided, protein will not be saved.

    ligand_save_path: str
        File path for saving prepared ligand. If ligand_save_path is not provided, ligand will not be saved.

    Returns
    -------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.

    ligand: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared ligand structure.
    """

    # pick alternate locations if defined in PDB file
    fact = oechem.OEAltLocationFactory(protein_ligand_complex)
    mol = oechem.OEGraphMol()
    fact.MakePrimaryAltMol(mol)

    # protonate complex
    print('Re-optimizing hydrogen positions...')
    opts = oechem.OEPlaceHydrogensOptions()
    details = oechem.OEPlaceHydrogensDetails()
    oechem.OEPlaceHydrogens(mol, details, opts)

    # create design units
    print('Identifying design units...')
    design_units = list(oespruce.OEMakeDesignUnits(mol))
    if len(design_units) == 1:
        design_unit = design_units[0]
    elif len(design_units) > 1:
        print('More than one design unit found---using first one')
        design_unit = design_units[0]
    else:
        raise Exception('No design units found')

    # get protein
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)

    # get ligand
    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)

    # save protein
    if protein_save_path is not None:
        with oechem.oemolostream(str(pathlib.Path(protein_save_path).absolute())) as ofs:
            oechem.OEWriteMolecule(ofs, protein)

    # save ligand
    if ligand_save_path is not None:
        with oechem.oemolostream(str(pathlib.Path(ligand_save_path).absolute())) as ofs:
            oechem.OEWriteMolecule(ofs, ligand)

    return protein, ligand


def prepare_protein(protein, protein_save_path=None):
    """
    Prepare a OEChem molecule holding a protein structure for docking.

    Parameters
    ----------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a structure with protein.

    protein_save_path: str
        File path for saving prepared protein. If protein_save_path is not provided, protein will not be saved.

    Returns
    -------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.
    """

    # pick alternate locations if defined in PDB file
    fact = oechem.OEAltLocationFactory(protein)
    mol = oechem.OEGraphMol()
    fact.MakePrimaryAltMol(mol)

    # protonate complex
    print('Re-optimizing hydrogen positions...')
    opts = oechem.OEPlaceHydrogensOptions()
    describe = oechem.OEPlaceHydrogensDetails()
    oechem.OEPlaceHydrogens(mol, describe, opts)

    # create bio design units
    print('Identifying bio design units...')
    bio_design_units = list(oespruce.OEMakeBioDesignUnits(mol))
    if len(bio_design_units) == 1:
        bio_design_unit = bio_design_units[0]
    elif len(bio_design_units) > 1:
        print('More than one design unit found---using first one')
        bio_design_unit = bio_design_units[0]
    else:
        raise Exception('No design units found')

    # get protein
    protein = oechem.OEGraphMol()
    bio_design_unit.GetProtein(protein)

    # save protein
    if protein_save_path is not None:
        with oechem.oemolostream(str(pathlib.Path(protein_save_path).absolute())) as ofs:
            oechem.OEWriteMolecule(ofs, protein)

    return protein


def create_complex_receptor(protein, ligand, receptor_save_path=None):
    """
    Prepare a OEChem molecule holding a protein structure for docking.

    Parameters
    ----------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.

    ligand: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared ligand structure.

    receptor_save_path: str
        File path for saving created receptor. If receptor_save_path is not provided, receptor will not be saved.

    Returns
    -------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a receptor with protein and ligand.
    """

    # create receptor
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)

    # save receptor
    if receptor_save_path is not None:
        oedocking.OEWriteReceptorFile(receptor, str(pathlib.Path(receptor_save_path).absolute()))

    return receptor


def write_mols(molecules, file_path):
    """
    Save molecules to file.

    Parameters
    ----------
    molecules: list of oechem.OEGraphMol
        A list of oechem.OEGraphMol objects holding molecules for writing.

    file_path: str
        File path for saving molecules.
    """
    with oechem.oemolostream() as ofs:
        if ofs.open(file_path):
            for molecules in molecules:
                oechem.OEWriteMolecule(ofs, molecules)
    return


def hybrid_docking(receptor, molecules, num_poses=1, docking_poses_save_path=None):
    """
    Dock molecules into a prepared receptor holding protein and ligand structure.

    Parameters
    ----------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding the prepared receptor.

    molecules: list of oechem.OEGraphMol
        A list of oechem.OEGraphMol objects holding prepared molecules for docking.

    num_poses: int
        Number of docking poses to generate per molecule.

    docking_poses_save_path: str
        File path for saving docking poses. If docking_poses_save_path is not provided, docking poses will not be saved.

    Returns
    -------
    docked_molecules: list of oechem.OEGraphMol
        A list of oechem.OEGraphMol objects holding the docked molecules.
    """

    def score(molecule, field='Hybrid2'):
        """Return the docking score."""
        value = oechem.OEGetSDData(molecule, field)
        return float(value)

    # initialize receptor
    dock_method = oedocking.OEDockMethod_Hybrid2
    dock_resolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dock_method, dock_resolution)
    dock.Initialize(receptor)

    docked_molecules = list()

    # dock molecules
    for molecule in molecules:
        # enumerate tautomers
        tautomer_options = oequacpac.OETautomerOptions()
        tautomer_options.SetMaxTautomersGenerated(4096)
        tautomer_options.SetMaxTautomersToReturn(16)
        tautomer_options.SetCarbonHybridization(True)
        tautomer_options.SetMaxZoneSize(50)
        tautomer_options.SetApplyWarts(True)
        pKa_norm = True
        tautomers = [oechem.OEMol(tautomer) for tautomer in
                     oequacpac.OEGetReasonableTautomers(molecule, tautomer_options, pKa_norm)]

        # set up omega
        omega_options = oeomega.OEOmegaOptions()
        omega_options.SetMaxSearchTime(60.0)  # time out
        omega = oeomega.OEOmega(omega_options)
        omega.SetStrictStereo(False)  # enumerate stereochemistry if uncertain

        docked_tautomers = list()
        # dock tautomers
        for mol in tautomers:
            docked_mol = oechem.OEMol()
            # expand conformers
            omega.Build(mol)

            # dock molecule
            return_code = dock.DockMultiConformerMolecule(docked_mol, mol, num_poses)
            if return_code != oedocking.OEDockingReturnCode_Success:
                print(f'Docking failed for molecule with title {mol.GetTitle()} with error code '
                      f'{oedocking.OEDockingReturnCodeGetName(return_code)}.')
                continue

            # store docking data
            sd_tag = oedocking.OEDockMethodGetName(dock_method)
            oedocking.OESetSDScore(docked_mol, dock, sd_tag)

            # expand conformations
            for conformation in docked_mol.GetConfs():
                docked_tautomers.append(conformation.CreateCopy())

        # sort all conformations of all tautomers by score
        docked_tautomers.sort(key=score)

        # keep number of conformations as specified by num_poses
        docked_molecules += docked_tautomers[:num_poses]

    if len(docked_molecules) == 0:
        return None

    # save docking poses
    if docking_poses_save_path is not None:
        write_mols(docked_molecules, str(pathlib.Path(docking_poses_save_path).absolute()))

    return docked_molecules