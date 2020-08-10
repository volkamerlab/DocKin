"""This module contains functions to dock small molecules into proteins using the OpenEye toolkit."""


def get_structure_from_pdb(pdb_id):
    """
    Retrieve a structure from PDB and convert it into an OpenEye molecule.

    Parameters
    ----------
    pdb_id: str
        PDB ID.

    Returns
    -------
    mol: oechem.OEGraphMol
        An oechem.OEGraphMol object holding the structure retrieved from PDB.
    """
    # Standard libraries
    import tempfile

    # External libraries
    from openeye import oechem
    import requests

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

            # Print error if temporary file is not available
            if not ifs.open(temp_file.name):
                mol = None
                print(f'Unable to open {temp_file.name} for PDB ID {pdb_id}. Returning None')

            # Print error if content of temporary file is not readable
            mol = oechem.OEGraphMol()
            if not oechem.OEReadMolecule(ifs, mol):
                mol = None
                print(f'Unable to read molecule from {temp_file.name} for PDB ID {pdb_id}. Returning None.')

    return mol


def get_electron_density_from_pdb(pdb_id):
    """
    Retrieve an electron density from PDB and convert it into an OpenEye grid.

    Parameters
    ----------
    pdb_id: str
        PDB ID.

    Returns
    -------
    electron_density: oegrid.OESkewGrid
        An oegrid.OESkewGrid object holding the electron density retrieved from PDB.
    """
    # Standard libraries
    import tempfile

    # External libraries
    from openeye import oegrid
    import requests

    # get structure
    url = f'https://edmaps.rcsb.org/coefficients/{pdb_id}.mtz'
    response = requests.get(url)

    # store structure in temporary file
    with tempfile.NamedTemporaryFile(suffix='.mtz') as temp_file:
        temp_file.write(response.content)

        # read electron_density from temporary file
        electron_density = oegrid.OESkewGrid()
        if not oegrid.OEReadMTZ(temp_file.name, electron_density, oegrid.OEMTZMapType_Fwt):
            electron_density = None
            print(f'Unable to electron density from {temp_file.name} for PDB ID {pdb_id}. Returning None.')

    return electron_density


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
    # External libraries
    from openeye import oechem

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
    # External libraries
    from openeye import oechem

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
    # External libraries
    from openeye import oechem

    # DocKin library
    from dockin.utils import protein_resnames

    # do not change input mol
    selection = mol.CreateCopy()

    allowed_resnames = protein_resnames + ['HOH', ligand_id]

    # delete other residues
    for atom in selection.GetAtoms():
        residue = oechem.OEAtomGetResidue(atom)
        if residue.GetName() not in allowed_resnames:
            selection.DeleteAtom(atom)

    return selection


def prepare_complex(protein_ligand_complex, electron_density=None, loop_db=None, protein_save_path=None,
                    ligand_save_path=None):
    """
    Prepare a OEChem molecule holding a protein ligand complex for docking.

    Parameters
    ----------
    protein_ligand_complex: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a structure with protein and ligand.

    electron_density: oegrid.OESkewGrid
        An oegrid.OESkewGrid object holding the electron density.

    loop_db: str
        File path for OpenEye Spruce loop database.

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
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oespruce

    # create design units
    structure_metadata = oespruce.OEStructureMetadata()
    design_unit_options = oespruce.OEMakeDesignUnitOptions()
    if loop_db is not None:
        design_unit_options.GetPrepOptions().GetBuildOptions().GetLoopBuilderOptions().SetLoopDBFilename(loop_db)
    if electron_density is None:
        design_units = list(oespruce.OEMakeDesignUnits(protein_ligand_complex, structure_metadata, design_unit_options))
    else:
        design_units = list(oespruce.OEMakeDesignUnits(protein_ligand_complex, electron_density, structure_metadata,
                                                       design_unit_options))
    if len(design_units) == 1:
        design_unit = design_units[0]
    elif len(design_units) > 1:
        print('More than one design unit found---using first one')
        design_unit = design_units[0]
    else:
        print('No design units found')
        return None, None

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


def prepare_protein(protein, loop_db=None, protein_save_path=None):
    """
    Prepare a OEChem molecule holding a protein structure for docking.

    Parameters
    ----------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a structure with protein.

    loop_db: str
        File path for OpenEye Spruce loop database.

    protein_save_path: str
        File path for saving prepared protein. If protein_save_path is not provided, protein will not be saved.

    Returns
    -------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.
    """
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oespruce

    # create bio design units
    structure_metadata = oespruce.OEStructureMetadata()
    design_unit_options = oespruce.OEMakeDesignUnitOptions()
    if loop_db is not None:
        design_unit_options.GetPrepOptions().GetBuildOptions().GetLoopBuilderOptions().SetLoopDBFilename(loop_db)
    bio_design_units = list(oespruce.OEMakeBioDesignUnits(protein, structure_metadata, design_unit_options))
    if len(bio_design_units) == 1:
        bio_design_unit = bio_design_units[0]
    elif len(bio_design_units) > 1:
        print('More than one design unit found---using first one')
        bio_design_unit = bio_design_units[0]
    else:
        print('No design units found')
        return None

    # get protein
    protein = oechem.OEGraphMol()
    bio_design_unit.GetProtein(protein)

    # save protein
    if protein_save_path is not None:
        with oechem.oemolostream(str(pathlib.Path(protein_save_path).absolute())) as ofs:
            oechem.OEWriteMolecule(ofs, protein)

    return protein


def create_hybrid_receptor(protein, ligand, receptor_save_path=None):
    """
    Create a receptor for hybrid docking.

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
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oedocking

    # create receptor
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)

    # save receptor
    if receptor_save_path is not None:
        oedocking.OEWriteReceptorFile(receptor, str(pathlib.Path(receptor_save_path).absolute()))

    return receptor


def create_hint_receptor(protein, hintx, hinty, hintz, receptor_save_path=None):
    """
    Create a hint receptor for docking.

    Parameters
    ----------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.

    hintx: float or int
        A number defining the hint x coordinate.

    hinty: float or int
        A number defining the hint y coordinate.

    hintz: float or int
        A number defining the hint z coordinate.

    receptor_save_path: str
        File path for saving created receptor. If receptor_save_path is not provided, receptor will not be saved.

    Returns
    -------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a receptor with defined binding site via hint coordinates.
    """
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oedocking

    # create receptor
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, hintx, hinty, hintz)

    # save receptor
    if receptor_save_path is not None:
        oedocking.OEWriteReceptorFile(receptor, str(pathlib.Path(receptor_save_path).absolute()))

    return receptor


def create_box_receptor(protein, xmax, ymax, zmax, xmin, ymin, zmin, receptor_save_path=None):
    """
    Create a box receptor for docking.

    Parameters
    ----------
    protein: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a prepared protein structure.

    xmax: float or int
        Maximal number in x direction.

    ymax: float or int
        Maximal number in y direction.

    zmax: float or int
        Maximal number in z direction.

    xmin: float or int
        Minimal number in x direction.

    ymin: float or int
        Minimal number in x direction.

    zmin: float or int
        Minimal number in z direction.

    receptor_save_path: str
        File path for saving created receptor. If receptor_save_path is not provided, receptor will not be saved.

    Returns
    -------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a receptor with defined binding site via box dimensions.
    """
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oedocking

    # create receptor
    box = oedocking.OEBox(xmax, ymax, zmax, xmin, ymin, zmin)
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, box)

    # save receptor
    if receptor_save_path is not None:
        oedocking.OEWriteReceptorFile(receptor, str(pathlib.Path(receptor_save_path).absolute()))

    return receptor


def get_mol_centroid(molecule):
    """
    Calculate the geometric center of a molecules heavy atom coordinates.

    Parameters
    ----------
    molecule: oechem.OEGraphMol
        An oechem.OEGraphMol object holding a molecule whos centroid should be calculated.

    Returns
    -------
    centroid: tuple of float
        The geometric centroid of the molecule based on heavy atom coordinates
    """

    coordinates = list(molecule.GetCoords().values())
    atoms = list(molecule.GetAtoms())

    x = [coordinate[0] for coordinate, atom in zip(coordinates, atoms) if atom.GetAtomicNum() != 1]
    y = [coordinate[1] for coordinate, atom in zip(coordinates, atoms) if atom.GetAtomicNum() != 1]
    z = [coordinate[2] for coordinate, atom in zip(coordinates, atoms) if atom.GetAtomicNum() != 1]

    centroid = (sum(x) / len(x), sum(y) / len(y), sum(z) / len(z))

    return centroid


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
    # External libraries
    from openeye import oechem

    with oechem.oemolostream(file_path) as ofs:
        for molecule in molecules:
            oechem.OEWriteMolecule(ofs, molecule)
    return


def run_docking(receptor, molecules, dock_method, num_poses=1, docking_poses_save_path=None):
    """
    Dock molecules into a prepared receptor.

    Parameters
    ----------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding the prepared receptor.

    molecules: list of oechem.OEGraphMol
        A list of oechem.OEGraphMol objects holding prepared molecules for docking.

    dock_method: int
        Constant defining the docking method.

    num_poses: int
        Number of docking poses to generate per molecule.

    docking_poses_save_path: str
        File path for saving docking poses. If docking_poses_save_path is not provided, docking poses will not be saved.

    Returns
    -------
    docked_molecules: list of oechem.OEGraphMol
        A list of oechem.OEGraphMol objects holding the docked molecules.
    """
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oedocking, oequacpac, oeomega

    # initialize receptor
    dock_resolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dock_method, dock_resolution)
    dock.Initialize(receptor)

    def score(molecule, dock=dock):
        """Return the docking score."""
        value = oechem.OEGetSDData(molecule, dock.GetName())
        return float(value)

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
            oedocking.OESetSDScore(docked_mol, dock, dock.GetName())

            # expand conformations
            for conformation in docked_mol.GetConfs():
                docked_tautomers.append(oechem.OEGraphMol(conformation))

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


def hybrid_docking(hybrid_receptor, molecules, num_poses=1, docking_poses_save_path=None):
    """
    Dock molecules into a prepared receptor holding protein and ligand structure.

    Parameters
    ----------
    hybrid_receptor: oechem.OEGraphMol
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
    # External libraries
    from openeye import oedocking

    dock_method = oedocking.OEDockMethod_Hybrid2
    docked_molecules = run_docking(hybrid_receptor, molecules, dock_method, num_poses, docking_poses_save_path)

    return docked_molecules


def chemgauss_docking(receptor, molecules, num_poses=1, docking_poses_save_path=None):
    """
    Dock molecules into a prepared receptor holding a protein structure.

    Parameters
    ----------
    receptor: oechem.OEGraphMol
        An oechem.OEGraphMol object holding the prepared hint or box receptor .

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
    # External libraries
    from openeye import oedocking

    dock_method = oedocking.OEDockMethod_Chemgauss4
    docked_molecules = run_docking(receptor, molecules, dock_method, num_poses, docking_poses_save_path)

    return docked_molecules


def superpose_proteins(reference_protein, fit_protein, superposed_protein_save_path=None):
    """
    Superpose a protein structure onto a reference protein.

    Parameters
    ----------
    reference_protein: oechem.OEGraphMol
        An oechem.OEGraphMol objects holding a protein structure which will be used as reference during superposition.

    fit_protein: oechem.OEGraphMol
        An oechem.OEGraphMol objects holding a protein structure which will be superposed onto the reference protein.

    superposed_protein_save_path: str
        File path for saving the superposed protein. If superposed_protein_save_path is not provided, the superpsoed
        protein will not be saved.

    Returns
    -------
    rmsd: float
        The RMSD between the superposed structures.

    superposed_protein: oechem.OEGraphMol
        An oechem.OEGraphMol objects holding the superposed protein structure.
    """
    # Standard libraries
    import pathlib

    # External libraries
    from openeye import oechem, oespruce

    # do not modify input
    superposed_protein = fit_protein.CreateCopy()

    # set superposition method
    options = oespruce.OESuperpositionOptions()
    options.SetSuperpositionType(oespruce.OESuperpositionType_Global)

    # perform superposition
    superposition = oespruce.OEStructuralSuperposition(reference_protein, superposed_protein, options)
    rmsd = superposition.GetRMSD()
    superposition.Transform(superposed_protein)

    # save superposed structure
    if superposed_protein_save_path is not None:
        with oechem.oemolostream(str(pathlib.Path(superposed_protein_save_path).absolute())) as ofs:
            oechem.OEWriteMolecule(ofs, superposed_protein)

    return rmsd, superposed_protein
