"""
DocKin
Docking into Kinases
"""

# Add imports here
from .assess import rmsds_from_sdfs
from .oe_docking import get_structure_from_pdb, select_chain, select_altloc, select_ligand, prepare_complex, \
    prepare_protein, create_hybrid_receptor, hybrid_docking, superpose_proteins
from .utils import print_progress, protein_resnames

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
