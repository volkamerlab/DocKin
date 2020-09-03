# Standard libraries
import json
import pathlib

# External libraries
import pandas as pd

# DocKin library
from dockin.assess import rmsds_from_sdfs


CWD = pathlib.Path(__file__).absolute().parent
kinase_df = pd.read_csv(CWD / 'data/re_docking_data.csv')
result_dict = dict()

for index, row in kinase_df.iterrows():
    pdb_id = row['pdb']
    ligand_id = row['ligand']
    ligand = False
    if (CWD / f'data/{pdb_id}_{ligand_id}.sdf').is_file():
        ligand = True
    receptor = False
    if (CWD / f'data/{pdb_id}_hybrid_receptor.oeb').is_file():
        receptor = True
    pose = False
    rmsd = None
    if (CWD / f'data/{pdb_id}_hybrid_re-docking.sdf').is_file():
        pose = True
        rmsds = rmsds_from_sdfs(str(CWD / f'data/{pdb_id}_{ligand_id}.sdf'),
                               str(CWD / f'data/{pdb_id}_hybrid_re-docking.sdf'))
        rmsd = rmsds[0]
    result_dict[pdb_id] = [ligand, receptor, pose, rmsd]

with open(CWD / 'data/oedocking_results.txt', 'w') as outfile:
    json.dump(result_dict, outfile)
