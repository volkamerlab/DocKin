# Standard libraries
import argparse
import math
import multiprocessing
import pathlib

# External libraries
import klifs_utils

# DocKin library
from dockin.oe_docking import get_structure_from_pdb, select_chain, select_altloc, select_ligand, prepare_complex, \
    create_hybrid_receptor, hybrid_docking


def redocking(jobs):
    """
    Redock co-crystalized ligands to their protein and redirect any warning or error to the specified log file.

    Parameters
    ----------
    jobs: pd.Dataframe
        Pandas dataframe containing the information to run the docking.
    """
    for i, job in jobs.iterrows():
        # perform docking and catch any error
        try:
            pdb_id = job['pdb']
            structure_complex = get_structure_from_pdb(pdb_id)
            structure_complex = select_chain(structure_complex, job['chain'])
            altloc = job['alt']
            if altloc != '':
                structure_complex = select_altloc(structure_complex, altloc)
            ligand_id = job['ligand']
            structure_complex = select_ligand(structure_complex, ligand_id)
            protein, ligand = prepare_complex(structure_complex,
                                              ligand_save_path=CWD / f'data/{pdb_id}_{ligand_id}.sdf')
            hybrid_receptor = create_hybrid_receptor(protein, ligand, CWD / f'data/{pdb_id}_hybrid_receptor.oeb')
            docking_poses = hybrid_docking(hybrid_receptor, [ligand],
                                           docking_poses_save_path=CWD / f'data/{pdb_id}_hybrid_re-docking.sdf')
        except:  # not PEP8, but want to catch any error
            continue
    return


def jobs_to_chunks(jobs_df, worker_number):
    chunks = []
    for i in range(worker_number):
        chunk_jobs_number = math.ceil(len(jobs_df) / (worker_number - i))
        chunks.append(jobs_df.iloc[0:chunk_jobs_number])
        jobs_df = jobs_df.iloc[chunk_jobs_number:]
    return chunks


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(prog='DocKin', description='Redock all ligands to kinases')
    parser.add_argument('-p', dest='num_processes', help='number of parallel processes', type=int, default=1)
    num_processes = int(parser.parse_args().num_processes)

    CWD = pathlib.Path(__file__).absolute().parent

    # retrieve data from KLIFS
    kinase_ids = klifs_utils.remote.kinases.kinase_names().kinase_ID.to_list()
    kinase_df = klifs_utils.remote.structures.structures_from_kinase_ids(kinase_ids)
    print('Number of PDB entries:', len(set(kinase_df['pdb'])))
    print('Number of KLIFS entries:', len(kinase_df))

    # filter for orthosteric ligand
    kinase_df = kinase_df[kinase_df['ligand'] != 0]
    print('Number of PDB entries with orthosteric ligand:', len(set(kinase_df['pdb'])))
    print('Number of KLIFS entries with orthosteric ligand:', len(kinase_df))

    # filter for structures with a single orthosteric ligand
    kinase_df = kinase_df.groupby('pdb').filter(lambda x: len(set(x['ligand'])) == 1)
    print('Number of PDB entries with a single orthosteric ligand:', len(set(kinase_df['pdb'])))
    print('Number of KLIFS entries  with a single orthosteric ligand:', len(kinase_df))

    # sort by alt, chain and quality score to pick representative structure in next step
    kinase_df = kinase_df.sort_values(by=['alt', 'chain', 'quality_score'], ascending=[True, True, False])

    # keep entry with highest quality score (alt 'A' preferred over alt 'B', chain 'A' preferred over 'B')
    kinase_df = kinase_df.groupby('pdb').head(1)
    print('Number of unique PDB entries with a single orthosteric ligand:', len(set(kinase_df['pdb'])))
    print('Number of unique KLIFS entries with a single orthosteric ligand:', len(kinase_df))

    # save data
    kinase_df.to_csv(CWD / 'data/re_docking_data.csv')

    # crating shareable jobs list
    chunks = jobs_to_chunks(kinase_df, num_processes)

    # starting jobs on separate processes
    processes = [multiprocessing.Process(target=redocking, args=(chunk,)) for chunk in chunks]
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    print('Finished!')
