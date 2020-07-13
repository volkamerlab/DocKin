# Standard libraries
import argparse
import multiprocessing
import pathlib

# External libraries
import klifs_utils
from openeye import oechem

# DocKin library
from dockin.assess import rmsds_from_sdfs
from dockin.oe_docking import get_structure_from_pdb, select_chain, select_altloc, select_ligand, prepare_complex, \
    create_hybrid_receptor, hybrid_docking


def redocking(jobs, log_file, jobs_number):
    """
    Redock co-crystalized ligands to their protein and redirect any warning or error to the specified log file.

    Parameters
    ----------
    jobs: list of pd.Series
        List of pandas series containing the information to run the docking.

    log_file: str
        Path to file for wrinting log messages.

    jobs_number:
        Total number of jobs to run
    """

    # redirect OpenEye logging
    ofs = oechem.oeofstream()
    ofs.open(log_file)
    oechem.OEThrow.SetOutputStream(ofs)

    while jobs:
        # run docking jobs until jobs list is empty
        try:
            job = jobs.pop(0)
        except IndexError:
            pass
        # perform docking and catch any error
        try:
            pdb_id = job['pdb']
            oechem.OEThrow.Info(f'Starting DocKin for {pdb_id}')
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
            rmsd = rmsds_from_sdfs(str(CWD / f'data/{pdb_id}_{ligand_id}.sdf'),
                                   str(CWD / f'data/{pdb_id}_hybrid_re-docking.sdf'))[0]
            oechem.OEThrow.Info(f'Redocking RMSD for {pdb_id}: {rmsd}')
            print(f'Redocking: Finished {jobs_number - len(jobs)}/{jobs_number}')
        except:  # not PEP8, but want to catch any error
            oechem.OEThrow.Info('DocKin error!')
            continue
    return


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(prog='DocKin', description='Redock all ligands to kinases',
                                     formatter_class=argparse.RawTextHelpFormatter)
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
    manager = multiprocessing.Manager()
    jobs = manager.list()
    for index, row in kinase_df.iterrows():
        jobs.append(row)

    # starting jobs on separate processes
    processes = [multiprocessing.Process(target=redocking,
                                         args=(jobs, str(CWD / f'data/oe_failure_{process_counter}.log'),
                                               len(jobs))) for process_counter in range(num_processes)]
    for process in processes:
        process.start()
    for process in processes:
        process.join()

    print('Finished!')
