""" This module contains utility functions used by other modules of DocKin."""


protein_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def print_progress(information_text, counter, max_counter):
    """
    Print progress of a process.

    Parameters
    ----------
    information_text: str
        Information about the process.

    counter: int
        Number describing the current process progress.

    max_counter: int
        Number describing the maximal value of counter, i.e. end of process.

    """
    progress = round(((counter + 1) / max_counter) * 100, 1)
    print(f'{information_text}: Progress {progress} %', end='\r')
    return
