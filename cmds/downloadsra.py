import os
import subprocess
import pathlib
from workspaces.dirpaths import BIODIR, SRADIR


def download_sra(acc_nums):
    """Download SRA runs to the configured sra-toolkit workspace for .sra files.
    :param acc_nums: iterator containing SRA accession numbers (e.g., 'SRR6664514')
    :return: None
    """
    for num in acc_nums:
        subprocess.run(['prefetch', num])
    return None


def check_dir_exists(directory):
    """Check if directory exists. If not, create it."""
    if not os.path.isdir(directory):
        print(str(directory) + "does not exist... Creating " + str(directory))
        os.mkdir(str(directory))


def convert_sra_to_fastq(
        fastq_dir=(pathlib.Path(BIODIR).joinpath('data', 'fastq')),
        select_files='',
        threads=1,
        delete=False):
    """Convert target .sra files to .fastq and store in fastq_dir.
    :param fastq_dir: string; path to directory where created .fastq files will
    be stored.
    :param select_files:
        Default value: Empty string (''); Target all .sra files from sra_dir.
        Input values: An iterable containing names of files in sra_dir as
        strings. Usually, these names are SRA accession numbers + '.sra'
        (e.g., SRR664513.sra).
    :param threads: str or int; the number of threads to create when performing
    parallel-fastq-dump.
    :param delete: #TODO: figure out the clean_sra_download function before implementing this.
    :return: None
    """
    # sra_dir = path to directory where target .sra files exist.
    sra_dir = pathlib.Path(SRADIR).joinpath('sra')

    check_dir_exists(fastq_dir)

    if select_files:
        # For each file in the sra Workspace Location...
        for f in sra_dir.iterdir():
            if f.name in select_files:
                sra_file_path = sra_dir.joinpath(f.name)
                subprocess.run(['parallel-fastq-dump', '--sra-id', sra_file_path,
                                '--threads', threads, '--outdir', fastq_dir])

    else:
        for f in sra_dir.iterdir():
            if f.suffix == '.sra':
                subprocess.run(['parallel-fastq-dump', '--sra-id', f,
                                '--threads', str(threads), '--outdir', fastq_dir])

    return None


def clean_sra_downloads(delete=False, storage_dir=os.path.join(BIODIR, 'data', 'sra.storage')):
    """Remove the .sra files from the default .sra download directory.
    :param delete: When False, moves all of the files from the .sra download location to the .sra storage folder. When
    True, simply deletes the files.
    :return: None
    """
    # check_dir_exists(storage_dir)
    # TODO: Figure out what direction to go with this function.
    pass
