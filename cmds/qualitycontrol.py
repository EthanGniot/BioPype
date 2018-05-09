import os
import subprocess
import pathlib
import tempfile
import BioPype.cmds.downloadsra as downloadsra
from BioPype.workspaces.dirpaths import DirPathsHelper
path_helper = DirPathsHelper()
_BIODIR = path_helper._BIODIR
_SRADIR = path_helper._SRADIR


#TODO: TURN THIS INTO A CLASS! Clearly would benefit from having class objects pass attributes to methods for each qc operation rather than having a milllion functions with similarly-named variables (short on time otherwise it would be class already)
def get_fastqc_reports(fastq_file_dir='', outdir='', select='', threads=4):
    """Use FastQC to perform quality control checks on fastq files.

    :param report_name: str; the file name for the report.
    :param fastq_file_dir: str; path to directory containing target fastq files.
        Default=''. Default path is _SRADIR > grouped_fastq
    :param outdir: str; path to the directory where the quality control report
        will be stored.
        Default=''. Default value will use _SRADIR > qcreports as the
        outdir path.
    :param select: iterable of specific fastq files OR subdirectories to analyze
        from fastq_file_dir. Iterable must contain only files or only
        subdirectories.
        Default=''. Default value will result in one report for each
        subdirectory of fastq_file_dir.
    :param threads: int or str; the number of threads that will be used to
        generate the quality control report.
    :return: None
    """
    if fastq_file_dir:
        fastq_dir = pathlib.Path(fastq_file_dir)
    else:
        fastq_dir = pathlib.Path(_SRADIR).joinpath('grouped_fastq')

    if outdir:
        qc_reports_dir = outdir
    else:
        qc_reports_dir = pathlib.Path(_SRADIR).joinpath('qcreports')
    downloadsra.check_dir_exists(str(qc_reports_dir))

    if select:
        paths = []
        for selected in select:
            # Search the subdirectories in fastq_file_dir for files
            # with names that match the user-given names in select_files.
            # The fastq_file_dir.rglob() method searches for file names that
            # contain the pattern specified by 'select_files' in both the root
            # directory (i.e., fastq_file_dir) and its subdirectories.
            # match = fastq_file_dir.rglob(str(selected))
            match = fastq_dir.rglob(str(selected))
            for item in match:
                item_path = item.resolve()
                paths.append(str(item_path))
        if os.path.isfile(paths[0]):
            _qc_report_selected_files(paths, qc_reports_dir, str(threads))
        if os.path.isdir(paths[0]):
            _qc_report_selected_dirs(fastq_dir, paths, qc_reports_dir, str(threads))

    # Create fastqc reports for all the files in each subdirectory, and
    # group the reports in subdirectories within qc_reports_dir
    else:
        fastq_folder_paths = [folder_path for folder_path in fastq_dir.iterdir()]

        _output_fastqc_files(fastq_folder_paths, qc_reports_dir, threads)

        # # For each folder...
        # for folder_path in fastq_folder_paths:
        #     # Create the name of the subdir within the qc_reports_dir where
        #     # the reports will be stored.
        #     qc_reports_subdir = qc_reports_dir.joinpath(folder_path.name)
        #
        #     # Need to create the subdir because fastqc can't output to a
        #     # non-existent directory.
        #     downloadsra.check_dir_exists(str(qc_reports_subdir))
        #
        #     # Get a list of all the files in the folder
        #     file_list = [str(file) for file in folder_path.iterdir()]
        #
        #     # Then run fastqc on the files, sending them to the outdir of
        #     # qc_reports_dir.
        #
        #     _run_fastqc(file_list, str(qc_reports_subdir), str(threads), unpack=True)


def _run_fastqc(target_file_paths, qc_reports_subdir, threads, unpack=False):
    """Execute fastqc with some pre-defined arguments.

    :param target_file_paths: iterable containing paths to the files that will
        be analyzed with fastqc.
    :param qc_reports_subdir: the name of the subdir within the qc_reports_dir where
        the reports will be stored.
    :param threads: number of threads to run fastqc with
    :return:
    """
    if unpack:
        subprocess.run(
            ['fastqc', *target_file_paths, '--outdir', str(qc_reports_subdir),
             '--extract', '-f', 'fastq', '--threads', str(threads)])
    else:
        subprocess.run(
            ['fastqc', target_file_paths, '--outdir', str(qc_reports_subdir),
             '--extract', '-f', 'fastq', '--threads', str(threads)])


def _output_fastqc_files(fastq_folder_paths_list, qc_reports_dir, threads):
    # For each folder...
    for folder_path in fastq_folder_paths_list:
        # Create the name of the subdir within the qc_reports_dir where
        # the reports will be stored.
        qc_reports_subdir = qc_reports_dir.joinpath(folder_path.name)
        downloadsra.check_dir_exists(str(qc_reports_subdir))

        # Get a list of all the files in the folder
        file_list = [str(file) for file in folder_path.iterdir()]

        # Then run fastqc on the files, sending them to the outdir of
        # qc_reports_dir.
        _run_fastqc(file_list, qc_reports_subdir, str(threads), unpack=True)


def _qc_report_selected_files(selected_files, qc_reports_dir, threads):
    """Run fastqc on a set of files.

    :param selected_files: an iterable of file paths.
    :param qc_reports_dir: path to the dir where ALL qc reports are stored.
    :return: None
    """
    def naming_recursion(name, count):
        count += 1
        # Create subdir name.
        n = name + str(count)
        # subdir_name = qc_reports_dir.joinpath(name, str(count))
        subdir_name = qc_reports_dir.joinpath(n)
        # Check if subdir already exists in qc_reports_dir.
        if subdir_name not in qc_reports_dir.iterdir():  # base case.
            return subdir_name
        else:
            return naming_recursion(name, count)

    # Define the destination directory for the files
    qc_reports_subdir = naming_recursion('selectqc_', 0)
    downloadsra.check_dir_exists(str(qc_reports_subdir))

    # Execute run_fastqc on the unpacked iterable of file paths, using
    # qc_reports_subdir as the outdir argument.
    _run_fastqc(selected_files, str(qc_reports_subdir), threads, unpack=True)


def _qc_report_selected_dirs(fastq_dir, selected_dirs, qc_reports_dir, threads):
    fastq_dir_subdirs = [subdir_path for subdir_path in fastq_dir.iterdir()]
    fastq_folder_paths = []
    for folder_path in selected_dirs:
        if pathlib.Path(folder_path) in fastq_dir_subdirs:
            fastq_folder_paths.append(pathlib.Path(folder_path))

    _output_fastqc_files(fastq_folder_paths, qc_reports_dir, threads)

    # For each folder...
    # for folder_path in fastq_folder_paths:
    #     # Create the name of the subdir within the qc_reports_dir where
    #     # the reports will be stored.
    #     qc_reports_subdir = qc_reports_dir.joinpath(folder_path.name)
    #     downloadsra.check_dir_exists(str(qc_reports_subdir))
    #
    #     # Get a list of all the files in the folder
    #     file_list = [file for file in folder_path.iterdir()]
    #
    #     # Then run fastqc on the files, sending them to the outdir of
    #     # qc_reports_dir.
    #     _run_fastqc(file_list, qc_reports_subdir, str(threads), unpack=True)


def summarize_reports(outfile_name):
    """Summarizes the

    :param outfile_name:
    :return:
    """
    report_summary_dir = os.path.join(_BIODIR, 'data', 'repsummaries')
    downloadsra.check_dir_exists(report_summary_dir)

    subprocess.run(['multiqc', ])
    pass
    # subprocess.run('multiqc?????')
    # fastq_folder_paths = [folder_path for folder_path in fastq_dir.iterdir()]
    # for folder in fastq_folder_paths:
    # # Make the qc reports feature the names of the folders they're
    # # summarizing.
    #     new_report_name = report_name + "_" + folder.name
    #     subprocess.run(
    #         ['fastqc', str(folder), '--outdir', qc_reports_dir,
    #         '--extract', '-f', 'fastq', '--threads', str(threads),
    #         '--filename', new_report_name])

    # with tempfile.TemporaryFile(mode='r+') as tf:
    #     tf.writelines(selected_files)

def trim_sequences():
    pass

