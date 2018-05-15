import os
import subprocess
import pathlib
import BioPype.cmds.downloadsra as downloadsra
from BioPype.workspaces.dirpaths import DirPathsHelper
path_helper = DirPathsHelper()
_BIODIR = path_helper._BIODIR
_SRADIR = path_helper._SRADIR


# TODO: TURN THIS INTO A CLASS! Clearly would benefit from having class objects pass attributes to methods for each qc operation rather than having a milllion functions with similarly-named variables (short on time otherwise it would be class already)
def get_fastqc_reports(input_dir='', outdir='', select='', threads=4):
    """Use FastQC to perform quality control checks on fastq files.

    :param input_dir: str; path to directory containing target fastq files.
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
    if input_dir:
        fastq_dir = pathlib.Path(input_dir)
    else:
        fastq_dir = pathlib.Path(_SRADIR).joinpath('grouped_fastq')

    if outdir:
        qc_reports_dir = outdir
    else:
        qc_reports_dir = pathlib.Path(_SRADIR).joinpath('qcreports')
    downloadsra.DownloadHelper.check_dir_exists(str(qc_reports_dir))

    if select:
        paths = []
        for selected in select:
            # Search the subdirectories in fastq_dir for files
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
        downloadsra.DownloadHelper.check_dir_exists(str(qc_reports_subdir))

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
    downloadsra.DownloadHelper.check_dir_exists(str(qc_reports_subdir))

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


def summarize_reports(summary_filename, input_dir=''):
    """Summarizes the fastqc reports contained in the target folders.

    :param summary_filename: str; the name of the multiqc file that will be
        created.
    :param input_dir: str; the name of the directory that contains the
        target FastQC reports.
    :return: None
    """
    # Define the target directory that will be searched for fastqc reports.
    if input_dir:
        p = os.path.join(_SRADIR, input_dir)
        downloadsra.DownloadHelper.check_dir_exists(p)
        fastqc_dir = pathlib.Path(p)
    else:
        p = os.path.join(_SRADIR, 'qcreports')
        downloadsra.DownloadHelper.check_dir_exists(p)
        fastqc_dir = pathlib.Path(p)

    # Define and/or create the destination directory.
    multiqc_dir = os.path.join(_SRADIR, 'multiqc_summaries')
    downloadsra.DownloadHelper.check_dir_exists(str(multiqc_dir))

    fastqc_folder_paths = [folder_path for folder_path in fastqc_dir.iterdir()]
    for folder in fastqc_folder_paths:
        print(folder)
        # Make the multiqc reports feature the names of the folders they're
        # summarizing.
        new_report_name = str(summary_filename + "_" + folder.name)
        _run_multiqc(str(folder), multiqc_dir, new_report_name)


def _run_multiqc(fastqc_files, outdir, multiqc_filename):
    # NEEDED TO COPY THE MULTIQC EXECUTABLE FILE FROM home/anaconda/envs/biopype/bin/multiqc TO home/anaconda/bin BECAUSE FOR WHATEVER REASON THE VIRTUAL ENVIRONMENT BIN WASN'T BEING CHECKED, EVEN THOUGH I ADDED IT TO $PATH AND $PYTHONPATH, AND KEPT THROWING AND ERROR THAT SAID 'multiqc is not a command'
    subprocess.run(['multiqc', fastqc_files, '--outdir', outdir, '-f', '--filename', multiqc_filename])


def run_trim_galore(target_dir, out_dir=''):
    # TODO: adapt this so that the FastQC reports get sent to somewhere different from the trimmed fastq files?

    # TODO: determine if the downloadsra.check_dir_exists() checks are necessary at this step. might not need to call until the later 'for-loop'
    if out_dir:
        outdir = pathlib.Path(out_dir)
        downloadsra.DownloadHelper.check_dir_exists(str(outdir))
    else:
        outdir = pathlib.Path(_SRADIR).joinpath('grouped_fastq.trimmed')
        downloadsra.DownloadHelper.check_dir_exists(str(outdir))

    home = pathlib.Path(target_dir)
    for fpath in home.iterdir():
        sfile = str(fpath)
        fin_outdir = outdir.joinpath(home.name)
        downloadsra.DownloadHelper.check_dir_exists(str(fin_outdir))
        subprocess.run(
            ['trim_galore', '--output_dir', str(fin_outdir), '--fastqc', '--fastqc_args', '"--threads 4"', sfile])
    return None
