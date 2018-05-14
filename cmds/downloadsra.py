import os
import shutil
import subprocess
import pathlib
from BioPype.workspaces.dirpaths import DirPathsHelper
import BioPype.cmds.runtable as runtable
path_helper = DirPathsHelper()
_BIODIR = path_helper._BIODIR
_SRADIR = path_helper._SRADIR


class DownloadHelper:
    def __init__(self, accession_numbers_list, select_x_samples, threads=1):
        self.acc_nums = accession_numbers_list
        self.threshold = select_x_samples
        self.threads = threads

        # Define attributes to use when tracking download of sra files.
        self.n_files_in_fastq_dir = 0
        self.download_tracker = []
        self.success_download_counter = 0

        # Create dict for successfully-downloaded paired-end files. Will
        # use IDs from this list to create metadata file and QIIME2 manifest
        # file after downloading is finished. Keys = sra_id. Value = Tuple.
        # Tuple contains 2 more tuples. Each subtuple contains 1) the absolute
        # file path to one of the paired-end fastq files, and 2) whether the
        # fastq file contains forward or reverse reads.
        self.fin_download_dict = {}

    def download_paired_end_data(self, downloaded_sra_folder, converted_fastq_subdir):
        print('THRESHOLD')
        print(self.threshold)
        print('NUM_FASTQ_FILES: ' + str(self.n_files_in_fastq_dir))
        cfs = pathlib.Path(_SRADIR).joinpath('grouped_fastq', converted_fastq_subdir)
        self.check_dir_exists(cfs)
        # Randomly select SRA ID(s) from the pool.
        sample = runtable.RunTable.random_sample_subset(self.acc_nums, 1)[0]

        # Check if the sample has already been downloaded
        if sample in self.download_tracker:
            # TODO: determine whether or not I should be using "return" here.
            return self.download_paired_end_data(downloaded_sra_folder, cfs)

        else:
            # Add the SRA ID(s) to self.download_tracker.
            self.download_tracker.append(sample)
            # Download the sample.
            self.download_sra(sample, downloaded_sra_folder)
            # Convert the sample to FASTQ format. (Should be placed in a subfolder within _SRADIR/grouped_fastq/ )
            destination_fastq_dirpath = pathlib.Path(_SRADIR).joinpath(
                'grouped_fastq',
                cfs)
            self.convert_sra_to_fastq(
                select_files=sample,
                threads=self.threads,
                _final_grouped_fastq_subdir=str(destination_fastq_dirpath))

            # Check if 2 fastq files were made (if paired-end sequencing was
            # used, the fastq-conversion should result in 2 fastq files.)
            expected_fastq_files = self.n_files_in_fastq_dir + 2
            print('EXPECTED_FASTQ_FILES: ' + str(expected_fastq_files))
            total_fastq_files = [file for file in destination_fastq_dirpath.iterdir()]
            print('TOTAL_FASTQ_FILES: ' + str(total_fastq_files))
            print(len(total_fastq_files))
            if len(total_fastq_files) == expected_fastq_files:
                self.n_files_in_fastq_dir = expected_fastq_files
                print(self.n_files_in_fastq_dir)
                self.success_download_counter += 1
                print(self.success_download_counter)

                # Add fastq absolute file paths and read-direction info to
                # self.fin_download_dict.
                print('CHECK------------------------')
                s = str(sample) + '_pass_[0-9].fastq'
                match = cfs.rglob(s)
                # matches = []
                for x in match:
                    # matches.append(x)
                    if 'pass_1' in str(x):
                        abs_path_1 = str(x)
                    if 'pass_2' in str(x):
                        abs_path_2 = str(x)
                # abs_path_1 = str(sorted(match)[0].resolve())
                # abs_path_1 = str(sorted_matches[0].resolve())
                # abs_path_2 = str(sorted(match)[1].resolve())
                # abs_path_2 = str(sorted_matches[1].resolve())
                self.fin_download_dict[sample] = ((abs_path_1, 'forward'), (abs_path_2, 'reverse'))
                print('POINT------------------------')

                # Check if the counter has reached the desired number of samples.
                print(self.success_download_counter, self.threshold)
                if self.success_download_counter >= self.threshold:
                    return None
                else:
                    return self.download_paired_end_data(downloaded_sra_folder, cfs)
            # If 2 fastq files were not made...
            else:
                # ...get the absolute paths to the fastq files that were just created...
                print('WE GOT A BAD ONE HERE!!!!!!!!!!!!!!!!!!')
                s = str(sample) + '_pass.fastq'
                match = cfs.rglob(s)
                # ... then delete the files.
                for file in match:
                    print('THERE ARE FILES BEING MATCHED!!!!!!!!!')
                    abs_path = file.resolve()
                    abs_path.unlink()
                    print('SHOULD BE DELETED NOW')
                return self.download_paired_end_data(downloaded_sra_folder, cfs)

    def create_manifest_file(self, filename):
        with open(os.path.join(_BIODIR, filename), 'w') as f:
            header = 'sample-id,absolute-filepath,direction\n'
            f.write(header)
            for key, value in self.fin_download_dict.items():
                for tup in value:
                    e = ','.join(tup)
                    entry = ','.join(key, e)
                    f.writelines(entry)
        return None

    def create_metadata_file(self, pdTableOfSamples, metadata_filepath):
        """Create a .tsv file of sample metadata for the downloaded .sra files.

        :param pdTableOfSamples: pandas DataFrame or BioPype.cmds.runtable.RunTable
        object. Contains metadata for all the samples (aka: 'Runs') of a
        study from the SRA database.
        :param metadata_filepath: full path to where the metadata file will be
        saved.
        :return: None
        """
        the_samples = [key for key in self.fin_download_dict.keys()]
        metadata = pdTableOfSamples.df.loc[pdTableOfSamples.df['Run'].isin(the_samples)]
        idx = 0
        new_col = [run for run in metadata['Run']]  # can be a list, a Series, an array or a scalar
        # Insert a column to the table because the metadata table needs an 'id' column,
        # and without this the ids end up in a weird column without a header, which
        # causes errors.
        metadata.insert(loc=idx, column='id', value=new_col)
        metadata.to_csv(metadata_filepath, sep='\t', index=False)


    def download_sra(self, acc_nums, grouped_folder_name):
        """Download SRA runs to the configured sra-toolkit workspace.
        :param acc_nums: iterator containing SRA accession numbers (e.g., ['SRR6664514'])
        :param grouped_folder_name: str; name of the folder that the downloaded
            .sra files will be stored in within data.grouped_sra
        :return: None
        """
        # Delete all files in the data.sra directory before downloading new ones.
        folder = os.path.join(_SRADIR, 'sra')
        self.check_dir_exists(folder)
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                # Un-comment the line beneath this one to delete directories too.
                # elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)

        # Download sra files using the given accession numbers.
        if not isinstance(acc_nums, str):
            for num in acc_nums:
                subprocess.run(['prefetch', num])
        else:
            subprocess.run(['prefetch', acc_nums])

        if not isinstance(acc_nums, str):
            # Define the path to the destination folder.
            grouped_dir_path = os.path.join(_SRADIR, 'grouped_sra', grouped_folder_name)
            # Create the destination folder and copy the .sra files into it.
            self.check_dir_exists(grouped_dir_path)
            self.copy_sra_to_group_dir(folder, grouped_dir_path)

            # (un-comment this if the recursive copy_sra_to_group_dir function fails)
            # shutil.copytree(folder, grouped_dir_path)

        else:
            # p is the path to the .sra file that was just downloaded.
            p = os.path.join(folder, (str(acc_nums) + '.sra'))
            d = os.path.join(_SRADIR, 'grouped_sra', grouped_folder_name)
            self.check_dir_exists(d)
            shutil.move(p, d)
        return None

    def copy_sra_to_group_dir(self, current_folder, destination_folder):
        try:
            shutil.copytree(current_folder, destination_folder)
        except FileExistsError:
            print('WARNING: DESTINATION DIRECTORY ALREADY EXISTS. SENDING FILES TO <destination_dir(copy)>\n')
            grouped_dir_path = str(destination_folder + '(copy)')
            # TODO: Determine if I need to move the recursion result to the return statement
            self.copy_sra_to_group_dir(current_folder, grouped_dir_path)
        return None

    def check_dir_exists(self, directory):
        """Check if directory exists. If not, create it."""
        if not os.path.isdir(str(directory)):
            print(str(directory) + "does not exist... Creating " + str(directory) + "\n")
            os.makedirs(str(directory))

    def parallel_fastq_dump(self, sra_file_path, threads, fastq_dir):
        """Execute a standardized parallel-fastq-dump."""
        # TODO: Change this to just accept an arbitrary number of args then change the rest of the code to account for the lack of standard args
        subprocess.run(['parallel-fastq-dump', '--sra-id', str(sra_file_path),
                        '--threads', str(threads), '--outdir', str(fastq_dir),
                        '--split-3', '--readids', '--skip-technical', '--clip',
                        '--read-filter', 'pass', '--dumpbase'])

    def convert_sra_to_fastq(self, select_files='', threads=1, _final_grouped_fastq_subdir=''):
        """Convert target .sra files to .fastq and store in fastq_dir.
        :param select_files: Must be an iterable containing names of files or subdirectories
            in data.grouped_sra as strings. Usually, these names are SRA accession
            numbers + '.sra' (e.g., SRR664513.sra) or names of folders
            (e.g,. control_sequences).
            Default=''. Default value targets all sub-directories from
            data.grouped_sra.
        :param threads: str or int; the number of threads to create when performing
        parallel-fastq-dump.
        :return: None
        """
        # sra_dir = path to the directory where target .sra files exist.
        grouped_sra_dir = pathlib.Path(_SRADIR).joinpath('grouped_sra')
        # check_dir_exists(str(temp_fastq_dir))

        # If the user only wants to convert a few specific .sra files/groups of
        # files...
        if select_files:
            paths = []
            if not isinstance(select_files, str):
                for selected in select_files:
                    print(selected)
                    # Search the subdirectories in grouped_sra_dir for files with names
                    # that match the user-given names in select_files.
                    # os.path.join('**', str(selected)) will create a string that looks
                    # like --> '**/name_of_file_in_select_files'
                    # The grouped_sra_dir.glob() method searches for file names that
                    # contain the pattern 'selected' in 'select_files'.
                    # The '**' tells grouped_sra_dir.glob() to look within
                    # subdirectories of grouped_sra_dir as well, rather than just
                    # grouped_sra_dir.
                    match = grouped_sra_dir.glob(os.path.join('**', str(selected)))
                    for item in match:
                        # print('+++++++++++++MATCH_ITEM: ' + str(item))
                        item_path = item.resolve()
                        # print('+++++++++++++ITEM_PATH: ' + str(item_path))
                        paths.append(str(item_path))
            else:
                s = select_files + '.*'
                match = grouped_sra_dir.glob(os.path.join('**', s))
                for item in match:
                    item_path = item.resolve()
                    paths.append(str(item_path))
            for path in paths:
                # If the path is a subdirectory in the data.grouped_sra
                # directory (i.e., they specified a whole group of .sra files) ...
                if os.path.isdir(str(path)):
                    subdir_path = pathlib.Path(path)
                    outdir = os.path.join(_SRADIR, 'grouped_fastq', subdir_path.name)
                    # ... for each file in the subdirectory...
                    for f in subdir_path.iterdir():
                        # ... get the path to the file...
                        sra_file_path = subdir_path.joinpath(f.name)
                        # ... then download the file from the SRA database.
                        self.parallel_fastq_dump(sra_file_path, threads, outdir)

                # If the path is a specific .sra file...
                if os.path.isfile(str(path)):
                    if _final_grouped_fastq_subdir:
                        # This is a path
                        fastq_outdir = _final_grouped_fastq_subdir
                    else:
                        fastq_outdir = os.path.join(_SRADIR, 'grouped_fastq')
                    sra_file_path = grouped_sra_dir.joinpath(path)
                    self.parallel_fastq_dump(sra_file_path, threads, fastq_outdir)

        else:  # For batch-conversion of all groups in the target grouped_sra dir.
            for folder in grouped_sra_dir.iterdir():
                outdir = os.path.join(_SRADIR, 'grouped_fastq', folder.name)
                self.check_dir_exists(outdir)
                for sra_file in folder.iterdir():
                    # This is to ignore any hidden os files (e.g., .DS_Store)
                    if sra_file.suffix == '.sra':
                        sra_file_path = grouped_sra_dir.joinpath(sra_file)
                        self.parallel_fastq_dump(sra_file_path, threads, outdir)
        return None

#
# def download_sra(acc_nums, grouped_folder_name):
#     """Download SRA runs to the configured sra-toolkit workspace for .sra files.
#     :param acc_nums: iterator containing SRA accession numbers (e.g., ['SRR6664514'])
#     :param grouped_folder_name: str; name of the folder that the downloaded
#         .sra files will be stored in within data.grouped_sra
#     :return: None
#     """
#     # Delete all files in the data.sra directory before downloading new ones.
#     folder = os.path.join(_SRADIR, 'sra')
#     check_dir_exists(folder)
#     for the_file in os.listdir(folder):
#         file_path = os.path.join(folder, the_file)
#         try:
#             if os.path.isfile(file_path):
#                 os.unlink(file_path)
#             # Un-comment the line beneath this one to delete directories too.
#             # elif os.path.isdir(file_path): shutil.rmtree(file_path)
#         except Exception as e:
#             print(e)
#
#     # Download sra files using the given accession numbers.
#     if not isinstance(acc_nums, str):
#         for num in acc_nums:
#             subprocess.run(['prefetch', num])
#     else:
#         subprocess.run(['prefetch', acc_nums])
#
#     # Define the path to the destination folder.
#     grouped_dir_path = os.path.join(_SRADIR, 'grouped_sra', grouped_folder_name)
#     # Create the destination folder and copy the .sra files into it.
#     copy_sra_to_group_dir(folder, grouped_dir_path)
#
#     # (un-comment this if the recursive copy_sra_to_group_dir function fails)
#     # shutil.copytree(folder, grouped_dir_path)
#
#     return None
#
#
# def copy_sra_to_group_dir(sra_folder, grouped_folder_name):
#     try:
#         shutil.copytree(sra_folder, grouped_folder_name)
#     except FileExistsError:
#         print('WARNING: DESTINATION DIRECTORY ALREADY EXISTS. SENDING FILES TO <destination_dir(copy)>\n')
#         grouped_dir_path = str(grouped_folder_name + '(copy)')
#         #TODO: Determine if I need to move the recursion result to the return statement
#         copy_sra_to_group_dir(sra_folder, grouped_dir_path)
#     return None
#
#
# def check_dir_exists(directory):
#     """Check if directory exists. If not, create it."""
#     if not os.path.isdir(str(directory)):
#         print(str(directory) + "does not exist... Creating " + str(directory) + "\n")
#         os.makedirs(str(directory))
#
#
# def parallel_fastq_dump(sra_file_path, threads, fastq_dir):
#     """Execute a standardized parallel-fastq-dump."""
#     # TODO: Change this to just accept an arbitrary number of args then change the rest of the code to account for the lack of standard args
#     subprocess.run(['parallel-fastq-dump', '--sra-id', str(sra_file_path),
#                     '--threads', str(threads), '--outdir', str(fastq_dir),
#                     '--split-3', '--readids', '--skip-technical', '--clip',
#                     '--read-filter', 'pass', '--dumpbase'])
#
#
# def convert_sra_to_fastq(select_files='', threads=1, _final_grouped_fastq_subdir=''):
#     """Convert target .sra files to .fastq and store in fastq_dir.
#     :param select_files: Must be an iterable containing names of files or subdirectories
#         in data.grouped_sra as strings. Usually, these names are SRA accession
#         numbers + '.sra' (e.g., SRR664513.sra) or names of folders
#         (e.g,. control_sequences).
#         Default=''. Default value targets all sub-directories from
#         data.grouped_sra.
#     :param threads: str or int; the number of threads to create when performing
#     parallel-fastq-dump.
#     :return: None
#     """
#     # sra_dir = path to the directory where target .sra files exist.
#     grouped_sra_dir = pathlib.Path(_SRADIR).joinpath('grouped_sra')
#     # check_dir_exists(str(temp_fastq_dir))
#
#     # If the user only wants to convert a few specific .sra files/groups of
#     # files...
#     if select_files:
#         paths = []
#         for selected in select_files:
#             # Search the subdirectories in grouped_sra_dir for files with names
#             # that match the user-given names in select_files.
#             # os.path.join('**', str(selected)) will create a string that looks
#             # like --> '**/name_of_file_in_select_files'
#             # The grouped_sra_dir.glob() method searches for file names that
#             # contain the pattern 'selected' in 'select_files'.
#             # The '**' tells grouped_sra_dir.glob() to look within
#             # subdirectories of grouped_sra_dir as well, rather than just
#             # grouped_sra_dir.
#             match = grouped_sra_dir.glob(os.path.join('**', str(selected)))
#             for item in match:
#                 item_path = item.resolve()
#                 paths.append(str(item_path))
#
#         for path in paths:
#             # If the path is a subdirectory in the data.grouped_sra
#             # directory (i.e., they specified a whole group of .sra files) ...
#             if os.path.isdir(str(path)):
#                 subdir_path = pathlib.Path(path)
#                 outdir = os.path.join(_SRADIR, 'grouped_fastq', subdir_path.name)
#                 # ... for each file in the subdirectory...
#                 for f in subdir_path.iterdir():
#                     # ... get the path to the file...
#                     sra_file_path = subdir_path.joinpath(f.name)
#                     # ... then download the file from the SRA database.
#                     parallel_fastq_dump(sra_file_path, threads, outdir)
#
#             # If the path is a specific .sra file...
#             if os.path.isfile(str(path)):
#                 if _final_grouped_fastq_subdir:
#                     # This is a path
#                     fastq_outdir = _final_grouped_fastq_subdir
#                 else:
#                     fastq_outdir = os.path.join(_SRADIR, 'grouped_fastq')
#                 sra_file_path = grouped_sra_dir.joinpath(path)
#                 parallel_fastq_dump(sra_file_path, threads, fastq_outdir)
#
#     else:  # For batch-conversion of all groups in the target grouped_sra dir.
#         for folder in grouped_sra_dir.iterdir():
#             outdir = os.path.join(_SRADIR, 'grouped_fastq', folder.name)
#             check_dir_exists(outdir)
#             for sra_file in folder.iterdir():
#                 # This is to ignore any hidden os files (e.g., .DS_Store)
#                 if sra_file.suffix == '.sra':
#                     sra_file_path = grouped_sra_dir.joinpath(sra_file)
#                     parallel_fastq_dump(sra_file_path, threads, outdir)
#
#     return None

# sra = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_sra/uc_old/SRR6180539.sra'
# testpath = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq'
# parallel_fastq_dump(sra, 4, testpath)

# convert_sra_to_fastq(threads=4)