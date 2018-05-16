import os
import shutil
import subprocess
import pathlib
import csv

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
        self.fin_download_dict = {}  # THESE ARE UNFILTERED

    def download_paired_end_data(self, downloaded_sra_folder, converted_fastq_subdir):
        cfs = pathlib.Path(_SRADIR).joinpath('grouped_fastq', converted_fastq_subdir)
        self.check_dir_exists(cfs)
        # Randomly select SRA ID(s) from the pool.
        # sample = runtable.RunTable.random_sample_subset(self.acc_nums, 1)[0]
        random_samples = runtable.RunTable.random_sample_subset(self.acc_nums, 5)
        print(random_samples)
        for sample in random_samples:
            print(sample)
            # Check which files are single-end data.
            with open(os.path.join(_SRADIR, 'files', 'single_end_sra_files.txt'), 'r') as f:
                lines = [line.strip() for line in f]

            # Check if the sample has already been downloaded
            if (sample in self.download_tracker) or (sample in lines):
                # TODO: determine whether or not I should be using "return" here.
                # return self.download_paired_end_data(downloaded_sra_folder, cfs)
                continue

            # Check if the SRA file contains paired-end data.
            elif not self.isPairedSRA(sample):
                print('{} does not contain paired-end data. Restarting search.'.format(str(sample)))
                self.download_tracker.append(sample)
                self.write_to_single_end_file(sample)
                # TODO: This (and other) recursive call to download_paired_end_data seems to be working fine with 'cfs' as an argument, but I feel like that should be causing problems... if there's a strange error, try changing 'cfs' to 'converted_fastq_subdir'
                # self.download_paired_end_data(downloaded_sra_folder, cfs)
                continue

            elif self.isPairedSRA(sample):
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

                # Check if 2 fastq files were made (if the sra data were
                # paired-end, the fastq-conversion should result in 2 fastq files
                # b/c it uses --split-files as an argument).
                expected_fastq_files = self.n_files_in_fastq_dir + 2
                total_fastq_files = [file for file in destination_fastq_dirpath.iterdir()]
                if len(total_fastq_files) == expected_fastq_files:
                    self.n_files_in_fastq_dir = expected_fastq_files
                    self.success_download_counter += 1

                    # Add fastq absolute file paths and read-direction info to
                    # self.fin_download_dict.
                    s = str(sample) + '_pass_[0-9].fastq'
                    match = cfs.rglob(s)
                    # matches = []
                    for x in match:
                        # matches.append(x)
                        if 'pass_1' in str(x):
                            abs_path_1 = str(x)
                        if 'pass_2' in str(x):
                            abs_path_2 = str(x)

                    # p = pathlib.Path(abs_path_1)
                    # desc_sample_name = str(p.parent.name) + '_' + sample

                    # self.fin_download_dict[sample] = ((abs_path_1, 'forward'), (abs_path_2, 'reverse'), desc_sample_name)
                    self.fin_download_dict[sample] = (
                    (abs_path_1, 'forward'), (abs_path_2, 'reverse'))
                    print('DOWNLOAD_DICT_ENTRY: ' + str(self.fin_download_dict[sample]))

                    # Check if the counter has reached the desired number of samples.
                    if self.success_download_counter >= self.threshold:
                        return None
                    else:
                        # return self.download_paired_end_data(downloaded_sra_folder, cfs)
                        continue

                # If 2 fastq files were not made...
                else:
                    # ...get the absolute paths to the fastq files that were just created...
                    s = str(sample) + '_pass.fastq'
                    match = cfs.rglob(s)
                    # ... then delete the files.
                    for file in match:
                        # TODO: add informative print statements or write to stdout about what's happening in this loop/else block.
                        abs_path = file.resolve()
                        abs_path.unlink()

                    # Add the SRA ID to a file that tracks single-end datafiles
                    # to help speed up future downloads.
                    self.write_to_single_end_file(sample)
                    # return self.download_paired_end_data(downloaded_sra_folder, cfs)
                    continue
            else:
                raise TypeError('something went wrong... sra data were neither paired-end nor not-paired-end... Neither if-condition was satisfied.')

        if self.n_files_in_fastq_dir < self.threshold:
            self.download_paired_end_data(downloaded_sra_folder, converted_fastq_subdir)



        self.fin_download_dict[sample] = ((abs_path_1, 'forward'), (abs_path_2, 'reverse'))

    def method_create_manifest_file(self, filename):
        man_path = os.path.join(_BIODIR, filename)
        with open(man_path, 'w', newline='') as csvfile:
            manifestwriter = csv.writer(csvfile, delimiter=',')
            header = ('sample-id', 'absolute-filepath', 'direction')
            manifestwriter.writerow(header)

            for key, value in self.fin_download_dict.items():
                for tup in value:
                    manifestwriter.writerow([key, tup[0], tup[1]])
        return man_path

    @staticmethod
    def create_manifest_file(filename, *args):
        """Create fastq manifest file by using several dictionaries.

        # TODO: make this function description more detailed.
        :param filename:
        :param args: dictionaries that contain sample IDs, absolute filepaths,
        and read direction for sequence-data-files.
        :return: path to the created manifest file.
        """
        man_path = os.path.join(_BIODIR, filename)
        with open(man_path, 'w', newline='') as csvfile:
            manifestwriter = csv.writer(csvfile, delimiter=',')
            header = ('sample-id', 'absolute-filepath', 'direction')
            manifestwriter.writerow(header)

            for dictionary in args:
                print('ITEMS: ' + str(dictionary.items()))
                for key, value in dictionary.items():
                    print('KEY:' + str(key))
                    print('VALUE: ' + str(value))
                    print('LEN: ' + str(len(value)))
                    tup1 = value[0]
                    tup2 = value[1]

                    abs_path_1 = tup1[0]
                    print('ABS_PATH1 = ' + str(abs_path_1))
                    p = pathlib.Path(abs_path_1)
                    desc_sample_name = str(p.parent.name) + '_' + key

                    # desc_sample_name = value[2]
                    manifestwriter.writerow([desc_sample_name, tup1[0], tup1[1]])
                    manifestwriter.writerow([desc_sample_name, tup2[0], tup2[1]])

                    # for tup1, tup2, string in value:
                    #     print(desc_sample_name, tup[0], tup[1])
                    #     manifestwriter.writerow([desc_sample_name, tup1[0], tup1[1]])
                    #     manifestwriter.writerow([desc_sample_name, tup2[0], tup2[1]])
        return man_path

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
        print(the_samples)
        metadata = pdTableOfSamples.df.loc[pdTableOfSamples.df['Run'].isin(the_samples)]
        print(metadata)
        idx = 0
        new_col = [run for run in metadata['Run']]  # can be a list, a Series, an array or a scalar
        print(new_col)

        # Insert a column to the table because the metadata table needs an 'id' column,
        # and without this the ids end up in a weird column without a header, which
        # causes errors.
        metadata.insert(loc=idx, column='id', value=new_col)
        metadata.to_csv(metadata_filepath, sep='\t', index=False)
        return metadata_filepath

    def download_sra(self, acc_nums, outdirectory_name):
        """Download SRA runs to the configured sra-toolkit workspace.
        :param acc_nums: iterator containing SRA accession numbers (e.g., ['SRR6664514'])
        :param outdirectory_name: str; name of the folder that the downloaded
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
            grouped_dir_path = os.path.join(_SRADIR, 'grouped_sra', outdirectory_name)
            # Create the destination folder and copy the .sra files into it.
            self.check_dir_exists(grouped_dir_path)
            self.copy_sra_to_group_dir(folder, grouped_dir_path)

            # (un-comment this if the recursive copy_sra_to_group_dir function fails)
            # shutil.copytree(folder, grouped_dir_path)

        else:
            # p is the path to the .sra file that was just downloaded.
            p = os.path.join(folder, (str(acc_nums) + '.sra'))
            d = os.path.join(_SRADIR, 'grouped_sra', outdirectory_name)
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

    @staticmethod
    def check_dir_exists(directory):
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
                        item_path = item.resolve()
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

    def isPairedSRA(self, filename):
        # From https://www.biostars.org/p/139422/
        # filename = os.path.abspath(filename)
        try:
            contents = subprocess.check_output(["fastq-dump", "-X", "1", "-Z", "--split-spot", filename], universal_newlines=True)
        except subprocess.CalledProcessError as e:
            raise Exception("Error running fastq-dump on", filename)

        if contents.count("\n") == 4:
            return False
        elif contents.count("\n") == 8:
            return True
        else:
            raise Exception("Unexpected output from fast-dump on ", filename)

    def write_to_single_end_file(self, appended_text):
        """Write SRA ID's/accession numbers to a file for tracking single-end data.

        :param appended_text: Text to append to file
        :return: None
        """
        filepath = os.path.join(_SRADIR, 'files', 'single_end_sra_files.txt')
        with open(filepath, "a") as f:
            f.write(appended_text + '\n')
