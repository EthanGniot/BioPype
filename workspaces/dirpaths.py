# Author: Ethan Gniot
import os
import warnings


class DirPathsHelper:
    def __init__(self):
        # Define the current directory as the path to the directory that
        # contains this dirpaths.py file. Need this variable so that this
        # class can find the text files that define the paths to the BioPype
        # and SRA Toolkit workspaces.
        current_dir = os.path.dirname(os.path.realpath(__file__))

        self.bio_dirname_file = os.path.join(current_dir, '__biopypeworkdir.txt')
        self.sra_dirname_file = os.path.join(current_dir, '__sratoolworkspace.txt')

        # Read the file to find out what is currently defined as the working directory
        # for BioPype.
        with open(self.bio_dirname_file, 'r') as f:
            # Use .strip() to remove the newline character at the end of the string.
            self._BIODIR = f.readline().strip()

        # Read the file to find out what is currently defined as the Workspace Location
        # for the SRA Toolkit.
        with open(self.sra_dirname_file, 'r') as f:
            self._SRADIR = f.readline().strip()

    def greeting(self):
        """
        This function is called by the BioPype package's __init__ file.
        :return: None
        """
        if self._BIODIR:
            self.work_dir_info('biopype', self._BIODIR)
        else:
            biopype_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
            self._BIODIR = biopype_dir
            self.work_dir_info('BioPype', biopype_dir)
            self.rewrite_workspace_file('biopype', biopype_dir, self.bio_dirname_file)

        if self._SRADIR:
            self.work_dir_info('sra', self._SRADIR)
        else:
            warnings.warn("***THERE IS CURRENTLY NO PATH DEFINED FOR THE SRA TOOLKIT WORKSPACE*** "
                          "BioPype relies on functionality from NCBI's "
                          "SRA Toolkit, which requires configuration of the "
                          "toolkit's Workspace Location. For information on "
                          "how to configure the SRA Toolkit Workspace"
                          " Location, please refer to either the BioPype manual "
                          "(Chapter 8: Software and Set-up), or the SRA Toolkit website "
                          "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")
        return None

    def work_dir_info(self, biopype_or_sra, work_dir):
        if biopype_or_sra.lower() == 'biopype':
            print("\n The current working directory for BioPype is: " + work_dir)
            print("\t>To change BioPype's working directory, please run: BioPype.path_helper.setup('biopype')")
            print("\t>The current working directory is where BioPype's "
                  "functions will look for target files and directories.")

        elif biopype_or_sra.lower() == 'sra':
            print("\nThe current Workspace Location for the SRA Toolkit is: " + work_dir)
            print( "\t>For information on how to configure the SRA Toolkit Workspace "
                  "Location, please refer to either the BioPype manual (Chapter 8: "
                  "Software and Set-up), or the SRA Toolkit website "
                  "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")
        return None

    def rewrite_workspace_file(self, biopype_or_sra, new_dir, workspace_file):
        """
        This function will replace the path currently given by a given
        workspace-name file with a new path. It should only be called within
        the DirPathsHelper.setup() method.
        :param biopype_or_sra: str; either 'biopype' or 'sra'.
        :param new_dir: string; the new path that will overwrite the old path.
        :param workspace_file: string; path to a file containing a workspace
        name. Will be a path to either __sratoolworkspace.txt or __biopypeworkdir.txt
        :return: None
        """
        if os.path.isdir(new_dir):
            with open(workspace_file, 'w') as w:
                w.write(new_dir
                        + "\n The first line of this file (the above line^) is the "
                          "path to the working directory of the BioPype module.")
            self._BIODIR = new_dir
            print("->" + biopype_or_sra + " directory/workspace set to: " + new_dir)

        else:
            raise AssertionError(
                new_dir
                + " was not found to be an existing directory. Please make sure the"
                  " path to the directory is formatted correctly, and that the "
                  "directory already exists. (For some assistance with file paths, "
                  "including information about absolute paths vs relative paths, "
                  "see: https://en.wikipedia.org/wiki/Path_(computing)")
        return None

    def setup(self, biopype_or_sra):
        """
        This function lets the user define the working directory for BioPype.
        :param biopype_or_sra: string; 'biopype' will let the user define the path for
        BioPype's working directory. 'sra' will let the user define the path for
        the SRA Toolkit's Workspace Location, which is where .sra files are downloaded.
        :return: None
        """

        if biopype_or_sra == 'biopype':
            working_dir = input("Please input a path to an existing folder to define BioPype's working directory: ")
            self.rewrite_workspace_file('biopype', working_dir, self.bio_dirname_file)
            self._BIODIR = working_dir

        if biopype_or_sra == 'sra':
            workspace_location = input("Please input the path for the SRA Toolkit's configured Workspace Location: ")
            self.rewrite_workspace_file('sra', workspace_location, self.sra_dirname_file)
            self._SRADIR = workspace_location
        return None
