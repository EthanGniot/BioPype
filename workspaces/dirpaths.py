import os
import warnings


class DirPathsHelper:
    def __init__(self):
        # Define the current directory as the path to the directory that
        # contains this dirpaths.py file. Need this variable so that this
        # class can find the text files that define the paths to the BioPype
        # and SRA Toolkit workspaces.
        # current_dir = os.path.join(os.path.dirname(__file__))
        current_dir = os.path.dirname(os.path.realpath(__file__))
        print("The current_dir of DirPathsHelper object is : " + str(current_dir))

        self.bio_dirname_file = os.path.join(current_dir, 'biopypeworkdir.txt')
        self.sra_dirname_file = os.path.join(current_dir, 'sratoolworkspace.txt')

        # Read the file to find out what is currently defined as the working directory
        # for BioPype.
        with open(self.bio_dirname_file, 'r') as f:
            # Use .strip() here to remove the newline character at the end of the string.
            BIODIR = f.readline().strip()

            if BIODIR:
                self.work_dir_info('biopype', BIODIR)

            else:
                self.work_dir_info('BioPype', os.getcwd())

        # Read the file to find out what is currently defined as the Workspace Location
        # for the SRA Toolkit.
        with open(self.sra_dirname_file, 'r') as f:
            SRADIR = f.readline().strip()

            if SRADIR:
                self.work_dir_info('sra', SRADIR)

            else:
                warnings.warn("***THERE IS CURRENTLY NO PATH DEFINED FOR THE SRA TOOLKIT WORKSPACE*** "
                              "BioPype relies on functionality from NCBI's "
                              "SRA Toolkit, which requires configuration of the "
                              "toolkit's Workspace Location after it is installed. For "
                              "information on how to configure the SRA Toolkit Workspace"
                              " Location, please refer to either the BioPype manual "
                              "(Chapter 8: Software and Set-up), or the SRA Toolkit website "
                              "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")

    def work_dir_info(self, biopype_or_sra, work_dir):
        if biopype_or_sra.lower() == 'biopype':
            print("\n The current working directory for BioPype is: " + work_dir)
            print("\t>To change BioPype's working directory, please run: path_helper.setup()")
            print("\t>The current working directory is where BioPype's "
                  "functions will look for target files and directories. It is "
                  "also where BioPype's functions will save data downloaded "
                  "from the SRA database. ")

        elif biopype_or_sra.lower() == 'sra':
            print("\nThe current Workspace Location for the SRA Toolkit is: " + work_dir)
            print("\t>For information on how to configure the SRA Toolkit Workspace "
                  "Location, please refer to either the BioPype manual (Chapter 8: "
                  "Software and Set-up), or the SRA Toolkit website "
                  "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")
        return None

    def setup(self, biopype_or_sra):
        """
        This function lets the user define the working directory for BioPype.
        : biopype_or_sra: string; 'BioPype' will let the user define the path for
        BioPype's working directory. 'sra' will let the user define the path for
        the SRA Toolkit's Workspace Location, which is where it downloads .sra files.
        :return: None
        """

        def rewrite_workspace_file(new_dir, workspace_file):
            """
            This function will replace the path currently given by a given
            workspace-name file with a new path.
            :param new_dir: string; the new path that will overwrite the old path.
            :param workspace_file: string; path to a file containing a workspace
            name. Will be a path to either sratoolworkspace.txt or biopypeworkdir.txt
            :return: None
            """
            if os.path.isdir(new_dir):
                with open(workspace_file, 'w') as w:
                    w.write(new_dir
                            + "\n The first line of this file (the above line^) is the "
                              "path to the working directory of the BioPype module.")
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

        # bio_dirname_file = os.path.join(current_dir, 'biopypeworkdir.txt')
        # sra_dirname_file = os.path.join(current_dir, 'sratoolworkspace.txt')

        if biopype_or_sra == 'biopype':
            working_dir = input("Please input a path to an existing folder to define BioPype's working directory: ")
            rewrite_workspace_file(working_dir, self.bio_dirname_file)
            BIODIR = working_dir

        if biopype_or_sra == 'sra':
            workspace_location = input("Please input the path for the SRA Toolkit's configured Workspace Location: ")
            rewrite_workspace_file(workspace_location, self.sra_dirname_file)
            SRADIR = workspace_location
        return None

if __name__ != "__main__":
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    biodirfile = os.path.join(cur_dir, 'biopypeworkdir.txt')
    with open(biodirfile, 'r') as f:
        # Use .strip() here to remove the newline character at the end of the string.
        BIODIR = f.readline().strip()

    SRADIR = 'pass'