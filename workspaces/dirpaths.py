import os


def setup(biopype_or_sra):
    """
    This function lets the user define the working directory for BioPype.
    : biopype_or_sra: string; 'biopype' will let the user define the path for
    BioPype's working directory. 'sra' will let the user define the path for
    the SRA Toolkit's Workspace Location, which is where it downloads .sra files.
    :return: None
    """

    def rewrite_workspace_file(new_dir, workspace_file):
        """
        This function will replace the path currently given by a given workspace
         name file with a new one.
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
        rewrite_workspace_file(working_dir, bio_dirname_file)
        BIODIR = working_dir

    if biopype_or_sra == 'sra':
        workspace_location = input("Please input the path for the SRA Toolkit's configured Workspace Location: ")
        rewrite_workspace_file(workspace_location, sra_dirname_file)
        SRADIR = workspace_location
    return None


current_dir = os.path.join(os.path.dirname(__file__))

bio_dirname_file = os.path.join(current_dir, 'biopypeworkdir.txt')
sra_dirname_file = os.path.join(current_dir, 'sratoolworkspace.txt')

with open(bio_dirname_file, 'r') as f:
    BIODIR = f.readline().strip()

with open(sra_dirname_file, 'r') as f:
    SRADIR = f.readline().strip()

