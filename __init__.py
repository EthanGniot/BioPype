# author: Ethan Gniot
import os
import warnings
import workspaces.dirpaths


def work_dir_info(current_working_dir):
    print("The current working directory for BioPype is: "
          + current_working_dir
          + "\n\t>To change BioPype's working directory, please run: "
            "BioPype.workspaces.setup() \n\t>The current working directory is where BioPype's functions "
            "will look for target files and directories. It is also where "
            "BioPype's functions will save data downloaded from the SRA "
            "database. ")
    return None


# Get the path to the directory where the BioPype code located.
# Need this variable so that setup() and this __init__ file can find the
# text files containing the paths to the BioPype and SRA Toolkit workspaces.
# biopype_code_dir = os.path.dirname(os.path.realpath(__file__))

# Get paths to the files containing the workspace names.
# bio_dirname_file = os.path.join(biopype_code_dir, 'biopypeworkdir.txt')
# sra_dirname_file = os.path.join(biopype_code_dir, 'sratoolworkspace.txt')

bio_dirname_file = os.path.join('workspaces', 'biopypeworkdir.txt')
sra_dirname_file = os.path.join('workspaces', 'sratoolworkspace.txt')

with open(bio_dirname_file, 'r') as f:
    # Use .strip() here to remove the newline character at the end of the string.
    directory = f.readline().strip()
    if directory:
        work_dir_info(directory)
    else:
        work_dir_info(os.getcwd())

with open(sra_dirname_file, 'r') as f:
    directory = f.readline().strip()
    if directory:
        print("The current Workspace Location for the SRA Toolkit is: " + directory)
        print("\t>For information on how to configure the SRA Toolkit Workspace Location, "
          "please refer to either the BioPype manual (Chapter 8: Software and "
          "Set-up), or the SRA Toolkit website "
          "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")
    else:
        warnings.warn("***THERE IS CURRENTLY NO PATH DEFINED FOR THE SRA TOOLKIT WORKSPACE*** "
                      "BioPype relies on functionality from NCBI's "
                      "SRA Toolkit, which requires configuration of the "
                      "toolkit's Workspace Location after it is installed. For "
                      "information on how to configure the SRA Toolkit Workspace"
                      " Location, please refer to either the BioPype manual "
                      "(Chapter 8: Software and Set-up), or the SRA Toolkit website "
                      "(https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std#s-4).")
