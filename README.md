# BioPype

_This tool is currently a work-in-progress and does not represent the final 
product. Files and features are in varying stages of completion and will likely
change throughout the development process._

This tutorial aims to improve a newcomer's general understanding of bioinformatics by...
  1) Defining technical terms commonly used in bioinformatics, 
  2) Providing a collection of useful resources for learning about bioinformatics, and 
  3) Demonstrating how Python can be used to answer specific research questions
   by combining existing bioinformatics tools.

+++++

Functions and classes for a gut microbiome analysis pipeline can be found under 
the `cmds` directory. The `workspaces` directory contains classes and files 
used for configuring the BioPype work environment.

**Current operational files:**
* `cmds.runtable.py`: Defines the `RunTable` class for handling RunInfo Tables 
downloaded from the Sequence Read Archive database. 
* `cmds.downloadsra.py`: Defines functions for handling the download and conversion 
of .sra files to .fastq files.
* `cmds.qiimehelper.py`: Defines the `QiimeHelper` class, which handles the
the .qza and .qzv qiime artifacts used in the relative abundance analysis.
* `cmds.qualitycontrol.py`: Defines functions for performing quality control
on the raw sequencing data. 
* `workspaces.dirpaths.py`: Defines the `DirPathsHelper` class, which handles
initialization code for BioPype as well as methods for defining the working
directories for BioPype and the SRA Toolkit.

The current state of the accompanying manual/tutorial can be found in the 
`manual` directory as `main_manual.pdf`.

**Current sections under development:**
* Chapter 1: Foreword
* Chapter 2: Microbiome Analysis
* Chapter 3: Software and Set-up
* Chapter 4: The Dataset
* Appendix A: Web Resources and Explanations
* Appendix B: Python Resources

**Future sections**
* How to Find Tools
* Relative Abundance Analysis with QIIME2
* Predict ORFs
* Create Non-redundant Gene Sets
* Align Genes
* Get GenBank Accession Numbers
* Find COG Functional Class
* Gene Ontology (GO) Classification


**(Initial) Tutorial Pipeline Structure:**
![biopype_flowchart](https://user-images.githubusercontent.com/30661548/39219342-6ab16580-47ef-11e8-8c0d-5770d20b9b97.png)

