import subprocess

from BioPype.workspaces.dirpaths import DirPathsHelper

path_helper = DirPathsHelper()
_BIODIR = path_helper._BIODIR
_SRADIR = path_helper._SRADIR

class QiimeHelper:
    def __init__(self):
        return None

    def create_data_artifact(self, semanticType, inputPath, outputFilePath, sourceFormat):
        """Create a qiime2 artifact from sequencing data.

        :param semanticType: the qiime2 semantic-type of the created artifact.
        :param inputPath: Path to the file or directory that should be
        imported. (e.g., a qiime2 Fastq Manifest file.)
        :param outputFilePath: Path where output artifact should be written.
        :param sourceFormat: The format of the data to be imported. If not
        provided, data must be in the format expected by the semantic type
        provided via semanticType.
        :return: The path to the output .qza artifact
        """
        subprocess.run([
            'qiime', 'tools', 'import', '--type', semanticType, '--input-path',
            inputPath, '--output-path', outputFilePath, '--source-format', sourceFormat])
        return outputFilePath

    def custom_create_data_artifact(self, *args):
        """Create a qiime2 artifact from sequencing data.

        This method has more flexibility in specifying arguments than
        QiimeHelper.create_data_artifact(), but requires that the user have
        greater knowledge of the arguments that qiime allows for the
        tools.import command.
        :param args:
        Options:
        --type TEXT                The semantic type of the artifact that will be
                                     created upon importing. Use --show-importable-
                                     types to see what importable semantic types are
                                     available in the current deployment.  [required]
        --input-path PATH          Path to file or directory that should be
                                     imported.  [required]
        --output-path PATH         Path where output artifact should be written.
                                     [required]
        --source-format TEXT       The format of the data to be imported. If not
                                     provided, data must be in the format expected by
                                     the semantic type provided via --type.
        --show-importable-types    Show the semantic types that can be supplied to
                                     --type to import data into an artifact.
        --show-importable-formats  Show formats that can be supplied to --source-
                                     format to import data into an artifact.
        --help                     Show this message and exit.
        :return: None
        """
        # TODO: Test this method. Need to verify if unpacking is happening correctly.
        subprocess.run(['qiime', args])
        return None

    def summarize_data_artifact(self, input_artifact_path, output_visualization_path):
        """ Generates a visuzalization to summarize a data artifact.

        :param input_artifact_path:
        :param output_visualization_path:
        :return:
        """
        subprocess.run(['qiime', 'demux', 'summarize',
                        '--i-data',input_artifact_path,
                        '--o-visualization', output_visualization_path])
        return output_visualization_path

    def view_qzv(self, qzv_path):
        """Views a visualization object.

        :param qzv_path:
        :return:
        """
        subprocess.run(['qiime', 'tools', 'view', qzv_path])
        return None

    def denoise_paired_data(self, *args):
        """
        Options for *args:
        --i-demultiplexed-seqs ARTIFACT PATH SampleData[PairedEndSequencesWithQuality]
                                          The paired-end demultiplexed sequences to be
                                          denoised.  [required]
        --p-trunc-len-f INTEGER         Position at which forward read sequences
                                          should be truncated due to decrease in
                                          quality. This truncates the 3' end of the of
                                          the input sequences, which will be the bases
                                          that were sequenced in the last cycles.
                                          Reads that are shorter than this value will
                                          be discarded. After this parameter is
                                          applied there must still be at least a 20
                                          nucleotide overlap between the forward and
                                          reverse reads. If 0 is provided, no
                                          truncation or length filtering will be
                                          performed  [required]
        --p-trunc-len-r INTEGER         Position at which reverse read sequences
                                          should be truncated due to decrease in
                                          quality. This truncates the 3' end of the of
                                          the input sequences, which will be the bases
                                          that were sequenced in the last cycles.
                                          Reads that are shorter than this value will
                                          be discarded. After this parameter is
                                          applied there must still be at least a 20
                                          nucleotide overlap between the forward and
                                          reverse reads. If 0 is provided, no
                                          truncation or length filtering will be
                                          performed  [required]
        --p-trim-left-f INTEGER         Position at which forward read sequences
                                          should be trimmed due to low quality. This
                                          trims the 5' end of the input sequences,
                                          which will be the bases that were sequenced
                                          in the first cycles.  [default: 0]
        --p-trim-left-r INTEGER         Position at which reverse read sequences
                                          should be trimmed due to low quality. This
                                          trims the 5' end of the input sequences,
                                          which will be the bases that were sequenced
                                          in the first cycles.  [default: 0]
        --p-max-ee FLOAT                Reads with number of expected errors higher
                                          than this value will be discarded.
                                          [default: 2.0]
        --p-trunc-q INTEGER             Reads are truncated at the first instance of
                                          a quality score less than or equal to this
                                          value. If the resulting read is then shorter
                                          than `trunc_len_f` or `trunc_len_r`
                                          (depending on the direction of the read) it
                                          is discarded.  [default: 2]
        --p-chimera-method [none|consensus|pooled]
                                          The method used to remove chimeras. "none":
                                          No chimera removal is performed. "pooled":
                                          All reads are pooled prior to chimera
                                          detection. "consensus": Chimeras are
                                          detected in samples individually, and
                                          sequences found chimeric in a sufficient
                                          fraction of samples are removed.  [default:
                                          consensus]
        --p-min-fold-parent-over-abundance FLOAT
                                          The minimum abundance of potential parents
                                          of a sequence being tested as chimeric,
                                          expressed as a fold-change versus the
                                          abundance of the sequence being tested.
                                          Values should be greater than or equal to 1
                                          (i.e. parents should be more abundant than
                                          the sequence being tested). This parameter
                                          has no effect if chimera_method is "none".
                                          [default: 1.0]
        --p-n-threads INTEGER           The number of threads to use for
                                          multithreaded processing. If 0 is provided,
                                          all available cores will be used.  [default:
                                          1]
        --p-n-reads-learn INTEGER       The number of reads to use when training the
                                          error model. Smaller numbers will result in
                                          a shorter run time but a less reliable error
                                          model.  [default: 1000000]
        --p-hashed-feature-ids / --p-no-hashed-feature-ids
                                          If true, the feature ids in the resulting
                                          table will be presented as hashes of the
                                          sequences defining each feature. The hash
                                          will always be the same for the same
                                          sequence so this allows feature tables to be
                                          merged across runs of this method. You
                                          should only merge tables if the exact same
                                          parameters are used for each run.  [default:
                                          True]
        --o-table ARTIFACT PATH FeatureTable[Frequency]
                                          The resulting feature table.  [required if
                                          not passing --output-dir]
        --o-representative-sequences ARTIFACT PATH FeatureData[Sequence]
                                          The resulting feature sequences. Each
                                          feature in the feature table will be
                                          represented by exactly one sequence, and
                                          these sequences will be the joined paired-
                                          end sequences.  [required if not passing
                                          --output-dir]
        --o-denoising-stats ARTIFACT PATH SampleData[DADA2Stats]
                                          [required if not passing --output-dir]
        --output-dir DIRECTORY          Output unspecified results to a directory
        --cmd-config PATH               Use config file for command options
        --verbose                       Display verbose output to stdout and/or
                                          stderr during execution of this action.
                                          [default: False]
        --quiet                         Silence output if execution is successful
                                          (silence is golden).  [default: False]
        --citations                     Show citations and exit.
        --help                          Show this message and exit.

        :return:
        """
        subprocess.run(['qiime', 'dada2', 'denoise-paired', args])
        return None



# qhelper = QiimeHelper()
# inp = '/Volumes/Gniot_Backup_Drive/repos/my_test/cd_young_manifest.csv'
# outp = '/Volumes/Gniot_Backup_Drive/repos/my_test/CDY-paired-end-demux.qza'
# qhelper.create_data_artifact('SampleData[SequencesWithQuality]', inp, outp, 'PairedEndFastqManifestPhred33')
