import os
import glob
import pandas as pd
import BioPype.cmds.runtable as rt
from BioPype.cmds.downloadsra import DownloadHelper

runinfotable_path = '/Volumes/Gniot_Backup_Drive/repos/my_test/runinfotable-human-micro-IBD.txt'
my_table = rt.RunTable(runinfotable_path)

# Narrow down samples
amp_samples = my_table.filter_data("Assay_Type", '==', 'AMPLICON')

uc_samples = amp_samples.filter_data("host_disease", '==', 'UC')
cd_samples = amp_samples.filter_data("host_disease", '==', "CD")

uc_young = uc_samples.filter_data("host_age", '<=', 21)
uc_old = uc_samples.filter_data("host_age", '>=', 40)
cd_young = cd_samples.filter_data("host_age", '<=', 21)
cd_old = cd_samples.filter_data("host_age", '>=', 40)


uc_young_nums = uc_young.get_accession_numbers()
uc_old_nums = uc_old.get_accession_numbers()
cd_young_nums = cd_young.get_accession_numbers()
cd_old_nums = cd_old.get_accession_numbers()


# Create DownloadHelper objects.
cd_young_obj = DownloadHelper(cd_young_nums, 5, threads=4)
cd_old_obj = DownloadHelper(cd_old_nums, 5, threads=4)
uc_old_obj = DownloadHelper(uc_old_nums, 5, threads=4)
uc_young_obj = DownloadHelper(uc_young_nums, 5, threads=4)

# Prefilter: Download SRA and convert to FASTQ
cd_young_obj.download_paired_end_data('cd_young', 'cd_young')
cd_old_obj.download_paired_end_data('cd_old', 'cd_old')
uc_young_obj.download_paired_end_data('uc_young', 'uc_young')
uc_old_obj.download_paired_end_data('uc_old', 'uc_old')

### ONLY THE CODE UP TO THIS POINT HAS BEEN DOCUMENTED IN THE MANUAL SO FAR ##

# Create manifest files for the untrimmed data
manifest = DownloadHelper.create_manifest_file(
    'sample-manifest.csv',
    cd_young_obj.fin_download_dict,
    cd_old_obj.fin_download_dict,
    uc_young_obj.fin_download_dict,
    uc_old_obj.fin_download_dict)

cdy_man = cd_young_obj.method_create_manifest_file('cd_young_manifest.csv')
cdo_man = cd_old_obj.method_create_manifest_file('cd_old_manifest.csv')
uco_man = uc_old_obj.method_create_manifest_file('uc_old_manifest.csv')
ucy_man = uc_young_obj.method_create_manifest_file('uc_young_manifest.csv')

# Create metadata files for the data.
cdy_meta = cd_young_obj.create_metadata_file(my_table, 'cd_young.tsv')
cdo_meta = cd_old_obj.create_metadata_file(my_table, 'cd_old.tsv')
ucy_meta = uc_young_obj.create_metadata_file(my_table, 'uc_young.tsv')
uco_meta = uc_old_obj.create_metadata_file(my_table, 'uc_old.tsv')

# Merge the metadata files into one.
results = pd.DataFrame([])

for counter, file in enumerate(glob.glob("*.tsv")):
    namedf = pd.read_csv(file, sep='\t', header=0,)
    results = results.append(namedf)
results.to_csv('test_metadata.csv', sep='\t', index=False)

# Create data artifacts using the manifest files.
from BioPype.cmds.qiimehelper import QiimeHelper
qhelper = QiimeHelper()

DownloadHelper.check_dir_exists('/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data')
DownloadHelper.check_dir_exists('/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data.summaries')
data_summaries = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data.summaries'

demux = qhelper.create_data_artifact(
    semanticType='SampleData[PairedEndSequencesWithQuality]',
    inputPath='sample-manifest.csv',
    outputFilePath='/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data/demux.qza',
    sourceFormat='PairedEndFastqManifestPhred33')

# Now produce a visualization of the demux’ing. This will produce a file called
#  demux.qza, which can be viewed on the QIIME2 viewer at https://view.qiime2.org/
o = os.path.join(data_summaries, 'demux.qzv')
demux_qzv = qhelper.summarize_data_artifact(demux, o)


# Denoise the data. Choose options for the denoise function based on the
# quality of the data we see in the .qzv visualizations. Right now, the output
# files will be created in the current working directory.
qhelper.denoise_paired_data(
    'qiime', 'dada2', 'denoise-paired',
    '--i-demultiplexed-seqs', demux,
    '--p-trunc-len-f', '0',
    '--p-trunc-len-r', '0',
    '--p-n-threads', '0',
    '--o-representative-sequences', 'rep-seqs.qza',
    '--o-table', 'feat-table.qza',
    '--o-denoising-stats', 'stats.qza', '--verbose')


# To visualize the stats.qza file, use the following command-line command:
#     qiime metadata tabulate \
#       --m-input-file stats-dada2.qza \
#       --o-visualization stats-dada2.qzv

# Now it's time to do a taxonomic analysis.
# Get the Naive Bayes classifier that has been trained to recognize these data.
# subprocess.run(['wget', '-O', '"gg-13-8-99-515-806-nb-classifier.qza"', '"https://data.qiime2.org/2018.4/common/gg-13-8-99-515-806-nb-classifier.qza"'])
classifier_path = "/Volumes/Gniot_Backup_Drive/repos/my_test/gg-13-8-99-515-806-nb-classifier.qza"

# Execute the taxonomic classification.
tc = qhelper.run_taxonomic_classification(
    classifier_path=classifier_path,
    rep_seqs_filepath="rep-seqs.qza",
    output_filepath="tax-class.qza")

# Produce a visualization of the tax-class.qza artifact.
vis_tax = qhelper.vis_tax_class(
    input_filepath='tax-class.qza',
    output_filepath='tax-class.qzv')

# View the taxonomy classification table.
# qhelper.view_qzv(qzv_path=cdy_vis_tax)

# Plot the taxonomy classification in a barplot.
bar = qhelper.create_barplot(input_feat_table='feat-table.qza', input_tax_class=tc, metadata_filepath='test_metadata.csv', output_filepath='taxa-bar-plots.qzv')

# View the taxa barplot
# qhelper.view_qzv(qzv_path=cdy_bar)


#This section could go before the manifest files, but for the sake of testing i'm putting it here.
##############
# Get quality control reports for the fastq files.
# import BioPype.cmds.qualitycontrol as qc
# qc.get_fastqc_reports()

# Get multiqc summary of the fastqc files for each experimental group.
# qc.summarize_reports('summary')
################

# Filter the fastq files for each experimental group.
# target1 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/cd_young'
# target2 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/cd_old'
# target3 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/uc_young'
# target4 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/uc_old'
# qc.run_trim_galore(target_dir=target1)
# qc.run_trim_galore(target2)
# qc.run_trim_galore(target3)
# qc.run_trim_galore(target4)

# Get multiqc summary of the filtered fastqc files for each experimental group
# qc.summarize_reports('summary.trimmed', 'grouped_fastq.trimmed')

######
# Take some time to confirm that the filtering improved the data quality of the fastq files.
######
