import os
import BioPype.cmds.runtable as rt
from BioPype.cmds.downloadsra import DownloadHelper

runinfotable_path = '/Volumes/Gniot_Backup_Drive/repos/my_test/runinfotable-human-micro-IBD.txt'
my_table = rt.RunTable(runinfotable_path)

print(my_table.column_names)
# Narrow down samples
wgs_samples = my_table.filter_data("Assay_Type", '==', 'AMPLICON')
print(wgs_samples.df['host_disease'])
#
# uc_samples = wgs_samples.filter_data("host_disease", '==', 'ulcerative colitis')
uc_samples = wgs_samples.filter_data("host_disease", '==', 'UC')
# cd_samples = wgs_samples.filter_data("host_disease", '==', "Crohn''s disease")
cd_samples = wgs_samples.filter_data("host_disease", '==', "CD")

uc_young = uc_samples.filter_data("host_age", '<=', 21)
uc_old = uc_samples.filter_data("host_age", '>=', 40)
cd_young = cd_samples.filter_data("host_age", '<=', 21)
cd_old = cd_samples.filter_data("host_age", '>=', 40)
#

uc_young_nums = uc_young.get_accession_numbers()
uc_old_nums = uc_old.get_accession_numbers()
cd_young_nums = cd_young.get_accession_numbers()
cd_old_nums = cd_old.get_accession_numbers()

print(len(wgs_samples.df), len(uc_samples.df), len(cd_samples.df))

###############################

# Create DownloadHelper objects.
uc_young_obj = DownloadHelper(uc_young_nums, 2, threads=4)
uc_old_obj = DownloadHelper(uc_old_nums, 2, threads=4)
cd_young_obj = DownloadHelper(cd_young_nums, 2, threads=4)
cd_old_obj = DownloadHelper(cd_old_nums, 2, threads=4)

# Prefilter: Download SRA and convert to FASTQ
uc_young_obj.download_paired_end_data('uc_young', 'uc_young')
uc_old_obj.download_paired_end_data('uc_old', 'uc_old')
cd_young_obj.download_paired_end_data('cd_young', 'cd_young')
cd_old_obj.download_paired_end_data('cd_old', 'cd_old')

# Create manifest files for the untrimmed data
cdy_man = cd_young_obj.create_manifest_file('cd_young_manifest.csv')
cdo_man = cd_old_obj.create_manifest_file('cd_old_manifest.csv')
uco_man = uc_old_obj.create_manifest_file('uc_old_manifest.csv')
ucy_man = uc_young_obj.create_manifest_file('uc_young_manifest.csv')

# Create data artifacts using the manifest files.
from BioPype.cmds.qiimehelper import QiimeHelper
qhelper = QiimeHelper

DownloadHelper.check_dir_exists('/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data')
DownloadHelper.check_dir_exists('/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data.summaries')
data_summaries = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data.summaries'

cdy_demux = qhelper.create_data_artifact(
    semanticType='SampleData[SequencesWithQuality]',
    inputPath=cdy_man,
    outputFilePath='/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data/cdy-demux.qza',
    sourceFormat='PairedEndFastqManifestPhred33')

cdo_demux = qhelper.create_data_artifact(
    semanticType='SampleData[SequencesWithQuality]',
    inputPath=cdo_man,
    outputFilePath='/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data/cdo-demux.qza',
    sourceFormat='PairedEndFastqManifestPhred33')

uco_demux = qhelper.create_data_artifact(
    semanticType='SampleData[SequencesWithQuality]',
    inputPath=uco_man,
    outputFilePath='/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data/uco-demux.qza',
    sourceFormat='PairedEndFastqManifestPhred33')

ucy_demux = qhelper.create_data_artifact(
    semanticType='SampleData[SequencesWithQuality]',
    inputPath=ucy_man,
    outputFilePath='/Volumes/Gniot_Backup_Drive/repos/my_test/data/qiime2/objects.data/ucy-demux.qza',
    sourceFormat='PairedEndFastqManifestPhred33')

# Now produce a visualization of the demux’ing. This will produce a file called
#  demux.qza, which can be viewed on the QIIME2 viewer at https://view.qiime2.org/
o = os.path.join(data_summaries, 'cdy-demux.qzv')
cdy_demux_qzv = qhelper.summarize_data_artifact(cdy_demux, o)
qhelper.view_qzv(cdy_demux_qzv)

o = os.path.join(data_summaries, 'cdo-demux.qzv')
cdo_demux_qzv = qhelper.summarize_data_artifact(cdo_demux, o)
qhelper.view_qzv(cdo_demux_qzv)

o = os.path.join(data_summaries, 'uco-demux.qzv')
uco_demux_qzv = qhelper.summarize_data_artifact(uco_demux, o)
qhelper.view_qzv(uco_demux_qzv)

o = os.path.join(data_summaries, 'ucy-demux.qzv')
ucy_demux_qzv = qhelper.summarize_data_artifact(ucy_demux, o)
qhelper.view_qzv(ucy_demux_qzv)
#


# Denoise the data. Choose options for the denoise function based on the
# quality of the data we see in the .qzv visualizations. Right now, the output
# files will be created in the current working directory.
qhelper.denoise_paired_data(
    '--i-demultiplexed-seqs', cdy_demux,
    '--p-trunc-len-f', '0',
    '--p-trunc-len-r', '0',
    '--p-n-threads', '0',
    '--o-representative-sequences', 'cdy-repseqs.qza',
    '--o-table', 'cdy-table.qza',
    '--o-denoising-stats', 'cdy-stats.qza')

qhelper.denoise_paired_data(
    '--i-demultiplexed-seqs', cdo_demux,
    '--p-trunc-len-f', '0',
    '--p-trunc-len-r', '0',
    '--p-n-threads', '0',
    '--o-representative-sequences', 'cdo-repseqs.qza',
    '--o-table', 'cdo-table.qza',
    '--o-denoising-stats', 'cdo-stats.qza')

qhelper.denoise_paired_data(
    '--i-demultiplexed-seqs', ucy_demux,
    '--p-trunc-len-f', '0',
    '--p-trunc-len-r', '0',
    '--p-n-threads', '0',
    '--o-representative-sequences', 'ucy-repseqs.qza',
    '--o-table', 'ucy-table.qza',
    '--o-denoising-stats', 'ucy-stats.qza')

qhelper.denoise_paired_data(
    '--i-demultiplexed-seqs', uco_demux,
    '--p-trunc-len-f', '0',
    '--p-trunc-len-r', '0',
    '--p-n-threads', '0',
    '--o-representative-sequences', 'uco-repseqs.qza',
    '--o-table', 'uco-table.qza',
    '--o-denoising-stats', 'uco-stats.qza')

# To visualize the ____-stats.qza file, use the following command-line command:
#     qiime metadata tabulate \
#       --m-input-file stats-dada2.qza \
#       --o-visualization stats-dada2.qzv



#This section could go before the manifest files, but for the sake of testing i'm putting it here.
##############
# Get quality control reports for the fastq files.
import BioPype.cmds.qualitycontrol as qc
qc.get_fastqc_reports()

# Get multiqc summary of the fastqc files for each experimental group.
qc.summarize_reports('summary')
################

# Filter the fastq files for each experimental group.
target1 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/cd_young'
target2 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/cd_old'
target3 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/uc_young'
target4 = '/Volumes/Gniot_Backup_Drive/repos/my_test/data/grouped_fastq/uc_old'
qc.run_trim_galore(target_dir=target1)
qc.run_trim_galore(target2)
qc.run_trim_galore(target3)
qc.run_trim_galore(target4)

# Get multiqc summary of the filtered fastqc files for each experimental group
qc.summarize_reports('summary.trimmed', 'grouped_fastq.trimmed')

######
# Take some time to confirm that the filtering improved the data quality of the fastq files.
######

cd_young_obj.create_metadata_file(my_table, 'cd_young.tsv')
cd_old_obj.create_metadata_file(my_table, 'cd_old_metadata.tsv')
uc_old_obj.create_metadata_file(my_table, 'uc_old_metadata.tsv')
uc_young_obj.create_metadata_file(my_table, 'uc_young_metadata.tsv')




