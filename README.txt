
INDEL CALLING PIPELINE STEPS
1. Software dependencies must be installed on your system and accessible through the $PATH variable (installing these through miniconda is recommended if your default samtools is version 1.0+. samtools v1.0+ will not work with this pipeline):

 bowtie2
 bowtie2-build
 samtools 0.1.19

2. These scripts/executables must be placed in a directory that you will specify when running the pipeline:

GenomAnalysisTK.jar
VarScan.v2.4.0.jar
AddOrReplaceReadGroups.jar
pull_indel_candidates_from_SAM_file.pl (in-house script)
create_custom_fasta_with_indel_081516.pl (in-house script)

3. Create the directory 'BAM' in your working directory (accessible through /full/path/to/job/dir/BAM). The BAM file for <sample_name> is assumed to be in /full/path/to/job/dir/BAM and is assumed to have the following name format: <sample_name>*rmdup*.bam"
4. Run indel_calling_pipeline.sh. Here is an example command: indel_calling_pipeline.sh  C29_H3K27ac /home/ars51/my_analysis /home/ars51/genomes/hg19/hg19 12 /home/ars51/scripts_for_pipeline

