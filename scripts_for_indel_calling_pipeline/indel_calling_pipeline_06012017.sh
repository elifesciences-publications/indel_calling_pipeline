#!/bin/sh




if [ "$#" -ne 6 ] 
  then
  echo "Usage: $0 sample_name full_path_to_job_dir mark full_path_to_ref_genome full_path_to_ref_genome_2bit min_cov_at_indel  " >&2
  echo -e "\nIMPORTANT: The BAM file for <sample_name> is assumed to be in <full_path_to_job_dir>/BAM and is assumed to have the following name format: <sample_name>_<mark>*rmdup*.bam" >&2
  echo -e "\n Software dependencies: bowtie2, bowtie2-build, samtools, GenomAnalysisTK.jar, AddOrReplaceReadGroups.jar, pull_indel_candidates_from_SAM_file.pl (in-house script), create_custom_fasta_with_indel_081516.pl (in-house script)" >&2
  exit 1
fi


###IMPORTANT-THESE PATHS WILL NEED TO BE MODIFIED IN ACCORDANCE WITH YOUR DIRECTORY STRUCTURE!!!###
path_to_gatk="/home/ars51/programs/executables/GenomeAnalysisTK.jar"
path_to_picard_reformatter="/home/ars51/programs/picard-tools-1.104/AddOrReplaceReadGroups.jar"
path_to_varscan="/home/ars51/programs/VarScan.v2.4.0.jar"
path_to_helper_scripts="/mnt/genetics/ScacheriLab/ars51/HOME1/SCRIPTS_REPO"
java_tmp_dir_for_varscan="/home/ars51/java_temp"
######


s=$1 ##sample name example:COLO205
job_dir=$2 ##working dir
mark=$3 ##example:H3K27ac
path_to_ref_genome=$4 #/home/ars51/genomes/hg19/hg19.fa
path_to_ref_genome_2bit=$5 #/home/ars51/genomes/hg19/hg19.2bit 
cov=$6 ##example:10

min_reads_at_alt=`echo $cov | perl -MList::Util -lane 'print List::Util::max(1,$F[0]*0.2)'`
##create directories for storing output files
mkdir -p $job_dir/{genome_files,text_files,fastq}

bam=`ls $job_dir/BAM/$s"_"$mark*rmdup*.bam`
if [[ ! -f $bam.bai ]]; then
    samtools index $bam
fi


step_size=`samtools view -F 4 $bam |awk '{print length($10)}' | sort -u -k 1 -n -r |head -1` ##get the read length
prefix=$s"_"$mark"_"$step_size"bp"
indels_file=$job_dir/text_files/$s"_"$mark"_candidate_indels_from_indel_realigned_bam.txt"


##reformat orig BAM file for further processing
java -Xmx4g -jar $path_to_picard_reformatter I=$bam O=$job_dir/BAM/$s"_"$mark.picard_reformatted.bam LB=$s"_"$mark PL=illumina PU=$s SM=$s
samtools index $job_dir/BAM/$s"_"$mark.picard_reformatted.bam

##create targets for indel realignment and perform realignment around candidate indels
java -Xmx4g -Djava.io.tmpdir=/home/ars51/java_temp -jar $path_to_gatk -T RealignerTargetCreator -R $path_to_ref_genome -I $job_dir/BAM/$s"_"$mark.picard_reformatted.bam -o $job_dir/text_files/forIndelRealigner.$s"_"$mark.intervals --allow_potentially_misencoded_quality_scores
java -Xmx4g -Djava.io.tmpdir=/home/ars51/java_temp -jar $path_to_gatk -T IndelRealigner -R $path_to_ref_genome -I $job_dir/BAM/$s"_"$mark.picard_reformatted.bam -targetIntervals $job_dir/text_files/forIndelRealigner.$s"_"$mark.intervals  -o $job_dir/BAM/indel_realigned.$s"_"$mark.bam --allow_potentially_misencoded_quality_scores
samtools view -o $job_dir/BAM/indel_realigned.$s"_"$mark.sam -F 4 $job_dir/BAM/indel_realigned.$s"_"$mark.bam

##aligned reads aligning with mismatches/indels are saved for later as a BAM file; reads that align perfectly to the reference genome are  saved in FASTQ format for subsequent realignment to the indel genome
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s"_"$mark.sam | grep -E  '.+AS:i:0.+' | perl -lane 'print "@".$F[0]."\n".$F[9]."\n+\n".$F[10]' > $job_dir/fastq/for_realignment_$s"_"$mark.fastq
awk '$6 !~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s"_"$mark.sam > $job_dir/BAM/for_mpileup_part1_indel_realigned.$s"_"$mark.sam

##pull candidate indels from re-aligned file
samtools view -b -o $job_dir/BAM/indel_realigned.$s"_"$mark.bam    -T $path_to_ref_genome $job_dir/BAM/for_mpileup_part1_indel_realigned.$s"_"$mark.sam
$path_to_helper_scripts/pull_indel_candidates_from_SAM_file.pl  < $job_dir/BAM/indel_realigned.$s"_"$mark.sam > $job_dir/text_files/$s"_"$mark"_candidate_indels_from_indel_realigned_bam.txt"


##create indel genome FASTA file
max_ext_on_right_side=`cat $indels_file| grep '-' | cut -f 3 | sed 's/-\(.*\)/\1/g' | awk '{print length($1)}' | sort -k 1 -n -r | head -1`
cat $indels_file | grep -v Chrom| sort | uniq -c | awk -v mr=$min_reads_at_alt '$1>=mr' | awk '{OFS="\t"}{print $2, $3, $4}' | awk -v s=$step_size -v e=$max_ext_on_right_side '{OFS="\t"}{print $1, $2-s, $2+s+e, $1"_"$2}' | awk '($2<0){$2=0}{OFS="\t"}{print $1, $2, $3, $4}'> $job_dir/text_files/$prefix.bed
cat $indels_file | grep -v Chrom | sort | uniq -c | awk -v mr=$min_reads_at_alt '$1>=mr' | awk '{OFS="\t"}{print $2, $3, $4}' | awk '{OFS="\t"}{print $1"_"$2,$3}' > $job_dir/text_files/for_grep_$prefix.bed
twoBitToFa  $path_to_ref_genome_2bit  -bed=$job_dir/text_files/$prefix.bed  $job_dir/genome_files/$prefix.fa
$path_to_helper_scripts/create_custom_fasta_with_indel_081516.pl   $job_dir/genome_files/$prefix.fa $step_size $job_dir/text_files/for_grep_$prefix.bed  $job_dir/genome_files/indel_genome_$prefix.fa
bowtie2-build $job_dir/genome_files/indel_genome_$prefix.fa  $job_dir/genome_files/indel_genome_$prefix

##align to indel genome
bowtie2 -x $job_dir/genome_files/indel_genome_$prefix -U $job_dir/fastq/for_realignment_$s"_"$mark.fastq -S $job_dir/BAM/indel_genome_$prefix.sam

##get IDs of reads that map perfectly to indel genome (those reads will be excluded from indel calling analyses because they also map perfectly to the reference genome)
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_genome_$prefix.sam | grep -E  '.+AS:i:0.+' | cut -f 1   > $job_dir/text_files/$prefix"_ids_of_reads_to_remove.txt"

##get IDs of reads that map perfectly to reference genome
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s"_"$mark.sam | grep -E  '.+AS:i:0.+' > $job_dir/BAM/indel_realigned_perfect_matches.$s"_"$mark.sam

##remove ambiguously mapping reads (those aligning to both reference and indel genome)
grep -F -vf $job_dir/text_files/$prefix"_ids_of_reads_to_remove.txt"  $job_dir/BAM/indel_realigned_perfect_matches.$s"_"$mark.sam > $job_dir/BAM/indel_realigned_filtered.$s"_"$mark.sam

##combine indel_realigned/filtered BAM and reads that do not perfectly align to reference genome (potential indel reads)
samtools view -b -o $job_dir/BAM/for_mpileup_part2_indel_realigned.$s"_"$mark.bam -T $path_to_ref_genome $job_dir/BAM/indel_realigned_filtered.$s"_"$mark.sam
samtools view -b -o $job_dir/BAM/for_mpileup_part1_indel_realigned.$s"_"$mark.bam -T $path_to_ref_genome $job_dir/BAM/for_mpileup_part1_indel_realigned.$s"_"$mark.sam
samtools merge -f  $job_dir/BAM/for_mpileup.indel_realigned_and_filtered.$s"_"$mark.bam $job_dir/BAM/for_mpileup_part1_indel_realigned.$s"_"$mark.bam $job_dir/BAM/for_mpileup_part2_indel_realigned.$s"_"$mark.bam

##run mpileup command and then run Varscan for indel calling on merged BAM
samtools mpileup -f $path_to_ref_genome -o $job_dir/text_files/$s"_"$mark"_mpileup" $job_dir/BAM/for_mpileup.indel_realigned_and_filtered.$s"_"$mark.bam
java -Xmx4g -Djava.io.tmpdir=$java_tmp_dir_for_varscan -jar $path_to_varscan mpileup2indel $job_dir/text_files/$s"_"$mark"_mpileup" --min-var-freq 0.2 --min-coverage $cov> $job_dir/text_files/$s"_"$mark"_mpileup2indel_varscan"

echo "Pipeline finished. The main output files of interest are: $job_dir/text_files/"$s"_"$mark"_mpileup2indel_varscan and  $job_dir/BAM/for_mpileup.indel_realigned_and_filtered."$s"_$mark.bam.  All other intermediate files were kept for troubleshooting purposes only."