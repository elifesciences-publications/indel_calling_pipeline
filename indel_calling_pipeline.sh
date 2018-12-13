#!/bin/sh

function error_exit {
    echo "$1" >&2   ## Send message to stderr. Exclude >&2 if you don't want it that way.
    exit "${2:-1}"  ## Return a code specified by $2 or 1 by default.
}

if [ "$#" -ne 5 ] 
  then
echo "Usage: $0 full_sample_name full/path/to/job/dir full/path/to/genome/prefix (requires the following file types: .fa, .2bit, and .bt2)  min_cov_at_indel full/path/to/scripts " >&2
echo -e "\nExample: indel_calling_pipeline.sh  C29_H3K27ac /home/ars51/my_analysis /home/ars51/genomes/hg19/hg19 12 /home/ars51/scripts_for_pipeline"
echo -e "\nMake sure the following executables/scripts are available in full/path/to/scripts:"
echo -e "bowtie2\nbowtie2-build\nsamtools\nGenomAnalysisTK.jar\nAddOrReplaceReadGroups.jar\npull_indel_candidates_from_SAM_file.pl (in-house script)\ncreate_custom_fasta_with_indel_081516.pl (in-house script)" >&2
echo -e "\nIMPORTANT: The BAM file for <sample_name> is assumed to be in /full/path/to/job/dir/BAM and is assumed to have the following name format: <sample_name>*rmdup*.bam" >&2
  exit 1
fi

###INPUT PARAMETERS###
s=$1 ##sample name example:COLO205
job_dir=$2 ##working dir
path_to_genome_file_prefix=$3 #/home/ars51/genomes/hg19/hg19.fa
cov=$4 ##example:10
full_path_to_scripts=$5 ##/home/ars51/my_scripts



path_to_gatk="$full_path_to_scripts/GenomeAnalysisTK.jar"
path_to_picard_reformatter="$full_path_to_scripts/AddOrReplaceReadGroups.jar"
path_to_varscan="$full_path_to_scripts/VarScan.v2.4.0.jar"
path_to_helper_scripts=$full_path_to_scripts
java_tmp_dir_for_varscan=$job_dir/java_temp
path_to_ref_genome_2bit=$path_to_genome_file_prefix.2bit
path_to_ref_genome=$path_to_genome_file_prefix.fa
[ ! -d $java_tmp_dir_for_varscan ] && mkdir -p $java_tmp_dir_for_varscan



##At least 20% of the total reads have to map to the alternate allele
min_reads_at_alt=`echo $cov | perl -MList::Util -lane 'print List::Util::max(1,$F[0]*0.2)'`

##create directories for storing output files
mkdir -p $job_dir/{genome_files,text_files,fastq}

[ `find $job_dir/BAM/ -name $s*rmdup*.bam | wc -l` ==  0 ] &&  error_exit "BAM file doesn't exist in $job_dir/BAM. Make sure the file name has the proper format. See README.txt"

echo "Indexing BAM file"
bam=`ls $job_dir/BAM/$s*rmdup*.bam`
if [[ ! -f $bam.bai ]]; then
	
    samtools index $bam
fi

echo "Reformatting original BAM file for further processing"
java -Xmx4g -jar $path_to_picard_reformatter I=$bam O=$job_dir/BAM/$s.picard_reformatted.bam LB=$s PL=illumina PU=$s SM=$s


echo "Indexing reformatted BAM file"
samtools index $job_dir/BAM/$s.picard_reformatted.bam

echo "Create targets for indel realignment and perform realignment around candidate indels"
java -Xmx4g -Djava.io.tmpdir=$java_tmp_dir_for_varscan -jar $path_to_gatk -T RealignerTargetCreator -R $path_to_ref_genome -I $job_dir/BAM/$s.picard_reformatted.bam -o $job_dir/text_files/forIndelRealigner.$s.intervals --allow_potentially_misencoded_quality_scores
java -Xmx4g -Djava.io.tmpdir=$java_tmp_dir_for_varscan -jar $path_to_gatk -T IndelRealigner -R $path_to_ref_genome -I $job_dir/BAM/$s.picard_reformatted.bam -targetIntervals $job_dir/text_files/forIndelRealigner.$s.intervals  -o $job_dir/BAM/indel_realigned.$s.bam --allow_potentially_misencoded_quality_scores
samtools view -o $job_dir/BAM/indel_realigned.$s.sam -F 4 $job_dir/BAM/indel_realigned.$s.bam

echo "Save reads that align perfectly to the ref genome as a FASTQ file for further realignmeht to the indel genome. Save reads aligning with mismatches/indels as a BAM file"
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s.sam | grep -E  '.+AS:i:0.+' | perl -lane 'print "@".$F[0]."\n".$F[9]."\n+\n".$F[10]' > $job_dir/fastq/for_realignment_$s.fastq
awk '$6 !~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s.sam > $job_dir/BAM/for_mpileup_part1_indel_realigned.$s.sam

echo "Pull candidate indels from IndelRealigner realigned file"
samtools view -b -o $job_dir/BAM/indel_realigned.$s.bam    -T $path_to_ref_genome $job_dir/BAM/for_mpileup_part1_indel_realigned.$s.sam
indels_file=$job_dir/text_files/$s"_candidate_indels_from_indel_realigned_bam.txt"
$path_to_helper_scripts/pull_indel_candidates_from_SAM_file.pl  < $job_dir/BAM/indel_realigned.$s.sam > $indels_file

echo "Get the read length in "$bam". This will be used to determine the size of each indel genome sequence"
step_size=`samtools view -F 4 $bam |awk '{print length($10)}' | sort -u -k 1 -n -r |head -1` 
prefix=$s"_"$step_size"bp"


echo "Create indel genome FASTA"
#max extention is equal to size of max indel
max_ext_on_right_side=`cat $indels_file| grep '-' | cut -f 3 | sed 's/-\(.*\)/\1/g' | awk '{print length($1)}' | sort -k 1 -n -r | head -1`
#exclude indels that don't have a minimum number of reads mapping to alt allele
cat $indels_file | grep -v Chrom| sort | uniq -c | awk -v mr=$min_reads_at_alt '$1>=mr' | awk '{OFS="\t"}{print $2, $3, $4}' | awk -v s=$step_size -v e=$max_ext_on_right_side '{OFS="\t"}{print $1, $2-s, $2+s+e, $1"_"$2}' | awk '($2<0){$2=0}{OFS="\t"}{print $1, $2, $3, $4}'> $job_dir/text_files/$prefix.bed
cat $indels_file | grep -v Chrom | sort | uniq -c | awk -v mr=$min_reads_at_alt '$1>=mr' | awk '{OFS="\t"}{print $2, $3, $4}' | awk '{OFS="\t"}{print $1"_"$2,$3}' > $job_dir/text_files/for_grep_$prefix.bed
##grab sequences at indel BED file coordinates
twoBitToFa  $path_to_ref_genome_2bit  -bed=$job_dir/text_files/$prefix.bed  $job_dir/genome_files/$prefix.fa
##create custom indel genome files but inserting indel sequences into fasta file from previous step, where appropriate
$path_to_helper_scripts/create_custom_fasta_with_indel_081516.pl   $job_dir/genome_files/$prefix.fa $step_size $job_dir/text_files/for_grep_$prefix.bed  $job_dir/genome_files/indel_genome_$prefix.fa
##build bowtie index of new indel genome file
bowtie2-build $job_dir/genome_files/indel_genome_$prefix.fa  $job_dir/genome_files/indel_genome_$prefix

echo "Aligning to custom indel genome"
bowtie2 -x $job_dir/genome_files/indel_genome_$prefix -U $job_dir/fastq/for_realignment_$s.fastq -S $job_dir/BAM/indel_genome_$prefix.sam

echo "Getting IDs of reads mapping perfectly to the indel genome. Those reads will be excluded from indel calling analyses because they map perfectly to the reference genome"
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_genome_$prefix.sam | grep -E  '.+AS:i:0.+' | cut -f 1   > $job_dir/text_files/$prefix"_ids_of_reads_to_remove.txt"

echo "Getting IDs of reads that map perfectly to reference genome"
awk '$6 ~ /^[0-9]+M$/' $job_dir/BAM/indel_realigned.$s.sam | grep -E  '.+AS:i:0.+' > $job_dir/BAM/indel_realigned_perfect_matches.$s.sam

echo "Removing ambiguously mapping reads (those aligning to both reference and indel genome)"
grep -F -vf $job_dir/text_files/$prefix"_ids_of_reads_to_remove.txt"  $job_dir/BAM/indel_realigned_perfect_matches.$s.sam > $job_dir/BAM/indel_realigned_filtered.$s.sam

echo "Combining indel_realigned/filtered BAM and reads that do not perfectly align to reference genome (potential indel reads)"
samtools view -b -o $job_dir/BAM/for_mpileup_part2_indel_realigned.$s.bam -T $path_to_ref_genome $job_dir/BAM/indel_realigned_filtered.$s.sam
samtools view -b -o $job_dir/BAM/for_mpileup_part1_indel_realigned.$s.bam -T $path_to_ref_genome $job_dir/BAM/for_mpileup_part1_indel_realigned.$s.sam
samtools merge -f  $job_dir/BAM/for_mpileup.indel_realigned_and_filtered.$s.bam $job_dir/BAM/for_mpileup_part1_indel_realigned.$s.bam $job_dir/BAM/for_mpileup_part2_indel_realigned.$s.bam

echo "Get mpileup output for merged final BAM" 
samtools mpileup -f $path_to_ref_genome   $job_dir/BAM/for_mpileup.indel_realigned_and_filtered.$s.bam > $job_dir/text_files/$s"_mpileup"
echo "Final step: run Varscan mpileup2indel on BAM"
java -Xmx4g -Djava.io.tmpdir=$java_tmp_dir_for_varscan -jar $path_to_varscan mpileup2indel $job_dir/text_files/$s"_mpileup" --min-var-freq 0.2 --min-coverage $cov> $job_dir/text_files/$s"_mpileup2indel_varscan"

echo "Pipeline finished. The main output files of interest are: $job_dir/text_files/"$s"_mpileup2indel_varscan and  $job_dir/BAM/for_mpileup.indel_realigned_and_filtered."$s".bam.  All other intermediate files were kept for troubleshooting purposes only."


