
# Genome Seq pipeline
# piper 
# https://github.com/bowhan/piper.git
# An integrated pipeline for piRNA and transposon analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome/RACE, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Weng lab
# Howard Hughes Medical Institute
# RNA Therapeutics Institute
# University of Massachusetts Medical School

##########
# Config #
##########
export GENOMESEQ_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

Genome Seq pipeline single library mode v$GENOMESEQ_VERSION from the $BOLD$PACKAGE_NAME$RESET
$GENOMIC_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thanks for using it. 

${UNDERLINE}usage${RESET}:
	piper dna \ 
		-l left.fq \ 
		-r right.fq \ 
		-g dm3 \ 
		-d 500 [100] \ 
		-o output_directory [current directory] \ 
		-c cpu [8] 
	
OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
	-l      Left reads from Paired-End sequencing; required
	-r      Right reads from Paired-End sequencing; required
	-g      Genome assembly name, like mm9 or dm3. required 
		 Check $PIPELINE_DIRECTORY/common/genome_supported.txt for genome assemblies currently installed; 
		 Use "install" to install new genome
${OPTIONAL}[ optional ]
	-d      VCF filtering depth, passed to "vcfutils.pl varFilter -D" and "retroseq.pl -call -depth", default: 100
	-o      Output directory, default: current directory $PWD
	-c      Number of CPUs to use, default: 8
	
EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hl:r:c:o:g:d:v" OPTION; do
	case $OPTION in
		h)	usage && exit 1 ;;
		l)	LEFT_FASTQ=`readlink -f $OPTARG` ;;
		r)	RIGHT_FASTQ=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		v)	echo2 "GENOMESEQ_VERSION: v$GENOMESEQ_VERSION" && exit 0 ;;		
		d)	VCFFILTER_DEPTH=$OPTARG ;;
		*)	usage && exit 1 ;;
	esac
done
# if INPUT_FASTQ or GENOME is undefined, print out usage and exit
[[ -z $LEFT_FASTQ ]] && usage && echo2 "Missing option -l for input fastq, left file " "error" 
[[ -z $RIGHT_FASTQ ]] && usage && echo2 "Missing option -r for input fastq, right file " "error" 
[[ -z $GENOME ]]  && usage && echo2 "Missing option -g for specifying which genome assembly to use" "error" 
# check whether the this genome is supported or not
check_genome $GENOME
[ ! -f $LEFT_FASTQ ] && echo2 "Cannot find input file $LEFT_FASTQ" "error"
[ ! -f $RIGHT_FASTQ ] && echo2 "Cannot find input file $RIGHT_FASTQ" "error"
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z $VCFFILTER_DEPTH ] || VCFFILTER_DEPTH=100 # 
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR

BOWTIE2_GENOMIC_MAPPING_DIR=bowtie2_output && mkdir -p $BOWTIE2_GENOMIC_MAPPING_DIR
BWA_GENOMIC_MAPPING_DIR=bwa_bcftools_output && mkdir -p $BWA_GENOMIC_MAPPING_DIR
MRFAST_GENOMIC_MAPPING_DIR=mrFast_VariationHunter_output && mkdir -p $MRFAST_GENOMIC_MAPPING_DIR
BREAKDANCER_DIR=break_dancer_out && mkdir -p $BREAKDANCER_DIR
RETROSEQ=retroSeq_discovering && mkdir -p $RETROSEQ
SUMMARY_DIR=summaries && mkdir -p $SUMMARY_DIR
BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR

########################
# running binary check #
########################
checkBin "md5sum"
checkBin "awk"
checkBin "perl"
checkBin "python"
checkBin "samtools"
checkBin "gs"
checkBin "Rscript"
checkBin "bowtie2"
checkBin "ParaFly"
checkBin "bedtools_piper"
checkBin "bedGraphToBigWig"
checkBin "express"
checkBin "mrfast"
checkBin "multiInd_SetCover"
checkBin "VH"

#############
# Variables #
#############
STEP=1
JOBUID=`echo ${LEFT_FASTQ} | md5sum | cut -d" " -f1`
LEFT_FASTQ_NAME=`basename $LEFT_FASTQ`
RIGHT_FASTQ_NAME=`basename $RIGHT_FASTQ`
PREFIX=`echo -e "${LEFT_FASTQ_NAME}\n${RIGHT_FASTQ_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${LEFT_FASTQ_NAME%.f[aq]} # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
# read lengtg
READ_LEN=`head -2 $LEFT_FASTQ | awk '{getline; printf "%d", length($1)}'`
# directories storing the common files for this organism
export COMMON_FOLDER=$PIPELINE_DIRECTORY/common/$GENOME
# assign different values to the generalized variables (same name for different GENOMEs) according to which GENOME fed
. $COMMON_FOLDER/variables
# fasta file for the genome
export GENOME_FA=$COMMON_FOLDER/${GENOME}.fa
# chrom information of this GENOME
CHROM=$COMMON_FOLDER/${GENOME}.ChromInfo.txt
# bowtie2 index
export BOWTIE2_INDEXES=$COMMON_FOLDER/Bowtie2Index
# bwa index
export BWA_INDEXES=$COMMON_FOLDER/BWAIndex
# mrFast indx
export MRFAST_INDEX=$COMMON_FOLDER/mrFastIndex

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] Genome-Seq pipeline version $GENOMESEQ_VERSION" 

###########################
# determine fastQ version #
###########################
echo2 "Determining the version of fastQ using SolexaQA"
# determine version of fastq used, using a modified SolexaQA.pl
PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piper.pl ${LEFT_FASTQ}`
case ${PHRED_SCORE} in
solexa)		bowtie2PhredOption="--solexa-quals" ;; # Solexa+64, raw reads typically (-5, 40)
illumina)	bowtie2PhredOption="--phred64" ;; # Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
sanger)		bowtie2PhredOption="--phred33" ;; # Phred+33,  raw reads typically (0, 40) (http://en.wikipedia.org/wiki/FASTQ_format)
*)			echo2 "unable to determine the fastq version. Using sanger..." "warning";;
esac

#####################################
# Align reads to genome: 1. Bowtie2 #
#####################################
# so far, no tool uses bowtie2 output
# echo2 "Mapping to genome ${GENOME} with Bowtie2"
# [ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2 ] && \
# bowtie2 -x genome \
# 	-1 ${LEFT_FASTQ} \
# 	-2 ${RIGHT_FASTQ} \
# 	-q \
# 	$bowtie2PhredOption \
# 	--very-sensitive-local \
# 	-X 800 \
# 	--no-mixed \
# 	-p $CPU \
# 	2> ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.log | \
# 	samtools view -uS -f0x2 - > ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.bam && \
# 	samtools sort -@ $CPU ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.bam ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.sorted && \
# 	samtools index ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.sorted.bam && \
# 	rm -rf ${BOWTIE2_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.b2.bam && \
# 	touch .${JOBUID}.status.${STEP}.genome_mapping_bowtie2
# STEP=$((STEP+1))

#################################
# Align reads to genome: 2. BWA #
#################################
echo2 "Mapping to genome ${GENOME} with BWA and calling varation by bcftools"
[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bwa ] && \
	bwa mem \
		-t $CPU \
		-c 1000000 \
		$BWA_INDEXES/genome \
		$LEFT_FASTQ \
		$RIGHT_FASTQ \
		2> $BWA_GENOMIC_MAPPING_DIR/${PREFIX}.${GENOME}.bwa.log | \
	samtools view -uS -f0x2 - > ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.bam && \
	samtools sort -@ $CPU ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.bam ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted && \
	samtools index ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted.bam && \
	samtools mpileup -uf $GENOME_FA ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted.bam | \
	bcftools view -bvcg - | tee ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.var.raw.bcf | \
	bcftools view - | \
	perl $PIPELINE_DIRECTORY/bin/vcfutils.pl varFilter -D $VCFFILTER_DEPTH > ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.var.flt.vcf && \
	rm -rf ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.bam && \
	touch .${JOBUID}.status.${STEP}.genome_mapping_bwa
STEP=$((STEP+1))

#################################
# Discovering SV by BreakDancer #
#################################
echo2 "Discovering Variation by BreakDancer using BWA output" 
[ ! -f .${JOBUID}.status.${STEP}.BreakDancer ] && \
	perl $PIPELINE_DIRECTORY/bin/bam2cfg_piper.pl \
		${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted.bam \
		> $BREAKDANCER_DIR/config && \
	breakdancer-max \
		-d $BREAKDANCER_DIR/${PREFIX}.breakdancer \
		-g $BREAKDANCER_DIR/${PREFIX}.breakdancer.SV.bed \
		$BREAKDANCER_DIR/config \
		2> $BREAKDANCER_DIR/${PREFIX}.breakdancer.log && \
touch .${JOBUID}.status.${STEP}.BreakDancer
STEP=$((STEP+1))

####################################
# Discovering deletion by retroSeq #
####################################
echo2 "Discovering deletions by retroSeq using BWA output"
[ ! -f .${JOBUID}.status.${STEP}.retroSeq ] && \
	echo $COMMON_FOLDER/UCSC.RepeatMask.bed > $RETROSEQ/exd.file && \
	perl $PIPELINE_DIRECTORY/bin/retroseq.pl -discover \
		-bam  ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted.bam \
		-eref $COMMON_FOLDER/${GENOME}.repBase.eref \
		-exd  $RETROSEQ/exd.file \
		-output $RETROSEQ/${PREFIX}.${GENOME}.bwa.sorted.retroSeq \
		-align \
		1>&2 2> $RETROSEQ/${PREFIX}.retroSeq.discover.log && \
	perl $PIPELINE_DIRECTORY/bin/retroseq.pl -call \
		-ref  $GENOME_FA \
		-bam  ${BWA_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.bwa.sorted.bam \
		-input $RETROSEQ/${PREFIX}.${GENOME}.bwa.sorted.retroSeq* \
		-output $RETROSEQ/${PREFIX}.${GENOME}.bwa.sorted.bam.vcf \
		-filter $COMMON_FOLDER/UCSC.RepeatMask.bed \
		-reads 10 \
		-depth $VCFFILTER_DEPTH \
		1>&2 2> $RETROSEQ/${PREFIX}.retroSeq.call.log && \
	touch .${JOBUID}.status.${STEP}.retroSeq
STEP=$((STEP+1))

#####################################
# Align reads to genome: 3. mrFast #
#####################################
echo2 "Mapping to genome ${GENOME} with mrFast. The memory used by mrFast is propotional to the size of input, so it might uses lots of memory. Please be ware of this."
mrFast_min=0
mrFast_max=800
[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_mrFast ] && \
	mrfast --search $MRFAST_INDEX/${GENOME}.fa --pe --discordant-vh --seq1 $LEFT_FASTQ --seq2 $RIGHT_FASTQ --min $mrFast_min --max $mrFast_max -o ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.sam 1>&2 2> ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.log && \
	samtools view -uS -f0x2 ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.sam > ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.bam && \
	samtools sort -@ $CPU ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.bam ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.sorted && \
	samtools index ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.sorted.bam && \
	rm -rf ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.sam ${MRFAST_GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.mrFast.bam unmapped && \
	touch .${JOBUID}.status.${STEP}.genome_mapping_mrFast
STEP=$((STEP+1))

########################################
# Variation calling by VariationHunter #
########################################
echo2 "Calling variation by VariationHunter"
[ ! -f .${JOBUID}.status.${STEP}.VariationHunter -a -f .${JOBUID}.status.${STEP}.genome_mapping_mrFast ] && \
	echo 1 > $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.sample.lib && \
	echo -e "${PREFIX}\t${PREFIX}\t${mrFast_min}\t${mrFast_max}\t${READ_LEN}" > $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.sample.lib && \
	VH -c $CHROM -i $PIPELINE_DIRECTORY/bin/VH_initInfo -l $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.sample.lib -r $COMMON_FOLDER/UCSC.RepeatMask.Satellite.bed -g $COMMON_FOLDER/${GENOME}.gap.bed -o $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.out -t $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.sample.name -x 500 -p 0.001 1>&2 2> $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.log       && \
	setCover -l $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.sample.lib -r $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.sample.name -c $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.out -t 1000000 -o $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.VHcluster.sample.Out.SV 1>&2 2> $MRFAST_GENOMIC_MAPPING_DIR/${PREFIX}.setCover.log && \
	touch .${JOBUID}.status.${STEP}.VariationHunter
STEP=$((STEP+1))

#############
# finishing #
#############
echo2 "Finished running ${PACKAGE_NAME} Genome-Seq pipeline version $GENOMESEQ_VERSION"
echo2 "---------------------------------------------------------------------------------"
touch .${GENOME}.GENOMESEQ_VERSION.${GENOMESEQ_VERSION}




