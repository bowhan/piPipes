
# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) 2014  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##########
# Config #
##########
export CHIPSEQ_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

ChIP Seq pipeline single library mode v$CHIPSEQ_VERSION from the $BOLD$PACKAGE_NAME$RESET
$CHIP_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs.
Thank you for using it.

==================< paired-end >==================
	${UNDERLINE}usage${RESET}:
		pipile	chip \
			-l left.IP.fq \
			-r right.IP.fq \
			-L left.INPUT.fq \
			-R right.INPUT.fq \
			-g mm9 \
			-c cpu[8] \
			-o output_dir[cwd] \
			-x 500 [1000] \
			-M /path/to/user/defined/region.bed1,/path/to/user/defined/region2.bed,/path/to/user/defined/region3.bed...

==================< single-end >==================
	${UNDERLINE}usage${RESET}:
		pipile	chip \
			-i IP.fq \
			-I INPUT.fq \
			-g mm9 \
			-c cpu[8] \
			-o output_dir[cwd] \
			-x 500 [1000] \
			-M /path/to/user/defined/region.bed1,/path/to/user/defined/region2.bed,/path/to/user/defined/region3.bed...

OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
========================< paired-end >========================
	-l      Left reads from Paired-End sequencing of ChIP-Seq Immunoprecipitation
	-r      Right reads from Paired-End sequencing of ChIP-Seq Immunoprecipitation
	-L      Left reads from Paired-End sequencing of ChIP-Seq Input
	-R      Right reads from Paired-End sequencing of ChIP-Seq Input
========================< single-end >========================
	-i      Reads from Single-End sequencing of ChIP-Seq Immunoprecipitation
	-I      Reads from Single-End sequencing of ChIP-Seq Input

	-g      Genome assembly name, like mm9 or dm3. required
		 Check $PIPELINE_DIRECTORY/common/genome_supported.txt for genome assemblies currently installed;
		 Use "piPipes install" to install new genome
${OPTIONAL}[ optional ]
	-B      Use "Broader" peak calling algorithm in MACS2; should be used for library like H3K9me3. default: off
	-Q      MAPQ value used as threshold to filter alignments from Bowtie2 SAM/BAM output. The higher the MAPQ, the more unique the alignment. default: 1
		 NOTE: Bowtie2 is ALWAYS called without -k or -a. Consequently, for each read, there will only be one best alignment given, no matter how many times it maps.
		 If all the alignments are equally good, Bowtie2 chooses one randomly (see Bowtie2 manual for more details).
		 Including multiple mappers ensures repetitive regions are not depleted. Keeping only one alignment avoids false positive signal on repetitive regions.
		*If you would like to only consider unique mapper, set this option to more than 1, like 10. This usually removes multiple mappers which have equal good alignments.
	-x      Length to extend up/downstream of each genomic features to draw the metagene plot. default: 1000
	-M      Path to BED files for meta-plot analysis on user-defined regions. One plot for each BED file and different BED files should be delimited by comma.
	     Do not use ~ to represent the home directory as ~ is expanded before passed into piPipes and it will not be expanded correctly unless repsent at the beginning of a string.
		 Use $HOME instead, like -M $HOME/path/to/gfp,$HOME/path/to/white.fa
	-o      Output directory, default: current directory: $PWD
	-c      Number of CPUs to use, default: 8
EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hl:r:L:R:c:o:g:Bvx:Q:i:I:M:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		v)	echo2 "CHIPSEQ_VERSION: v$CHIPSEQ_VERSION" && exit 0 ;;
		l)	LEFT_IP_FASTQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
		r)	RIGHT_IP_FASTQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
		L)	LEFT_INPUT_FASTQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
		R)	RIGHT_INPUT_FASTQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
		i)	IP_FASTQ=`readlink -f $OPTARG`; SE_MODE=1 ;;
		I)	INPUT_FASTQ=`readlink -f $OPTARG`; SE_MODE=1 ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		M)	export USER_DEFINED_BED_FILES=$OPTARG ;;
		c)	export CPU=$OPTARG ;;
		g)	export GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		B)	export MACS2_BROAD_OPT="--broad" ;;
		Q)	export MINIMAL_MAPQ=$OPTARG ;;
		x)	export EXT_LEN=$OPTARG ;;
		*)	usage && exit 1 ;;
	esac
done

if [[ -z $PE_MODE && -z $SE_MODE ]]; then usage ; echo2 "Please specify the input file!" "error"; fi
if [[ -n $PE_MODE && -n $SE_MODE ]]; then usage ; echo2 "Please only choose single-end OR paired-end, but not both" "error"; fi

if [[ -n $PE_MODE ]]; then
	[[ -z $LEFT_IP_FASTQ ]] && usage && echo2 "Missing option -l for IP fastq of left file, or file does not exist " "error"
	[[ -z $RIGHT_IP_FASTQ ]] && usage && echo2 "Missing option -r for IP fastq of right file, or file does not exist " "error"
	[[ -z $LEFT_INPUT_FASTQ ]] && usage && echo2 "Missing option -L for INPUT fastq of left file, or file does not exist " "error"
	[[ -z $RIGHT_INPUT_FASTQ ]] && usage && echo2 "Missing option -R for INPUT fastq of right file, or file does not exist " "error"
	[[ ! -f $LEFT_IP_FASTQ ]] && usage && echo2 "Missing option -l for IP fastq of left file, or file does not exist " "error"
	[[ ! -f $RIGHT_IP_FASTQ ]] && usage && echo2 "Missing option -r for IP fastq of right file, or file does not exist " "error"
	[[ ! -f $LEFT_INPUT_FASTQ ]] && usage && echo2 "Missing option -L for INPUT fastq of left file, or file does not exist " "error"
	[[ ! -f $RIGHT_INPUT_FASTQ ]] && usage && echo2 "Missing option -R for INPUT fastq of right file, or file does not exist " "error"
fi

if [[ -n $SE_MODE ]]; then
	[[ -z $IP_FASTQ ]] && usage && echo2 "Missing option -i for IP fastq, or file does not exist " "error"
	[[ -z $INPUT_FASTQ ]] && usage && echo2 "Missing option -I for IP fastq, or file does not exist " "error"
	[[ ! -f $IP_FASTQ ]] && usage && echo2 "Missing option -i for IP fastq, or file does not exist " "error"
	[[ ! -f $INPUT_FASTQ ]] && usage && echo2 "Missing option -I for IP fastq, or file does not exist " "error"
fi

[[ -z $GENOME ]] && usage && echo2 "Missing option -g for specifying which genome assembly to use" "error"
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")

# check whether the this genome is supported or not
check_genome $GENOME
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8; export CPU
[ ! -z $MINIMAL_MAPQ ] || MINIMAL_MAPQ=1; export MINIMAL_MAPQ
[ ! -z "${EXT_LEN##*[!0-9]*}" ] || EXT_LEN=1000; export EXT_LEN
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
# READS_DIR=input_read_files && mkdir -p $READS_DIR
GENOMIC_MAPPING_DIR=genome_mapping && mkdir -p $GENOMIC_MAPPING_DIR
PEAKS_CALLING_DIR=macs2_peaks_calling && mkdir -p $PEAKS_CALLING_DIR
# SUMMARY_DIR=summaries && mkdir -p $SUMMARY_DIR
BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR
AGG_DIR=aggregate_output && mkdir -p $AGG_DIR

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
checkBin "bedtools_piPipes"
checkBin "bedGraphToBigWig"
checkBin "express"
checkBin "macs2"

#############
# Variables #
#############
STEP=1
if [[ -n $SE_MODE ]]; then
	JOBUID=`echo $IP_FASTQ | md5sum | cut -d" " -f1`
	IP_FASTQ_NAME=`basename $IP_FASTQ`
	INPUT_FASTQ_NAME=`basename $INPUT_FASTQ`
	export PREFIX=`echo -e "${IP_FASTQ_NAME}\n${INPUT_FASTQ_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` # && export PREFIX=${PREFIX%.*}
	[ -z "${PREFIX}" ] && export PREFIX=${IP_FASTQ_NAME%.f[aq]*}
else
	JOBUID=`echo ${LEFT_IP_FASTQ} | md5sum | cut -d" " -f1`
	LEFT_IP_FASTQ_NAME=`basename $LEFT_IP_FASTQ`
	RIGHT_IP_FASTQ_NAME=`basename $RIGHT_IP_FASTQ`
	export PREFIX=`echo -e "${LEFT_IP_FASTQ_NAME}\n${RIGHT_IP_FASTQ_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` # && export PREFIX=${PREFIX%.*}
	[ -z "${PREFIX}" ] && export PREFIX=${LEFT_FASTQ_NAME%.f[aq]*} # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
fi

# table to store the basic statistics of the library (genomic mappability)
TABLE=${PREFIX}.basic_stats
# directories storing the common files for this organism
export COMMON_FOLDER=$PIPELINE_DIRECTORY/common/$GENOME
# fasta file for the genome
export GENOME_FA=$COMMON_FOLDER/${GENOME}.fa
# chrom information of this GENOME
CHROM=$COMMON_FOLDER/${GENOME}.ChromInfo.txt
# bowtie2 index
export BOWTIE2_INDEXES=$COMMON_FOLDER/Bowtie2Index

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] ChIP-Seq pipeline version $CHIPSEQ_VERSION"

###########################
# determine fastQ version #
###########################
echo2 "Determining the version of fastQ using SolexaQA"
# determine version of fastq used, using a modified SolexaQA.pl
if [[ -n $SE_MODE ]]; then
	PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piPipes.pl $IP_FASTQ`
	READ_LEN=`head -2 $IP_FASTQ | tail -1 | awk '{print length($1)}'`
else
	PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piPipes.pl ${LEFT_IP_FASTQ}`
	READ_LEN=`head -2 $LEFT_IP_FASTQ | tail -1 | awk '{print length($1)}'`
fi
echo $READ_LEN > .READ_LEN
case ${PHRED_SCORE} in
solexa)		bowtie2PhredOption="--solexa-quals" ;; # Solexa+64, raw reads typically (-5, 40)
illumina)	bowtie2PhredOption="--phred64" ;; # Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
sanger)		bowtie2PhredOption="--phred33" ;; # Phred+33,  raw reads typically (0, 40) (http://en.wikipedia.org/wiki/FASTQ_format)
*)			echo2 "unable to determine the fastq version. Using sanger..." "warning";;
esac

############################
# Align IP reads to genome #
############################
echo2 "Mapping IP reads to genome ${GENOME} with Bowtie2"
if [[ -n $SE_MODE ]]; then
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP ] && \
	bowtie2 -x genome \
		-U $IP_FASTQ \
		-q \
		$bowtie2PhredOption \
		--very-sensitive-local \
		-p $CPU \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.log | \
		samtools view -uS -F0x4 -q ${MINIMAL_MAPQ} - > ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted.bam && \
		rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam && \
		touch .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP ] && echo2 "Failed in mapping IP to genome" "error"
else
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP ] && \
	bowtie2 -x genome \
		-1 ${LEFT_IP_FASTQ} \
		-2 ${RIGHT_IP_FASTQ} \
		-q \
		$bowtie2PhredOption \
		--very-sensitive-local \
		-X 800 \
		--no-mixed \
		-p $CPU \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.log | \
		samtools view -uS -f 0x2 -q ${MINIMAL_MAPQ} - > ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted.bam && \
		rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.bam && \
		touch .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2IP ] && echo2 "Failed in mapping IP to genome" "error"
fi
STEP=$((STEP+1))

###############################
# Align Input reads to genome #
###############################
echo2 "Mapping Input reads to genome ${GENOME} with Bowtie2"
if [[ -n $SE_MODE ]]; then
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_Input ] && \
	bowtie2 -x genome \
		-U $INPUT_FASTQ \
		-q \
		$bowtie2PhredOption \
		--very-sensitive-local \
		-p $CPU \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.log | \
		samtools view -uS -F0x4 -q ${MINIMAL_MAPQ} - > ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted.bam && \
		rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam && \
		touch .${JOBUID}.status.${STEP}.genome_mapping_bowtie2Input
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_bowtie2Input ] && echo2 "Failed in mapping Input to genome" "error"
else
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_Input ] && \
	bowtie2 -x genome \
		-1 ${LEFT_INPUT_FASTQ} \
		-2 ${RIGHT_INPUT_FASTQ} \
		-q \
		$bowtie2PhredOption \
		-X 800 \
		--very-sensitive-local \
		--no-mixed \
		-p $CPU \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.log | \
		samtools view -uS -f 0x2 -q ${MINIMAL_MAPQ} - > ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted.bam && \
		rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.bam && \
		touch .${JOBUID}.status.${STEP}.genome_mapping_Input
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping_Input ] && echo2 "Failed in mapping input to genome" "error"
fi
STEP=$((STEP+1))

#######################################
# Call peaks using MACS2, with --SPMR #
#######################################
echo2 "Calling peaks with MACS2 and make enrichment bedGraph and bigWig"
GENOME_SIZE=`awk '{a+=$2}END{print a}' $CHROM`
if [[ -n $SE_MODE ]]; then
	[ ! -f .${JOBUID}.status.${STEP}.peak_calling_with_macs2 ] && \
		macs2 callpeak \
			-f BAM \
			-t ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted.bam \
			-c ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted.bam \
			-g $GENOME_SIZE \
			$MACS2_BROAD_OPT \
			--outdir $PEAKS_CALLING_DIR \
			-n ${PREFIX} \
			-B --SPMR \
			2> $PEAKS_CALLING_DIR/${PREFIX}.callpeak.SPMR.log && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.bdg -m ppois 1> ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.stderr && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.bdg $CHROM $BW_OUTDIR/${PREFIX}.ppois.bigWig && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.FE.bdg -m FE 1> ${PEAKS_CALLING_DIR}/${PREFIX}.FE.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.FE.stderr  && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.FE.bdg $CHROM $BW_OUTDIR/${PREFIX}.FE.bigWig && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.bdg -m logLR -p 0.00001 1> ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.stderr  && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.bdg $CHROM $BW_OUTDIR/${PREFIX}.logLR.bigWig && \
	touch  .${JOBUID}.status.${STEP}.peak_calling_with_macs2
else
	[ ! -f .${JOBUID}.status.${STEP}.peak_calling_with_macs2 ] && \
		macs2 callpeak \
			-f BAMPE \
			-t ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.IP.b2.sorted.bam \
			-c ${GENOMIC_MAPPING_DIR}/${PREFIX}.${GENOME}.Input.b2.sorted.bam \
			-g $GENOME_SIZE \
			$MACS2_BROAD_OPT \
			--outdir $PEAKS_CALLING_DIR \
			-n ${PREFIX} \
			-B --SPMR \
			2> $PEAKS_CALLING_DIR/${PREFIX}.callpeak.SPMR.log && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.bdg -m ppois 1> ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.stderr && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.ppois.bdg $CHROM $BW_OUTDIR/${PREFIX}.ppois.bigWig && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.FE.bdg -m FE 1> ${PEAKS_CALLING_DIR}/${PREFIX}.FE.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.FE.stderr && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.FE.bdg $CHROM $BW_OUTDIR/${PREFIX}.FE.bigWig && \
		macs2 bdgcmp -t $PEAKS_CALLING_DIR/${PREFIX}_treat_pileup.bdg -c $PEAKS_CALLING_DIR/${PREFIX}_control_lambda.bdg -o ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.bdg -m logLR -p 0.00001 1> ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.stdout 2> ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.stderr && \
		bedGraphToBigWig ${PEAKS_CALLING_DIR}/${PREFIX}.logLR.bdg $CHROM $BW_OUTDIR/${PREFIX}.logLR.bigWig && \
	touch  .${JOBUID}.status.${STEP}.peak_calling_with_macs2
fi
STEP=$((STEP+1))

############################################
# draw figures for genomic features (mega) #
############################################
echo2 "Aggregating signal on each genomic features"
[ ! -f .${JOBUID}.status.${STEP}.aggregate_beds ] && \
	bash $DEBUG piPipes_aggregate_bw_on_beds.sh \
	$AGG_DIR \
	$EXT_LEN \
	$BW_OUTDIR/${PREFIX}.ppois.bigWig,$BW_OUTDIR/${PREFIX}.FE.bigWig,$BW_OUTDIR/${PREFIX}.logLR.bigWig && \
	touch .${JOBUID}.status.${STEP}.aggregate_beds
STEP=$((STEP+1))

#############
# finishing #
#############
echo2 "Finished running ${PACKAGE_NAME} ChIP-Seq pipeline version $CHIPSEQ_VERSION"
echo2 "---------------------------------------------------------------------------------"
echo $MACS2_BROAD_OPT > .MACS2_BROAD_OPT # only output this option when the run finishs
if [[ -n $SE_MODE ]]; then
	echo "SE" > .SE_OR_PE
else
	echo "PE" > .SE_OR_PE
fi
touch .${GENOME}.CHIPSEQ_VERSION.${CHIPSEQ_VERSION}
