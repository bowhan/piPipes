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
export CHIPSEQ2_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

ChIP Seq pipeline single library mode v$CHIPSEQ2_VERSION from the $BOLD$PACKAGE_NAME$RESET
$CHIP2_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs.
Thank you for using it.

${UNDERLINE}usage${RESET}:
	piPipes chip2	\ 
		-a chipseq_pipiline_output_dir1 \ 
		-b chipseq_pipiline_output_dir2 \ 
		-g dm3 \ 
		-c 24 [8] \ 
		-o output_dir [cwd] \ 
		-x 2000 [1000] \ 
		-A piwi_heterozygous [basename of -a] \ 
		-B piwi_mutant [basename of -b]

OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
	-a      Directory to the folder with the output of single library mode, for sample A (wild-type)
	-b      Directory to the folder with the output of single library mode, for sample B (mutant)
	-g      Genome assembly name, like mm9 or dm3
	        Check "$PIPELINE_DIRECTORY/common/genome_supported.txt" for genome assemblies currently installed;
	        Use "install" to install new genome
${OPTIONAL}[ optional ]
	-x      Length (bp) to extend up/downstream of each genomic features to draw the metagene plot. default: 1000
	-c      Number of CPUs to use, default: 8
	-o      Output directory, default: current directory $PWD
	-A      Name to use for Sample A, default: using the basename of -a
	-B      Name to use for Sample B, default: using the basename of -b

The pipeline will automatically detect the version and options of the single library run for the two samples and ensure the consistence.

EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hva:b:g:c:o:A:B:x:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		v)	echo2 "SMALLRNA2_VERSION: v$SMALLRNA2_VERSION" && exit 0 ;;
		a)	SAMPLE_A_DIR=`readlink -f $OPTARG` ;;
		b)	SAMPLE_B_DIR=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=${OPTARG};;
		A)  export SAMPLE_A_NAME=$OPTARG ;;
		B)  export SAMPLE_B_NAME=$OPTARG ;;
		x)	EXT_LEN=$OPTARG ;;
		*)	usage && exit 1 ;;
	esac
done
# if INPUT_FASTQ or GENOME is undefined, print out usage and exit
[[ -z $SAMPLE_A_DIR ]] && usage && echo2 "Missing option -a for input dir for sample A (wild-type) file " "error"
[[ ! -d $SAMPLE_A_DIR ]] && echo2 "Cannot find input directory $SAMPLE_A_DIR" "error"
[[ -z $SAMPLE_B_DIR ]] && usage && echo2 "Missing option -b for input dir for sample B (mutant) file " "error"
[[ ! -d $SAMPLE_B_DIR ]] && echo2 "Cannot find input directory $SAMPLE_B_DIR" "error"
[[ -z $GENOME ]]  && usage && echo2 "Missing option -g for specifying which genome assembly to use" "error"
check_genome $GENOME
[[ -z $SAMPLE_A_NAME ]] && export SAMPLE_A_NAME=`basename $SAMPLE_A_DIR`
[[ -z $SAMPLE_B_NAME ]] && export SAMPLE_B_NAME=`basename $SAMPLE_B_DIR`
PREFIX=`echo -e "${SAMPLE_A_NAME}\n${SAMPLE_B_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${SAMPLE_B_NAME}_${SAMPLE_B_NAME}
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "${EXT_LEN##*[!0-9]*}" ] || EXT_LEN=1000
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################
# Version Check #
#################
SAMPLE_A_VERSION=`ls -a $SAMPLE_A_DIR | grep CHIPSEQ_VERSION`
SAMPLE_B_VERSION=`ls -a $SAMPLE_B_DIR | grep CHIPSEQ_VERSION`
[ -z $SAMPLE_A_VERSION ] && echo2 "Please make sure the pipeline finishes for $SAMPLE_A_DIR" "error"
[ -z $SAMPLE_B_VERSION ] && echo2 "Please make sure the pipeline finishes for $SAMPLE_B_DIR" "error"
[ "$SAMPLE_A_VERSION" != "$SAMPLE_B_VERSION" ] && echo2 "It appears that the two runs were not done by the same assemly or same version." "error"
BROAD_A=`cat $SAMPLE_A_DIR/.MACS2_BROAD_OPT`
BROAD_B=`cat $SAMPLE_B_DIR/.MACS2_BROAD_OPT`
[ "$BROAD_A" != "$BROAD_B" ] && echo2 "The two run were not done by the same -b option" "error"
SE_OR_PE_A=`cat $SAMPLE_A_DIR/.SE_OR_PE`
SE_OR_PE_B=`cat $SAMPLE_B_DIR/.SE_OR_PE`
[ "$SE_OR_PE_A" != "$SE_OR_PE_B" ] && echo2 "The two run were not done by the type of SE/PE" "error"
if [[ "${SE_OR_PE_A}" == "SE" ]]; then MACS2_f="BAM"; TN="tags"; else MACS2_f="BAMPE"; TN="fragments"; fi

########################
# running binary check #
########################
checkBin "gs"
checkBin "Rscript"
checkBin "macs2"

#################################
# creating output files/folders #
#################################
TABLE=${PREFIX}.basic_stats
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
PEAKS_CALLING_DIR_A=macs2_peaks_calling_no_normalization_${SAMPLE_A_NAME} && mkdir -p $PEAKS_CALLING_DIR_A
PEAKS_CALLING_DIR_B=macs2_peaks_calling_no_normalization_${SAMPLE_B_NAME} && mkdir -p $PEAKS_CALLING_DIR_B
BDGDIFF_DIR=differential_peaks_calling && mkdir -p $BDGDIFF_DIR
# BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR
AGG_DIR=aggregate_output && mkdir -p $AGG_DIR
DIRECT_MAPPING=direct_mapping && mkdir -p $DIRECT_MAPPING

#############
# Variables #
#############
# step counter
STEP=1
# job uid
JOBUID=`echo ${PREFIX} | md5sum | cut -d" " -f1`
# directories storing the common files for this organism
export COMMON_FOLDER=$PIPELINE_DIRECTORY/common/$GENOME
# fasta file for the genome
export GENOME_FA=$COMMON_FOLDER/${GENOME}.fa
# chrom information of this GENOME
CHROM=$COMMON_FOLDER/${GENOME}.ChromInfo.txt
GENOME_SIZE=`awk '{a+=$2}END{print a}' $CHROM`

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] ChIPSeq pipeline dual library mode version $CHIPSEQ2_VERSION"

########################################################
# Call peaks using MACS2, without --SPMR, for sample A #
########################################################
# in single library mode, the peak calling were done with normalization (--SPMR)
# in order to call differential binding events, the two libraries need to be normalized here
# ref: https://github.com/taoliu/MACS/wiki/Call-differential-binding-events
SAMPLE_A_IP_BAM=$SAMPLE_A_DIR/genome_mapping/*IP.sorted.bam
SAMPLE_A_INPUT_BAM=$SAMPLE_A_DIR/genome_mapping/*Input.sorted.bam
echo2 "Running MACS2 without normalization for $SAMPLE_A_NAME"
[ ! -f .${JOBUID}.status.${STEP}.peak_calling_A ] && \
	macs2 callpeak \
		-f $MACS2_f \
		-t ${SAMPLE_A_IP_BAM} \
		-c ${SAMPLE_A_INPUT_BAM} \
		-g $GENOME_SIZE \
		$BROAD_A \
		--outdir $PEAKS_CALLING_DIR_A \
		-n ${SAMPLE_A_NAME} \
		-B \
		2> $SAMPLE_A_DIR/${SAMPLE_A_NAME}.callpeak.log && \
	touch .${JOBUID}.status.${STEP}.peak_calling_A
STEP=$((STEP+1))
EFFECTIVE_DEPTH_A=`grep "$TN after filtering in" $PEAKS_CALLING_DIR_A/*_peaks.xls | awk 'BEGIN{FS=" "; getline;m=$NF}{if (m>$NF) {m=$NF}} END{print m}'`

########################################################
# Call peaks using MACS2, without --SPMR, for sample B #
########################################################
SAMPLE_B_IP_BAM=$SAMPLE_B_DIR/genome_mapping/*IP.sorted.bam
SAMPLE_B_INPUT_BAM=$SAMPLE_B_DIR/genome_mapping/*Input.sorted.bam
echo2 "Running MACS2 without normalization for $SAMPLE_B_NAME"
[ ! -f .${JOBUID}.status.${STEP}.peak_calling_B ] && \
	macs2 callpeak \
		-f $MACS2_f \
		-t ${SAMPLE_B_IP_BAM} \
		-c ${SAMPLE_B_INPUT_BAM} \
		-g $GENOME_SIZE \
		$BROAD_B \
		--outdir $PEAKS_CALLING_DIR_B \
		-n ${SAMPLE_B_NAME} \
		-B \
		2> $SAMPLE_B_DIR/${SAMPLE_B_NAME}.callpeak.log && \
	touch .${JOBUID}.status.${STEP}.peak_calling_B
STEP=$((STEP+1))
EFFECTIVE_DEPTH_B=`grep "$TN after filtering in" $PEAKS_CALLING_DIR_B/*_peaks.xls | awk 'BEGIN{FS=" "; getline;m=$NF}{if (m>$NF) {m=$NF}} END{print m}'`

#######################################
# Call differentially expressed peaks #
#######################################
echo2 "Calling differential binding events with MACS2"
[ ! -f .${JOBUID}.status.${STEP}.bdgdiff ] && \
macs2 bdgdiff \
	--t1 $PEAKS_CALLING_DIR_A/${SAMPLE_A_NAME}_treat_pileup.bdg \
	--c1 $PEAKS_CALLING_DIR_A/${SAMPLE_A_NAME}_control_lambda.bdg \
	--t2 $PEAKS_CALLING_DIR_B/${SAMPLE_B_NAME}_treat_pileup.bdg \
	--c2 $PEAKS_CALLING_DIR_B/${SAMPLE_B_NAME}_control_lambda.bdg \
	--d1 $EFFECTIVE_DEPTH_A \
	--d2 $EFFECTIVE_DEPTH_B \
	-g `cat $SAMPLE_A_DIR/.READ_LEN` \
	--outdir $BDGDIFF_DIR \
	--o-prefix ${PREFIX}.${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME} \
	2> $BDGDIFF_DIR/${PREFIX}.${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.bdgdiff.log && \
	TO_BE_FIX=`ls $BDGDIFF_DIR/*common.bed` && awk 'BEGIN{FS=OFS="\t";getline}{$5=int($5);print $0,"+"}' $TO_BE_FIX > ${TO_BE_FIX}1 && mv ${TO_BE_FIX}1 ${TO_BE_FIX} && \
	TO_BE_FIX=`ls $BDGDIFF_DIR/*_cond1.bed` && awk 'BEGIN{FS=OFS="\t";getline}{$5=int($5);print $0,"+"}' $TO_BE_FIX > ${TO_BE_FIX}1 && mv ${TO_BE_FIX}1 ${TO_BE_FIX} && \
	TO_BE_FIX=`ls $BDGDIFF_DIR/*_cond2.bed` && awk 'BEGIN{FS=OFS="\t";getline}{$5=int($5);print $0,"+"}' $TO_BE_FIX > ${TO_BE_FIX}1 && mv ${TO_BE_FIX}1 ${TO_BE_FIX} && \
	touch .${JOBUID}.status.${STEP}.bdgdiff
STEP=$((STEP+1))

###########################################
# Annotate peaks; count transposon copies #
###########################################
echo2 "Annotating peaks and counting transposon copies"
if [ ! -f .${JOBUID}.status.${STEP}.annotate_peaks ]; then
	COMMON_PEAKS_BED=`ls $BDGDIFF_DIR/*common.bed`
	COND1_PEAKS_BED=`ls $BDGDIFF_DIR/*_cond1.bed`
	COND2_PEAKS_BED=`ls $BDGDIFF_DIR/*_cond2.bed`
	. $COMMON_FOLDER/genomic_features
	if [[ "$GENOME" == "dm3" ]]; then 
		TRANSPOSON_BED=$COMMON_FOLDER/dm3.transposon.rm.bed #TODO renew $Trn and Zamore groups...
		GENE_BED=$flyBase_Gene
	elif [[ "$GENOME" == "dm6" ]]; then
		TRANSPOSON_BED=$flyBase_repeatMasker
		GENE_BED=$Genes
	elif [[ "$GENOME" == "mm9" ]]; then
		TRANSPOSON_BED=$repeatMasker
		GENE_BED=$refSeq_GENE
	fi
	
	for i in COMMON_PEAKS_BED COND1_PEAKS_BED COND2_PEAKS_BED; do
		BED=${!i}
		bedtools_piPipes intersect -wa -u -a ${BED} -b $TRANSPOSON_BED | \
			bedtools_piPipes intersect -v -a stdin -b $MASK  | \
			bedtools_piPipes intersect -v -a stdin -b $GENE_BED > ${BED%bed}trn-no-gene.bed && \
		bedtools_piPipes intersect -wo -a ${BED%bed}trn-no-gene.bed -b $TRANSPOSON_BED | \
			sort -k1,1 -k2,2 -k3,3 -k4,4 -k6,6 | \
			bedtools_piPipes groupby -i stdin -g 1,2,3,4,6 -c 10 -o distinct | \
			awk '{OFS="\t"; print $1,$2,$3,$4,$6,$5}' > \
		${BED%bed}trn-no-gene.annotated.bed
	done
	touch .${JOBUID}.status.${STEP}.annotate_peaks
fi

############################################
# draw figures for genomic features (mega) #
############################################
echo2 "Aggregating signal on each genomic features"
[ ! -f .${JOBUID}.status.${STEP}.aggregate_beds ] && \
	bash $DEBUG piPipes_aggregate_bw_on_beds.sh \
	$AGG_DIR \
	$EXT_LEN \
	`ls $SAMPLE_A_DIR/bigWig/*ppois.bigWig`,`ls $SAMPLE_A_DIR/bigWig/*FE.bigWig`,`ls $SAMPLE_A_DIR/bigWig/*logLR.bigWig`,`ls $SAMPLE_B_DIR/bigWig/*ppois.bigWig`,`ls $SAMPLE_B_DIR/bigWig/*FE.bigWig`,`ls $SAMPLE_B_DIR/bigWig/*logLR.bigWig` \
	`ls $BDGDIFF_DIR/*cond1.bed` \
	`ls $BDGDIFF_DIR/*cond2.bed`  && \
	touch .${JOBUID}.status.${STEP}.aggregate_beds
STEP=$((STEP+1))

##############################
# draw direct mapping figure #
##############################
DIRECT_BS_A=$SAMPLE_A_DIR/direct_mapping/*sorted.unique.bigWig.summary
DIRECT_BS_B=$SAMPLE_B_DIR/direct_mapping/*sorted.unique.bigWig.summary
DEPTH_A=`cat $SAMPLE_A_DIR/.NormScale`
DEPTH_B=`cat $SAMPLE_B_DIR/.NormScale`
[ ! -f .${JOBUID}.status.${STEP}.draw_bigWig_summary ] && \
	paste $DIRECT_BS_A $DIRECT_BS_B | cut -f1,2,3,4,7,8 > $DIRECT_MAPPING/${SAMPLE_A_NAME}_${SAMPLE_B_NAME}.direct_mapping.sum && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summaryN.R $DIRECT_MAPPING/${SAMPLE_A_NAME}_${SAMPLE_B_NAME}.direct_mapping.sum $DIRECT_MAPPING/${SAMPLE_A_NAME}_${SAMPLE_B_NAME}.direct_mapping.sum $CPU 1,1,1,1 ${SAMPLE_A_NAME} ${SAMPLE_B_NAME} 1>&2 && \
	PDFs=$DIRECT_MAPPING/*pdf && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${SAMPLE_A_NAME}_${SAMPLE_B_NAME}.direct_mapping.unique.pdf ${PDFs} && \
	rm -rf ${PDFs} && \
touch .${JOBUID}.status.${STEP}.draw_bigWig_summary

