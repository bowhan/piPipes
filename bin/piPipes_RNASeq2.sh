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
export RNASEQ2_VERSION=1.0.1

#########
# USAGE #
#########
usage () {
cat << EOF

RNASeq pipeline dual library mode v$RNASEQ2_VERSION from the $BOLD$PACKAGE_NAME$RESET
$RNASEQ2_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

${UNDERLINE}usage${RESET}:
	piPipes rna2	\ 
		-a rnaseq_pipiline_output_sample_A_rep1,rnaseq_pipiline_output_sample_A_rep2,rnaseq_pipiline_output_sample_A_rep3 \ 
		-b rnaseq_pipiline_output_sample_B_rep1,rnaseq_pipiline_output_sample_B_rep2,rnaseq_pipiline_output_sample_B_rep3 \ 
		-g dm3 \ 
		-c 24 [8] \ 
		-o output_dir [cwd] \ 
		-A piwi_heterozygous [basename of -a] \ 
		-B piwi_mutant [basename of -b] \ 
		-n 100 [50]
		
OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
	-a      Directory to the folder with the output of single library mode, for sample A (wild-type)
	-b      Directory to the folder with the output of single library mode, for sample B (mutant)
	-g      Genome assembly name, like mm9 or dm3. 
	        Check "$PIPELINE_DIRECTORY/common/genome_supported.txt" for genome assemblies currently installed; 
	        Use "install" to install new genome
${OPTIONAL}[ optional ]
	-c      Number of CPUs to use, default: 8
	-o      Output directory, default: current directory $PWD
	-A      Name to use for Sample A, default: using the basename of -a
	-B      Name to use for Sample B, default: using the basename of -b
	-n      The top n genes to draw heatmap and bar plot in cummeRbund, sorted by q-value and fold-change, default: 50

The pipeline will automatically detect the version and options of the single library run for the two samples and ensure the consistence. 

EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hva:b:g:c:o:A:B:n:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		v)	echo2 "SMALLRNA2_VERSION: v$SMALLRNA2_VERSION" && exit 0 ;;
		a)	SAMPLE_As=$OPTARG ;;
		b)	SAMPLE_Bs=$OPTARG ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=${OPTARG};;
		A)  SAMPLE_A_NAME=$OPTARG ;;
		B)  SAMPLE_B_NAME=$OPTARG ;;
		n)	NUM_GENE_CUMM=$OPTARG ;;
		*)	usage && exit 1 ;;
	esac
done

# if INPUT_FASTQ or GENOME is undefined, print out usage and exit
[[ -z $SAMPLE_As ]] && usage && echo2 "Missing option -a for input directories for sample A (wild-type) file, use comma to separete replicates " "error" 
[[ -z $SAMPLE_Bs ]] && usage && echo2 "Missing option -b for input directories for sample B (mutant) file, use comma to separete replicates" "error" 
declare -a SAMPLE_A_DIR SAMPLE_B_DIR
eval `echo $SAMPLE_As | awk 'BEGIN{FS=","}{printf "export SAMPLE_A_DIR_RELATIVE=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}'`
for DIR in "${SAMPLE_A_DIR_RELATIVE[@]}"; do SAMPLE_A_DIR+=(`readlink -f $DIR`); done
eval `echo $SAMPLE_Bs | awk 'BEGIN{FS=","}{printf "export SAMPLE_B_DIR_RELATIVE=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}'`
for DIR in "${SAMPLE_B_DIR_RELATIVE[@]}"; do SAMPLE_B_DIR+=(`readlink -f $DIR`); done
echo "${SAMPLE_B_DIR[@]}"
for DIR in "${SAMPLE_A_DIR[@]}" "${SAMPLE_B_DIR[@]}"; do [[ ! -d $DIR ]] && echo2 "Cannot find input directory $DIR" "error"; done
[[ -z $GENOME ]]  && usage && echo2 "Missing option -g for specifying which genome assembly to use" "error" 
check_genome $GENOME
[[ -z $SAMPLE_A_NAME ]] && SAMPLE_A_NAME=`basename $SAMPLE_A_DIR`
[[ -z $SAMPLE_B_NAME ]] && SAMPLE_B_NAME=`basename $SAMPLE_B_DIR`
PREFIX=`echo -e "${SAMPLE_A_NAME}\n${SAMPLE_B_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${SAMPLE_A_NAME}_${SAMPLE_B_NAME}
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ -z "${NUM_GENE_CUMM##*[!0-9]*}" ] && NUM_GENE_CUMM=50
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################
# Version Check #
#################
for DIR in "${SAMPLE_A_DIR[@]}" "${SAMPLE_B_DIR[@]}"; do 
	VERSION=`ls -a $DIR | grep RNASEQ_VERSION` 
	[[ -z $VERSION ]] && echo2 "RNASeq single sample pipeline not finished in $DIR" "error"
	echo $VERSION >> .RNASEQ_VERSION
done
case `sort -u .RNASEQ_VERSION | wc -l` in 
	1) ;;
	*) echo2 "Not all the directories were ran under the same version or condition of the single sample pipeline" "error";;
esac

########################
# running binary check #
########################
checkBin "gs"
checkBin "Rscript"
checkBin "cuffdiff"

#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
CUFFDIFF_DIR=cuffdiff_output && mkdir -p $CUFFDIFF_DIR

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
# Transcriptome GTF
TRANSCRIPTOME_GTF=$COMMON_FOLDER/${GENOME}.genes.gtf
# get depth for each of the two libraries
declare -a SAMPLE_A_NORMFACTOR SAMPLE_B_NORMFACTOR
SAMPLE_A_BAMS="" # in case enviromental variable has been set
SAMPLE_B_BAMS=""
SAMPLE_A_EXPRESS=""
SAMPLE_B_EXPRESS=""

EXPRESS_DIR_NAME=direct_transcriptome_mapping

for DIR in "${SAMPLE_A_DIR[@]}" ; do 
	echo2 $DIR
	NORMFACTOR=`cat $DIR/.*.cufflinks_depth`
	SAMPLE_A_NORMFACTOR+=( "$NORMFACTOR" )
	SAMPLE_A_BAMS=`find ${DIR}/genome_mapping/ -name "*${GENOME}.sorted.bam" `","${SAMPLE_A_BAMS}
	SAMPLE_A_EXPRESS=`find ${DIR}/$EXPRESS_DIR_NAME/ -name "*results.xprs.normalized" `" "${SAMPLE_A_EXPRESS}
done
for DIR in "${SAMPLE_B_DIR[@]}" ; do 
	NORMFACTOR=`cat $DIR/.*.cufflinks_depth`
	SAMPLE_B_NORMFACTOR+=("$NORMFACTOR")
	SAMPLE_B_BAMS=`find ${DIR}/genome_mapping/ -name "*${GENOME}.sorted.bam" `","${SAMPLE_B_BAMS}
	SAMPLE_B_EXPRESS=`find ${DIR}/$EXPRESS_DIR_NAME/ -name "*results.xprs.normalized" `" "${SAMPLE_B_EXPRESS}
done

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] RNASeq pipeline dual library mode version $RNASEQ2_VERSION"

###################################################
# using cuffdiff to perform differential analysis #
###################################################
echo2 "Running cuffdiff to call differentially expressed gene"
[ ! -f .${JOBUID}.status.${STEP}.cuffdiff ] && \
cuffdiff \
	-o $CUFFDIFF_DIR \
	-L ${SAMPLE_A_NAME},${SAMPLE_B_NAME} \
	-p $CPU \
	-u \
	--compatible-hits-norm \
	-b $GENOME_FA \
	--library-type `cat ${SAMPLE_A_DIR}/.LIBRARY_TYPE` \
	--library-norm-method geometric \
	--no-update-check \
	$TRANSCRIPTOME_GTF \
	${SAMPLE_A_BAMS} \
	${SAMPLE_B_BAMS} \
	2> $CUFFDIFF_DIR/${PREFIX}.cuffdiff.log && \
touch .${JOBUID}.status.${STEP}.cuffdiff
STEP=$((STEP+1))

####################################
# perform analysis with cummeRbund #
####################################
echo2 "Analyzing cuffdiff output with cummeRbund"
[ ! -f .${JOBUID}.status.${STEP}.cummeRbund ] && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_cummeRbund.R \
	$CUFFDIFF_DIR \
	$PDF_DIR/${PREFIX} \
	$NUM_GENE_CUMM \
	2> $PDF_DIR/${PREFIX}.cummeRbund.log && \
touch .${JOBUID}.status.${STEP}.cummeRbund
STEP=$((STEP+1))

############################################
# count transposon & piRNA cluster & genes #
############################################
echo2 "Drawing scatterplot for eXpress counting of mRNA, transposon for flies"
[ ! -f .${JOBUID}.status.${STEP}.draw_eXpress ] && \
echo -e "target_id\teff_counts" > ${SAMPLE_A_NAME}.results.xprs && \
echo -e "target_id\teff_counts" > ${SAMPLE_B_NAME}.results.xprs && \
cat $SAMPLE_A_EXPRESS | cut -f2,8 | grep -v id | sort -k1,1 | bedtools_piPipes groupby -i stdin -g 1 -c 2 -o mean >> ${SAMPLE_A_NAME}.results.xprs && \
cat $SAMPLE_B_EXPRESS | cut -f2,8 | grep -v id | sort -k1,1 | bedtools_piPipes groupby -i stdin -g 1 -c 2 -o mean >> ${SAMPLE_B_NAME}.results.xprs && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_scatter_plot_eXpress_counts.R \
	${SAMPLE_A_NAME}.results.xprs \
	${SAMPLE_B_NAME}.results.xprs \
	$SAMPLE_A_NAME \
	$SAMPLE_B_NAME \
	$PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.gene_transposon.abundance && \
	touch .${JOBUID}.status.${STEP}.draw_eXpress
STEP=$((STEP+1))
