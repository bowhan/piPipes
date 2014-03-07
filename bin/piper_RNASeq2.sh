
# RNA-Seq pipeline dual mode
# pipipe 
# https://github.com/bowhan/pipipe.git
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
export RNASEQ2_VERSION=1.0.0

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
	pipipe rna2	\ 
		-a rnaseq_pipiline_output_dir1 \ 
		-b rnaseq_pipiline_output_dir2 \ 
		-g dm3 \ 
		-c 24 [8] \ 
		-o output_dir [cwd] \ 
		-A piwi_heterozygous [basename of -a] \ 
		-B piwi_mutant [basename of -b]
		
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

The pipeline will automatically detect the version and options of the single library run for the two samples and ensure the consistence. 

EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hva:b:g:c:o:A:B:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		v)	echo2 "SMALLRNA2_VERSION: v$SMALLRNA2_VERSION" && exit 0 ;;
		a)	SAMPLE_A_DIR=`readlink -f $OPTARG` ;;
		b)	SAMPLE_B_DIR=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		A)  SAMPLE_A_NAME=$OPTARG ;;
		B)  SAMPLE_B_NAME=$OPTARG ;;
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
[[ -z $SAMPLE_A_NAME ]] && SAMPLE_A_NAME=`basename $SAMPLE_A_DIR`
[[ -z $SAMPLE_B_NAME ]] && SAMPLE_B_NAME=`basename $SAMPLE_B_DIR`
PREFIX=`echo -e "${SAMPLE_A_NAME}\n${SAMPLE_B_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${SAMPLE_B_NAME} # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################
# Version Check #
#################
SAMPLE_A_VERSION=`ls -a $SAMPLE_A_DIR | grep RNASEQ_VERSION`
SAMPLE_B_VERSION=`ls -a $SAMPLE_B_DIR | grep RNASEQ_VERSION`
[ "$SAMPLE_A_VERSION" != "$SAMPLE_B_VERSION" ] && echo2 "It appears that the two runs were not done by the same assemly or same version of single library mode pipeline." "error"

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
SAMPLE_A_NORMFACTOR=`cat $SAMPLE_A_DIR/.*.cufflinks_depth`
SAMPLE_B_NORMFACTOR=`cat $SAMPLE_B_DIR/.*.cufflinks_depth`

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
	${SAMPLE_A_DIR}/genome_mapping/*${GENOME}.sorted.bam \
	${SAMPLE_B_DIR}/genome_mapping/*${GENOME}.sorted.bam \
	2> $CUFFDIFF_DIR/${PREFIX}.cuffdiff.log && \
touch .${JOBUID}.status.${STEP}.cuffdiff
STEP=$((STEP+1))

####################################
# perform analysis with cummeRbund #
####################################
echo2 "Analyzing cuffdiff output with cummeRbund"
[ ! -f .${JOBUID}.status.${STEP}.cummeRbund ] && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/pipipe_cummeRbund.R \
	$CUFFDIFF_DIR \
	$PDF_DIR/${PREFIX} \
	2> $PDF_DIR/${PREFIX}.cummeRbund.log && \
touch .${JOBUID}.status.${STEP}.cummeRbund
STEP=$((STEP+1))

############################################
# count transposon & piRNA cluster & genes #
############################################
echo2 "Drawing scatterplot for eXpress counting of mRNA, transposon and cluster"
[ ! -f .${JOBUID}.status.${STEP}.draw_eXpress ] && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/pipipe_draw_scatter_plot_eXpress_counts.R \
	$SAMPLE_A_DIR/gene_transposon_cluster_direct_mapping/results.xprs \
	$SAMPLE_B_DIR/gene_transposon_cluster_direct_mapping/results.xprs \
	$SAMPLE_A_NAME \
	$SAMPLE_B_NAME \
	$PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.gene_transposon_cluster.abundance && \
	touch .${JOBUID}.status.${STEP}.draw_eXpress
STEP=$((STEP+1))









