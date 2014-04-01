
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
export SMALLRNA2_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

small RNA Seq pipeline dual library mode v$SMALLRNA2_VERSION from the $BOLD$PACKAGE_NAME$RESET
$SMALLRNA2_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

${UNDERLINE}usage${RESET}:
	piPipes small2 \	
		-a small_rna_pipeline_single_mode_output_dir1 \ 
		-b small_rna_pipeline_single_mode_output_dir2 \ 
		-g dm3 \ 
		-N siRNA \ 
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
	-g      Genome assembly name, like mm9 or dm3
		 Check "$PIPELINE_DIRECTORY/common/genome_supported.txt" for genome assemblies currently installed; 
		 Use "install" to install new genome
${OPTIONAL}[ optional ]
	-N      Normalization method, choose from " unique | uniqueXmiRNA | all | allXmiRNA | miRNA | siRNA"
		 unique:	use non-rRNA genomic unique mappers <default>.
		 siRNA:	(only for dm3) use cis-NATs and structural loci siRNA (transposon siRNAs are EXCLUDED because it might contain piRNA degradation fragments). normalized to: reads per millions of siRNA <for oxidized library in fly>.
		 uniqueXmiRNA:	use non-rRNA genomic unique mappers excluding microRNAs <for oxidized library in general>.
		 all:	use non-rRNA genomic all mappers including microRNAs.
		 allXmiRNA:	use non-rRNA genomic all mappers excluding microRNAs.
		 miRNA:	use microRNAs. normalized to: reads per millions of miRNA <for unoxidized library that can assume no change on miRNAs>.
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
while getopts "hva:b:g:c:o:A:B:N:" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		v)	echo2 "SMALLRNA2_VERSION: v$SMALLRNA2_VERSION" && exit 0 ;;
		a)	SAMPLE_A_DIR=`readlink -f $OPTARG` ;;
		b)	SAMPLE_B_DIR=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		A)  export SAMPLE_A_NAME=$OPTARG ;;
		B)  export SAMPLE_B_NAME=$OPTARG ;;
		N)	export NORMMETHOD=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
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
[[ -z $NORMMETHOD ]] && export NORMMETHOD="unique";
PREFIX=`echo -e "${SAMPLE_A_NAME}\n${SAMPLE_B_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${SAMPLE_B_NAME}. # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

#################
# Version Check #
#################
SAMPLE_A_VERSION=`ls -a $SAMPLE_A_DIR | grep SMALLRNA_VERSION`
SAMPLE_B_VERSION=`ls -a $SAMPLE_B_DIR | grep SMALLRNA_VERSION`
[ $SAMPLE_A_VERSION != $SAMPLE_B_VERSION ] && echo2 "It appears that the two runs were not done by the same assemly or same version of single library mode pipeline." "error"

#################################
# creating output files/folders #
#################################
HAIRPIN_DIR=hairpin_compare && mkdir -p $HAIRPIN_DIR
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
TRN_DIR=transposon_abundance && mkdir -p $TRN_DIR
CLUSTER_DIR=piRNA_cluster_abundance && mkdir -p $CLUSTER_DIR

########################
# running binary check #
########################
checkBin "gs"
checkBin "Rscript"

#############
# Variables #
#############
# step counter
STEP=1
# job uid
JOBUID=`echo ${PREFIX} | md5sum | cut -d" " -f1`
# directories storing the common files for this organism
export COMMON_FOLDER=$PIPELINE_DIRECTORY/common/$GENOME
# assign different values to the generalized variables (same name for different GENOMEs) according to which GENOME fed
. $COMMON_FOLDER/variables
# fasta file for the genome
export GENOME_FA=$COMMON_FOLDER/${GENOME}.fa
# chrom information of this GENOME
CHROM=$COMMON_FOLDER/${GENOME}.ChromInfo.txt
# bowtie index directory
export BOWTIE_INDEXES=$COMMON_FOLDER/BowtieIndex
# reading the information to intersect with, as well as some other annotation files
. $COMMON_FOLDER/genomic_features 
# normalization method
# unique | uniqueXmiRNA | all | allXmiRNA | miRNA | siRNA
case "$NORMMETHOD" in
unique)
	SAMPLE_A_NORMFACTOR=`head -6 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -6 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
uniquexmirna)
	SAMPLE_A_NORMFACTOR=`head -7 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -7 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
all)
	SAMPLE_A_NORMFACTOR=`head -4 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -4 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
allxmirna)
	SAMPLE_A_NORMFACTOR=`head -5 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -5 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
mirna)
	SAMPLE_A_NORMFACTOR=`head -3 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -3 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
sirna)
	case $GENOME in
	dm3)
		SAMPLE_A_NORMFACTOR=`grep -P "structural_loci|cisNATs" $SAMPLE_A_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		SAMPLE_B_NORMFACTOR=`grep -P "structural_loci|cisNATs" $SAMPLE_B_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
	;;
	*)
		echo2 "The annotation for siRNA in ${GENOME} is poor. Please choose a different normalization method. \nIf unox library, choose \"miRNA\". If ox, choose \"uniquexmirna\"" "error"
		# SAMPLE_A_NORMFACTOR=`grep -P "refSeq_GENE" $SAMPLE_A_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/# a}'`
		# SAMPLE_B_NORMFACTOR=`grep -P "refSeq_GENE" $SAMPLE_B_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/a}'`	
	;;
	esac
	
;;
*)
	echo2 "unrecognized normalization option: $NORMMETHOD; using the default method" "warning"
	SAMPLE_A_NORMFACTOR=`head -6 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	SAMPLE_B_NORMFACTOR=`head -6 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
esac

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] small RNA pipeline dual library mode version $SMALLRNA2_VERSION"

#############################
# draw microRNA balloonplot #
#############################
echo2 "Drawing microRNA balloonplot"
SAMPLE_A_MIR_BED2=$SAMPLE_A_DIR/hairpins_mapping/*bed2.relative
SAMPLE_B_MIR_BED2=$SAMPLE_B_DIR/hairpins_mapping/*bed2.relative
[ ! -f .${JOBUID}.status.${STEP}.miRNA_balloonplot.normalized_by_$NORMMETHOD ] && \
	mergeMirBed2 $SAMPLE_A_MIR_BED2 $SAMPLE_B_MIR_BED2 $SAMPLE_A_NORMFACTOR $SAMPLE_B_NORMFACTOR > $HAIRPIN_DIR/${PREFIX}.miRNA.relative.abundance.normalized_by_$NORMMETHOD && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_balloonPlot.R $HAIRPIN_DIR/${PREFIX}.miRNA.relative.abundance.normalized_by_$NORMMETHOD ${SAMPLE_A_NAME} ${SAMPLE_B_NAME} $CPU $HAIRPIN_DIR 1>&2 2>$HAIRPIN_DIR/${PREFIX}.miRNA.balloonplot.normalized_by_$NORMMETHOD.log && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.miRNAballoon.normalized_by_$NORMMETHOD.pdf $HAIRPIN_DIR/*miRNAballoonPlot.pdf && \
	rm -rf $HAIRPIN_DIR/*miRNAballoonPlot.pdf && \
	touch .${JOBUID}.status.${STEP}.miRNA_balloonplot.normalized_by_$NORMMETHOD 
STEP=$((STEP+1))

###################################################
# draw piRNA abundance for each transposon family #
###################################################
ALL_BED2_A=`ls $SAMPLE_A_DIR/genome_mapping/*${GENOME}*all.piRNA.bed2`
ALL_BED2_B=`ls $SAMPLE_B_DIR/genome_mapping/*${GENOME}*all.piRNA.bed2`
echo2 "Calculating piRNA abundances for each transposon family"
case $GENOME in
dm3)
	[ ! -f .${JOBUID}.status.${STEP}.transposon_abundance ] && \
	para_file=$TRN_DIR/${RANDOM}${RANDOM}.para && \
	echo "bedtools_piPipes intersect -split -wo -f 0.99 -a $ALL_BED2_A -b $Trn | awk -v depth=$SAMPLE_A_NORMFACTOR '{split(\$11,name,\".\"); if (\$6==\$13) a[name[1]]+=\$4/\$NF/\$5; else b[name[1]]+=\$4/\$NF/\$5; total[name[1]]=1; class[name[1]]=\$12;}END{for (c in total) { printf \"%s\t%s\t%.2f\t%.2f\n\", c, class[c], (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $TRN_DIR/${SAMPLE_A_NAME}.transposon.abundance " >> $para_file  && \
	echo "bedtools_piPipes intersect -split -wo -f 0.99 -a $ALL_BED2_B -b $Trn | awk -v depth=$SAMPLE_B_NORMFACTOR '{split(\$11,name,\".\"); if (\$6==\$13) a[name[1]]+=\$4/\$NF/\$5; else b[name[1]]+=\$4/\$NF/\$5; total[name[1]]=1; class[name[1]]=\$12;}END{for (c in total) { printf \"%s\t%s\t%.2f\t%.2f\n\", c, class[c], (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $TRN_DIR/${SAMPLE_B_NAME}.transposon.abundance " >> $para_file  && \
	ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	Rscript --slave $PIPELINE_DIRECTORY/bin/piPipes_draw_scatter_plot.R  $TRN_DIR/${SAMPLE_A_NAME}.transposon.abundance $TRN_DIR/${SAMPLE_B_NAME}.transposon.abundance  $SAMPLE_A_NAME $SAMPLE_B_NAME $PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.transposon.abundance && \
	touch .${JOBUID}.status.${STEP}.transposon_abundance
;;
*)	
	[ ! -f .${JOBUID}.status.${STEP}.transposon_abundance ] && \
	para_file=$TRN_DIR/${RANDOM}${RANDOM}.para && \
	echo "bedtools_piPipes intersect -wo -f 0.99 -a $ALL_BED2_A -b $repeatMasker  | awk -v depth=$SAMPLE_A_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $TRN_DIR/${SAMPLE_A_NAME}.transposon.abundance " >> $para_file  && \
	echo  "bedtools_piPipes intersect -wo -f 0.99 -a $ALL_BED2_B -b $repeatMasker | awk -v depth=$SAMPLE_B_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $TRN_DIR/${SAMPLE_B_NAME}.transposon.abundance " >> $para_file  && \
	ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	Rscript --slave $PIPELINE_DIRECTORY/bin/piPipes_draw_scatter_plot.R  $TRN_DIR/${SAMPLE_A_NAME}.transposon.abundance $TRN_DIR/${SAMPLE_B_NAME}.transposon.abundance  $SAMPLE_A_NAME $SAMPLE_B_NAME $PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.transposon.abundance && \
	touch .${JOBUID}.status.${STEP}.transposon_abundance
;;
esac
STEP=$((STEP+1))

###############################################
# draw piRNA abundance for each piRNA cluster #
###############################################
echo2 "Calculating piRNA abundances for each piRNA cluster family"
[ ! -f .${JOBUID}.status.${STEP}.piRNA_cluster_abundance.normalized_by_$NORMMETHOD ] && \
para_file=$CLUSTER_DIR/${RANDOM}${RANDOM}.para && \
	echo "bedtools_piPipes intersect -split -wo -f 0.99 -a $ALL_BED2_A -b $piRNA_Cluster | awk -v depth=$SAMPLE_A_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $CLUSTER_DIR/${SAMPLE_A_NAME}.cluster.abundance.normalized_by_$NORMMETHOD " >> $para_file  && \
	echo "bedtools_piPipes intersect -split -wo -f 0.99 -a $ALL_BED2_B -b $piRNA_Cluster | awk -v depth=$SAMPLE_B_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > $CLUSTER_DIR/${SAMPLE_B_NAME}.cluster.abundance.normalized_by_$NORMMETHOD " >> $para_file  && \
	ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	Rscript --slave $PIPELINE_DIRECTORY/bin/piPipes_draw_scatter_plot.R  $CLUSTER_DIR/${SAMPLE_A_NAME}.cluster.abundance.normalized_by_$NORMMETHOD $CLUSTER_DIR/${SAMPLE_B_NAME}.cluster.abundance.normalized_by_$NORMMETHOD  $SAMPLE_A_NAME $SAMPLE_B_NAME $PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.piRNAcluster.abundance.normalized_by_$NORMMETHOD && \
	touch .${JOBUID}.status.${STEP}.piRNA_cluster_abundance.normalized_by_$NORMMETHOD
STEP=$((STEP+1))

############################################
# count transposon & piRNA cluster & genes #
############################################
echo2 "Drawing scatterplot for eXpress counting of mRNA, transposon and cluster"
[ ! -f .${JOBUID}.status.${STEP}.draw_eXpress.normalized_by_$NORMMETHOD ] && \
awk -v depth=$SAMPLE_A_NORMFACTOR 'BEGIN{FS=OFS="\t"; getline; print }{$8*=depth; print }' $SAMPLE_A_DIR/eXpress_quantification_no_normalization/results.xprs > $SAMPLE_A_DIR/eXpress_quantification_no_normalization/results.xprs.normalized_by_${NORMMETHOD} && \
awk -v depth=$SAMPLE_B_NORMFACTOR 'BEGIN{FS=OFS="\t"; getline; print }{$8*=depth; print }' $SAMPLE_B_DIR/eXpress_quantification_no_normalization/results.xprs > $SAMPLE_B_DIR/eXpress_quantification_no_normalization/results.xprs.normalized_by_${NORMMETHOD} && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_scatter_plot_eXpress_counts.R \
	$SAMPLE_A_DIR/eXpress_quantification_no_normalization/results.xprs.normalized_by_${NORMMETHOD} \
	$SAMPLE_B_DIR/eXpress_quantification_no_normalization/results.xprs.normalized_by_${NORMMETHOD} \
	$SAMPLE_A_NAME \
	$SAMPLE_B_NAME \
	$PDF_DIR/${SAMPLE_A_NAME}_vs_${SAMPLE_B_NAME}.gene_transposon_cluster.normalized_by_${NORMMETHOD}.abundance && \
	touch .${JOBUID}.status.${STEP}.draw_eXpress.normalized_by_$NORMMETHOD
STEP=$((STEP+1))






