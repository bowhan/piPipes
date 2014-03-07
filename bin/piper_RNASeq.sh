
# RNASeq pipeline single library mode
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
export RNASEQ_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

RNASeq pipeline single library mode v$RNASEQ_VERSION from the $BOLD$PACKAGE_NAME$RESET
$RNASEQ_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

${UNDERLINE}usage${RESET}:
	piper rna \ 
		-l left.fq \ 
		-r right.fq \ 
		-g dm3 \ 
		-o output_directory [current working directory] \ 
		-c cpu [8] 
	
Currently, only Paired-End input is accepted. By default, the pipeline assumes dUTR based method 
and \2 is in the same direction as the transcript, opposite to ligation-based degradome/cage-Seq. 
Use option -L to change the behavior.

OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
	-l      Left reads from Paired-End sequencing
	-r      Right reads from Paired-End sequencing
	-g      Genome assembly name, like mm9 or dm3. required
		 Check $PIPELINE_DIRECTORY/common/genome_supported.txt for genome assemblies currently installed; 
		 Use "install" to install new genome
${OPTIONAL}[ optional ]
	-L      Ligation based library preperation method; 
		 Left reads (\1) being in the same direction as the transcripts; 
		 default: off (dUTR based, \2 reads being in the same direction)
	-o      Output directory, default: current directory $PWD
	-c      Number of CPUs to use, default: 8
	
EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hl:r:c:o:g:vL" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		l)	LEFT_FASTQ=`readlink -f $OPTARG` ;;
		r)	RIGHT_FASTQ=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		v)	echo2 "RNASEQ_VERSION: v$RNASEQ_VERSION" && exit 0 ;;
		L)	LIGATIONLIB=1 ;; # ligation based
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
[ ! -z $OUTDIR ] || OUTDIR=$PWD # if -o is not specified, use current directory
[ "$OUTDIR" != `readlink -f $PWD` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "warning")
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
touch .writting_permission && rm -rf .writting_permission || (echo2 "Cannot write in directory ${OUTDIR}... Exiting..." "error")

if [ "${LIGATIONLIB}" == 1 ]; 
then 
	LIBRARY_TYPE="fr-secondstrand";
	END_TO_REVERSE_STRAND=2
	SENSE_HTSEQ_OPT="yes"; 
	ANTISENSE_HTSEQ_OPT="reverse"; 
	EXPRESS_OPTION="--fr-stranded"
else 
	LIBRARY_TYPE="fr-firststrand"
	END_TO_REVERSE_STRAND=1
	SENSE_HTSEQ_OPT="reverse"; 
	ANTISENSE_HTSEQ_OPT="yes"; 
	EXPRESS_OPTION="--rf-stranded"
fi
echo $LIBRARY_TYPE > .LIBRARY_TYPE

#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
READS_DIR=input_read_files && mkdir -p $READS_DIR 
rRNA_DIR=rRNA_mapping && mkdir -p $rRNA_DIR
GENOMIC_MAPPING_DIR=genome_mapping && mkdir -p $GENOMIC_MAPPING_DIR
CUFFLINKS_DIR=cufflinks_output && mkdir -p $CUFFLINKS_DIR
HTSEQ_DIR=htseq_count && mkdir -p $HTSEQ_DIR
DIRECTMAPPING_DIR=gene_transposon_cluster_direct_mapping && mkdir -p $DIRECTMAPPING_DIR
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
checkBin "STAR"
checkBin "ParaFly"
checkBin "bedtools_piper"
checkBin "bedGraphToBigWig"
checkBin "express"
checkBin "cufflinks"
checkBin "htseq-count" # the user need to install this separately

#############
# Variables #
#############
# step counter
STEP=1
# job uid
JOBUID=`echo ${INPUT_FASTQ} | md5sum | cut -d" " -f1`
LEFT_FASTQ_NAME=`basename $LEFT_FASTQ`
RIGHT_FASTQ_NAME=`basename $RIGHT_FASTQ`
PREFIX=`echo -e "${LEFT_FASTQ_NAME}\n${RIGHT_FASTQ_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && export PREFIX=${LEFT_FASTQ_NAME%.f[aq]} # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
# table to store the basic statistics of the library (genomic mappability)
TABLE=${PREFIX}.basic_stats
# directories storing the common files for this organism
export COMMON_FOLDER=$PIPELINE_DIRECTORY/common/$GENOME
# assign different values to the generalized variables (same name for different GENOMEs) according to which GENOME fed
. $COMMON_FOLDER/variables
# fasta file for the genome
export GENOME_FA=$COMMON_FOLDER/${GENOME}.fa
# chrom information of this GENOME
CHROM=$COMMON_FOLDER/${GENOME}.ChromInfo.txt
# Transcriptome GTF
TRANSCRIPTOME_GTF=$COMMON_FOLDER/${GENOME}.genes.gtf
# exporting BOWTIE2_INDEXES
export BOWTIE2_INDEXES=$COMMON_FOLDER/Bowtie2Index
# STAR index for the genome
STARINDEX=$COMMON_FOLDER/STARIndex

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] RNA-Seq pipeline version $RNASEQ_VERSION" 

###########################
# determine fastQ version #
###########################
echo2 "Determining the version of fastQ using SolexaQA"
# determine version of fastq used, using a modified SolexaQA.pl
PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA.pl ${LEFT_FASTQ}`
case ${PHRED_SCORE} in
solexa)		bowtie2PhredOption="--solexa-quals" && STARoutQSconversionAdd="-31" ;; # Solexa+64, raw reads typically (-5, 40)
illumina)	bowtie2PhredOption="--phred64" && STARoutQSconversionAdd="-31" ;; # Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
sanger)		bowtie2PhredOption="--phred33" && STARoutQSconversionAdd="0" ;; # Phred+33,  raw reads typically (0, 40) (http://en.wikipedia.org/wiki/FASTQ_format)
*)			echo2 "unable to determine the fastq version. Using sanger..." "warning";;
esac

############################
# map to rRNA with bowtie2 #
############################
echo2 "Mapping input reads to rRNA with Bowtie2"
[ ! -f .${JOBUID}.status.${STEP}.rRNA_mapping ] && \
bowtie2 \
	-x rRNA \
	-1 $LEFT_FASTQ \
	-2 $RIGHT_FASTQ \
	-q \
	$bowtie2PhredOption \
	--very-fast \
	-k 1 \
	--no-mixed \
	--no-discordant \
	--un-conc ${READS_DIR}/${PREFIX}.x_rRNA.fq \
	-p $CPU \
	-S /dev/null \
	2> ${rRNA_DIR}/${PREFIX}.rRNA.log && \
touch .${JOBUID}.status.${STEP}.rRNA_mapping
[ ! -f .${JOBUID}.status.${STEP}.rRNA_mapping ] && echo2 "Failed to map to rRNA" "error"
STEP=$((STEP+1))
InputReads=`head -2 ${rRNA_DIR}/${PREFIX}.rRNA.log | tail -1 | awk '{print $1}'`
rRNAReads=`head -4 ${rRNA_DIR}/${PREFIX}.rRNA.log | tail -1 | awk '{print $1}'`
echo -e "total_input_reads:\t${InputReads}" > $TABLE
echo -e "rRNA_reads:\t${rRNAReads}" >> $TABLE

###########################
# map to genome with STAR #
###########################
echo2 "Mapping non-rRNA reads to genome $GENOME with STAR"
xrRNA_LEFT_FQ=${READS_DIR}/${PREFIX}.x_rRNA.1.fq && \
xrRNA_RIGHT_FQ=${READS_DIR}/${PREFIX}.x_rRNA.2.fq && \
[ ! -f .${JOBUID}.status.${STEP}.genome_mapping ] && \
STAR \
	--runMode alignReads \
	--genomeDir $STARINDEX \
	--readFilesIn ${xrRNA_LEFT_FQ} ${xrRNA_RIGHT_FQ} \
	--runThreadN $CPU \
	--outFilterScoreMin 0 \
	--outFilterScoreMinOverLread 0.72 \
	--outFilterMatchNmin 0 \
	--outFilterMatchNminOverLread 0.72 \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax -1 \
	--outFilterMismatchNmax 10 \
	--outFilterMismatchNoverLmax 0.05 \
	--alignIntronMax 0 \
	--alignIntronMin 21 \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--genomeLoad NoSharedMemory \
	--outFileNamePrefix $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.${GENOME}. \
	--outSAMunmapped None \
	--outReadsUnmapped Fastx \
	--outSJfilterReads Unique \
	--seedSearchStartLmax 20 \
	--seedSearchStartLmaxOverLread 1.0 \
	--chimSegmentMin 0 2>&1 1> $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.${GENOME}.STAR.log && \
touch .${JOBUID}.status.${STEP}.genome_mapping
[ ! -f .${JOBUID}.status.${STEP}.genome_mapping ] && echo2 "Failed to map to genome" "error"
STEP=$((STEP+1))

# getting statistics
InputReads=`grep 'Number of input reads' $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.${GENOME}.Log.final.out | awk '{print $NF}'`
UniquReads=`grep 'Uniquely mapped reads number' $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.${GENOME}.Log.final.out | awk '{print $NF}'`
MultiReads=`grep 'Number of reads mapped to multiple loci' $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.${GENOME}.Log.final.out | awk '{print $NF}'`
AllMapReads=$((UniquReads+MultiReads))
UnMapReads=$((InputReads-UniquReads-MultiReads))
echo -e "genomie_mapper_reads:\t${AllMapReads}" >> $TABLE
echo -e "genomie_unique_mapper_reads:\t${UniquReads}" >> $TABLE
echo -e "genomie_multiple_mapper_reads:\t${MultiReads}" >> $TABLE
echo -e "genomie_unmappable_reads:\t${UnMapReads}" >> $TABLE

#######################
# Processing sam file #
#######################
echo2 "Processing mapping results"
[ ! -f .${JOBUID}.status.${STEP}.genome_bam_processing ] && \
	samtools view -bS ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam 2>/dev/null | tee ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam | bedtools_piper bamtobed -bedpe -mate1 -i - > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.bedpe && \
	samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted && \
	samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam && \
	awk 'BEGIN{FS=OFS="\t"}{if (ARGIND==1) {++ct[$7]} else {$8=1.0/ct[$7]; print $0 > "/dev/stdout"; if (ct[$7]==1) print $0 > "/dev/stderr"}}' ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.bedpe ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.bedpe \
		1> ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.all.normalized.bedpe \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.unique.normalized.bedpe && \
	rm -rf  ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam  ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.bedpe && \
	touch .${JOBUID}.status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

###################################################
# Transcripts abundance estimation with Cufflinks #
###################################################
echo2 "Estimating transripts abundance with cufflinks --compatible-hits-norm"
# the reason we use "--compatible-hits-norm" is because that transposon reads could sometimes change the sequencing depth
# we assume that the mRNA transcritome is less influenced by the piRNA pathway
# to mRNA comptatible reads better reflect the depth of the library and the change of transposon RNA
[ ! -f .${JOBUID}.status.${STEP}.quantification_by_cuff_compatible_hits_norm ] && \
cufflinks \
	-o $CUFFLINKS_DIR \
	-p $CPU \
	-G ${TRANSCRIPTOME_GTF} \
	-b $GENOME_FA \
	-u \
	--library-type $LIBRARY_TYPE \
	--compatible-hits-norm \
	--no-update-check \
	${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam \
	2> $CUFFLINKS_DIR/${PREFIX}.cufflinks.log && \
	touch .${JOBUID}.status.${STEP}.quantification_by_cuff_compatible_hits_norm
STEP=$((STEP+1))
MapMass=`grep 'Normalized Map Mass' $CUFFLINKS_DIR/${PREFIX}.cufflinks.log | cut -d' ' -f4`
export NormScale=`echo $MapMass | awk '{printf "%f",1000000.0/$1}'`
echo -ne "$NormScale" > .${JOBUID}.cufflinks_depth

####################################
# Making bigWig for Genome Browser #
####################################
# in order to make bigWig for unique mappers, we need to reverse the strand of one of the end
echo2 "Making bigWig from sorted bam"
[ ! -f .${JOBUID}.status.${STEP}.make_bigWig ] && \
	bedtools_piper bamtobed -bed12 -tag NH -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
	awk -v strand=$END_TO_REVERSE_STRAND 'BEGIN{FS=OFS="\t"}{e=substr($4,length($4)); if (e==strand) $6=($6=="+"?"-":"+"); if ($5==1) print $0; }' > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 && \
	para_file=${BW_OUTDIR}/${RANDOM}${RANDOM}.makeBigWigPE.para && \
	echo "bedtools_piper genomecov -scale $NormScale -split -bg -strand + -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bigWig"  >  $para_file && \
	echo "bedtools_piper genomecov -scale $NormScale -split -bg -strand - -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bigWig" >> $para_file && \
	ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	rm -rf ${para_file} ${para_file}.completed && \
	touch .${JOBUID}.status.${STEP}.make_bigWig
STEP=$((STEP+1))

#############################################
# genomic feature counting with htSeq-count #
#############################################
echo2 "Quantifying genomic features from genomic mapping using HTSeq-count"
[ ! -f .${JOBUID}.status.${STEP}.htseq_count ] && \
	. $COMMON_FOLDER/genomic_features && \
	[ ! -z $HTSEQ_TARGETS ] && \
	para_file=$HTSEQ_DIR/${RANDOM}${RANDOM}.para && \
	for t in "${HTSEQ_TARGETS[@]}"; do \
		echo "htseq-count -m intersection-strict -s $SENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.strict.S.out" >> $para_file
		echo "htseq-count -m intersection-strict -s $ANTISENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.strict.AS.out" >> $para_file
		echo "htseq-count -m union -s $SENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.union.S.out" >> $para_file
		echo "htseq-count -m union -s $ANTISENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.union.AS.out" >> $para_file
done && \
ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam && \
touch .${JOBUID}.status.${STEP}.htseq_count
STEP=$((STEP+1))

##################################################
# direct mapping and quantification with eXpress #
##################################################
echo2 "Mapping to genes, transposon and piRNA cluster directly with Bowtie2"
. $COMMON_FOLDER/genomic_features
[ ! -f .${JOBUID}.status.${STEP}.direct_mapping ] && \
bowtie2 -x gene+cluster+repBase \
	-1 ${xrRNA_LEFT_FQ} \
	-2 ${xrRNA_RIGHT_FQ} \
	-q \
	$bowtie2PhredOption \
	-a \
	-X 800 \
	--no-mixed \
	--quiet \
	-p $CPU \
	2> ${DIRECTMAPPING_DIR}/${PREFIX}.gene+cluster+repBase.log | \
	express -o $DIRECTMAPPING_DIR --no-update-check --library-size ${MapMass%.*} $COMMON_FOLDER/${GENOME}.gene+cluster+repBase.fa 1>&2 2> $DIRECTMAPPING_DIR/${PREFIX}.gene+cluster+repBase.eXpress.log && \
	touch .${JOBUID}.status.${STEP}.direct_mapping
STEP=$((STEP+1))

#############
# finishing #
#############
echo2 "Finished running ${PACKAGE_NAME} RNA-Seq pipeline version $RNASEQ_VERSION"
echo2 "---------------------------------------------------------------------------------"
touch .${GENOME}.${LIGATIONLIB}.RNASEQ_VERSION.${RNASEQ_VERSION}








