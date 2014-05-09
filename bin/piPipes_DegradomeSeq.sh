
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
export DEG_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

Degradome Seq pipeline single library mode v$DEG_VERSION from the $BOLD$PACKAGE_NAME$RESET
$CAGE_DEG_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

${UNDERLINE}usage${RESET}:
	piPipes deg \ 
		-l left.fq \ 
		-r right.fq \ 
		-g dm3 \ 
		-s small_RNA_pipeline_output \ 
		-o output_directory [current directory] \ 
		-c cpu [8] 
	
Currently, only Paired-End input is accepted. And \1 should be in the same direction as the transcripts, 
opposite to dUTR based RNASeq. 

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
	-s      small RNA pipeline output; if this option is provided, the pipeline check the Ping-Pong signiture between the two libraries.
	-o      Output directory, default: current directory $PWD
	-c      Number of CPUs to use, default: 8
	

EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hl:r:c:o:g:s:v" OPTION; do
	case $OPTION in
		h)	usage && exit 1 ;;
		l)	LEFT_FASTQ=`readlink -f $OPTARG` ;;
		r)	RIGHT_FASTQ=`readlink -f $OPTARG` ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
		v)	echo2 "DEG_VERSION: v$DEG_VERSION" && exit 0 ;;
		s)	export SRA_LIB_DIR=`readlink -f $OPTARG` ;;
		*)	usage && exit 1 ;;
	esac
done
# if INPUT_FASTQ or GENOME is undefined, print out usage and exit
[[ -z $LEFT_FASTQ ]] && usage && echo2 "Missing option -l for input fastq of left file, or file does not exist " "error" 
[[ -z $RIGHT_FASTQ ]] && usage && echo2 "Missing option -r for input fastq of right file, or file does not exist " "error" 
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

if [ ! -z "$SRA_LIB_DIR" ]; then
	[ ! -d "$SRA_LIB_DIR" ] && echo2 "directory $SRA_LIB_DIR not exist" "error"
	ls -a $SRA_LIB_DIR | grep SMALLRNA_VERSION 1>/dev/null 2>/dev/null || echo2 "seems that the small RNA pipeline was not finished normally. Either re-run it or run this degradome pipeline without -s option" "error"
	SRA_UNIQ_BED2=`find $SRA_LIB_DIR/genome_mapping -name "*unique.bed2"`
	[ -z "$SRA_UNIQ_BED2" ] && echo2 "failed to locate the \"unique.bed2\" file in the $SRA_LIB_DIR/genome_mapping directory." "error"
	[[ $SRA_UNIQ_BED2 = *$'\n'* ]] && echo2 "more than one files been found to match *unique.bed2. Seems that you have modified the folder $SRA_LIB_DIR/genome_mapping. Please remove additinal file and only keep one *unique.bed2" "error"
	export SRA_UNIQ_BED2
fi

# degradome options
SENSE_HTSEQ_OPT="yes"; 
ANTISENSE_HTSEQ_OPT="reverse"; 
	
#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
READS_DIR=input_read_files && mkdir -p $READS_DIR 
rRNA_DIR=rRNA_mapping && mkdir -p $rRNA_DIR
GENOMIC_MAPPING_DIR=genome_mapping && mkdir -p $GENOMIC_MAPPING_DIR
CUFFLINKS_DIR=cufflinks_output && mkdir -p $CUFFLINKS_DIR
HTSEQ_DIR=htseq_count && mkdir -p $HTSEQ_DIR
BEDTOOLS_DIR=bedtools_count && mkdir -p $BEDTOOLS_DIR
DIRECTMAPPING_DIR=gene_transposon_cluster_direct_mapping && mkdir -p $DIRECTMAPPING_DIR
SUMMARY_DIR=summaries && mkdir -p $SUMMARY_DIR
BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR
INDEX_OUTDIR=bowtie_index && mkdir -p $INDEX_OUTDIR

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
checkBin "bedtools_piPipes"
checkBin "bedGraphToBigWig"
checkBin "express"
checkBin "cufflinks"
checkBin "htseq-count" # the user need to install this separately
checkBin "piPipes_filter_CIGAR"

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
READLEN=`head -2 $LEFT_FASTQ | awk '{getline; print length($0)}'`
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
echo2 "Beginning running [${PACKAGE_NAME}] Degradome-Seq pipeline version $DEG_VERSION" 

###########################
# determine fastQ version #
###########################
echo2 "Determining the version of fastQ using SolexaQA"
# determine version of fastq used, using a modified SolexaQA.pl
PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piPipes.pl ${LEFT_FASTQ}`
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
NormScale=`echo ${UniquReads} | awk '{printf "%f",1000000.0/$1}'`

#######################
# Processing sam file #
#######################
echo2 "Processing mapping results"
# we dump the ones with softclipping on the 5' end
[ ! -f .${JOBUID}.status.${STEP}.genome_bam_processing ] && \
	piPipes_filter_CIGAR -5 -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam | \
	samtools view -uS -f0x2 - 2>/dev/null | \
	samtools sort -o -@ $CPU - foo > \
	${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam && \
	samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam && \
	touch .${JOBUID}.status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

####################################
# Making bigWig for Genome Browser #
####################################
# in order to make bigWig for unique mappers, we need to reverse the strand of one of the end
# we only takes the \1 from the bam 
echo2 "Making bigWig from sorted bam \1 reads without 5' soft-clipping"
[ ! -f .${JOBUID}.status.${STEP}.make_bigWig ] && \
	samtools view -huf0x40 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
	bedtools bamtobed -bed12 -tag NH -i - | \
	awk 'BEGIN{FS=OFS="\t"}{ if ($5==1) print $0 > "/dev/stderr"; $4=1; $5=1.0/$5; print $0 > "/dev/stdout"}' \
		2> ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.bed12 | \
	sort -k1,3 -k6,12 --parallel=$CPU | \
	bedtools_piPipes groupby -i - -g 1,2,3,6,7,8,9,10,11,12 -c 4,5 -o sum,sum | \
	awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$11,$12,$4,$5,$6,NR,$8,$9,$10}' > \
		${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.all.bed12 && \
	paraFile=${BW_OUTDIR}/${RANDOM}${RANDOM}.makeBigWig.para && \
	echo "bedtools_piPipes genomecov -5 -scale $NormScale -bg -strand + -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.bed12 -g $CHROM > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Watson.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Watson.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Watson.bigWig"  >  $paraFile && \
	echo "bedtools_piPipes genomecov -5 -scale $NormScale -bg -strand - -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.bed12 -g $CHROM | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Crick.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Crick.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.unique.Crick.bigWig" >> $paraFile && \
	ParaFly -c $paraFile -CPU $CPU -failed_cmds ${paraFile}.failedCommands 1>&2 && \
	rm -rf ${paraFile} ${paraFile}.completed && \
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

##########################################
# genomic feature counting with bedtools #
##########################################
echo2 "Quantifying genomic features from genomic mapping using BEDTools"
[ ! -f .${JOBUID}.status.${STEP}.bedtools ] && \
	bash $DEBUG piPipes_intersect_degradome_with_genomic_features.sh \
		${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.f0x40.noS.all.bed12 \
		$SUMMARY_DIR/${PREFIX}.summary \
		$CPU \
		$BEDTOOLS_DIR && \
	touch .${JOBUID}.status.${STEP}.bedtools
STEP=$((STEP+1))

#################################################
#direct mapping and quantification with eXpress #
#################################################
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
samtools view -bS - \
	> ${DIRECTMAPPING_DIR}/${PREFIX}.gene+cluster+repBase.bam && \
touch .${JOBUID}.status.${STEP}.direct_mapping

[ ! -f .${JOBUID}.status.${STEP}.eXpress_quantification ] && \
express \
	-B 21 \
	-o $DIRECTMAPPING_DIR \
	--no-update-check \
	--library-size ${AllMapReads} \
	$COMMON_FOLDER/${GENOME}.gene+cluster+repBase.fa \
	${DIRECTMAPPING_DIR}/${PREFIX}.gene+cluster+repBase.bam \
	1>&2 2> $DIRECTMAPPING_DIR/${PREFIX}.gene+cluster+repBase.eXpress.log && \
touch .${JOBUID}.status.${STEP}.eXpress_quantification
STEP=$((STEP+1))

##########################################################
# building bowtie index from left read for piRNA mapping #
##########################################################
[ ! -f .${JOBUID}.status.${STEP}.generate_bowtie_index ] && \
	echo ">${PREFIX}" > $INDEX_OUTDIR/${PREFIX}.r1.unique.fa && \
	samtools view -hb -f0x40 -F0x100 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
	samtools bam2fq - | \
	awk '{getline; if (!printed[$0]) { print $0; printed[$0]=1;} getline; getline }' >> $INDEX_OUTDIR/${PREFIX}.r1.unique.fa && \
	bowtie-build $INDEX_OUTDIR/${PREFIX}.r1.unique.fa $INDEX_OUTDIR/${PREFIX} 1>&2 && \
	touch .${JOBUID}.status.${STEP}.generate_bowtie_index
STEP=$((STEP+1))

if [ -n $SRA_UNIQ_BED2 ]; then 
	[ ! -f .${JOBUID}.status.${STEP}.map_smRNA_to_deg_index ] && \
		awk '{if (!printed[$7]) {print ">"$7"_"$4"\n"$7; printed[$7]=1}}' $SRA_UNIQ_BED2 | \
		bowtie -f -a --best --strata -v 0 -p $CPU -S $INDEX_OUTDIR/${PREFIX} - 2> $INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.log | \
		samtools view -bS - | \
		bedtools_piPipes bamtobed -i > \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed1 && \
		awk 'BEGIN{FS=OFS="\t"}{if (ARGIND==1) ++ct[$4]; else {split ($4,ar,"_");  print $1,$2,$3,ar[2],ct[$4],$6,ar[1]}}' \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed1 $INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed1 > \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed2 && \
		rm -rf $INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed1 && \
		awk -v read_len=$READLEN '{if ($6=="+") {ct[$2%read_len]+=1.0/$5}}   END{for (a=0;a<=(read_len-1);++a){printf "%d\t%.2f\n",a,ct[a]}}' \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed2 > \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.species.5end && \
		awk -v read_len=$READLEN '{if ($6=="+") {ct[$2%read_len]+=1.0*$4/$5}}END{for (a=0;a<=(read_len-1);++a){printf "%d\t%.2f\n",a,ct[a]}}' \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.bed2 > \
			$INDEX_OUTDIR/`basename $SRA_UNIQ_BED2`.map_to.${PREFIX}.reads.5end && \
		touch .${JOBUID}.status.${STEP}.map_smRNA_to_deg_index
	STEP=$((STEP+1))
fi

#############
# finishing #
#############
echo2 "Finished running ${PACKAGE_NAME} Degradome pipeline version $DEG_VERSION"
echo2 "---------------------------------------------------------------------------------"
touch .${GENOME}.DEG_VERSION.${DEG_VERSION}
