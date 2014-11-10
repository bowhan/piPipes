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
export RNASEQ_VERSION=1.0.3

#########
# USAGE #
#########
usage () {
cat << EOF

RNASeq pipeline single library mode v$RNASEQ_VERSION from the $BOLD$PACKAGE_NAME$RESET
$RNASEQ_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

==================< paired-end >==================
${UNDERLINE}usage${RESET}:
	piPipes rna \ 
		-l left.fq[.gz] \ 
		-r right.fq[.gz] \ 
		-g dm3 \ 
		-o output_directory [current working directory] \ 
		-c cpu [8] \ 
		-B 38 [21]
	
==================< single-end >==================
${UNDERLINE}usage${RESET}:
	piPipes rna \ 
		-i input.fq[.gz] \ 
		-g dm3 \ 
		-o output_directory [current working directory] \ 
		-c cpu [8] \ 
		-B 38 [21]

By default, the pipeline assumes dUTR based method. For Paired-End sample: \2 is in the same 
direction as the transcript, opposite to ligation-based degradome/CAGE-seq. 
Use option -L to change the behavior.


OPTIONS:
	-h      Show this message
	-v      Print out the version
${REQUIRED}[ required ]
	-l      Left reads from Paired-End sequencing
	-r      Right reads from Paired-End sequencing
	-i      Input reads from Single-End sequence
	-g      Genome assembly name, like mm9 or dm3. required
	        Check $PIPELINE_DIRECTORY/common/genome_supported.txt for genome assemblies currently installed; 
	        Use "install" to install new genome
${OPTIONAL}[ optional ]
	-L      Ligation based library preperation method; Left reads (\1) being in the same direction as the transcripts. default: off (dUTR based, \2 reads being in the same direction)
	-o      Output directory, default: current directory $PWD
	-c      Number of CPUs to use, default: 8
	-B      How many rounds of batch algorithm to run for eXpress, default: 21
	-D      Delete large bed/bam files after pipeline finishes to save space (this step can also be ran separately), default: false
EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hl:r:i:c:o:g:B:vLD" OPTION; do
	case $OPTION in
		h)	usage && exit 0 ;;
		l)	LEFT_FASTQ=`readlink -f $OPTARG`;  PE_MODE=1 ;;
		r)	RIGHT_FASTQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
		i)	INPUT_FASTQ=`readlink -f $OPTARG`; SE_MODE=1 ;;
		o)	OUTDIR=`readlink -f $OPTARG` ;;
		c)	CPU=$OPTARG ;;
		g)	export GENOME=${OPTARG};;
		v)	echo2 "RNASEQ_VERSION: v$RNASEQ_VERSION" && exit 0 ;;
		L)	LIGATIONLIB=1 ;; # ligation based
		B)	eXpressBATCH=$OPTARG ;;
		D)	CLEAN=1;;
		*)	usage && exit 1 ;;
	esac
done

if [[ -z $PE_MODE && -z $SE_MODE ]]; then usage ; echo2 "Please specify the input file!" "error"; fi
if [[ -n $PE_MODE && -n $SE_MODE ]]; then usage ; echo2 "Please only choose single-end OR paired-end, but not both" "error"; fi

# if INPUT_FASTQ or GENOME is undefined, print out usage and exit
if [[ -n $PE_MODE ]]; then
	[[ -z $LEFT_FASTQ ]] && usage && echo2 "Missing option -l for input fastq of left file, or file does not exist " "error"
	[[ -z $RIGHT_FASTQ ]] && usage && echo2 "Missing option -r for input fastq of right file, or file does not exist " "error"
	[ ! -f $LEFT_FASTQ ] && echo2 "Cannot find input file $LEFT_FASTQ" "error"
	[ ! -f $RIGHT_FASTQ ] && echo2 "Cannot find input file $RIGHT_FASTQ" "error"
fi
if [[ -n $SE_MODE ]]; then
	[[ -z $INPUT_FASTQ ]] && usage && echo2 "Missing option -i for input fastq, or file does not exist " "error"
	[ ! -f $INPUT_FASTQ ] && echo2 "Cannot find input file $INPUT_FASTQ" "error"
fi
[[ -z $GENOME ]]  && usage && echo2 "Missing option -g for specifying which genome assembly to use" "error"

# check whether the this genome is supported or not
check_genome $GENOME
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "${eXpressBATCH##*[!0-9]*}" ] || eXpressBATCH=21
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
	EXPRESS_OPTION_PE="--fr-stranded"
	EXPRESS_OPTION_SE="--f-stranded"
else 
	LIBRARY_TYPE="fr-firststrand"
	END_TO_REVERSE_STRAND=1
	SENSE_HTSEQ_OPT="reverse"; 
	ANTISENSE_HTSEQ_OPT="yes"; 
	EXPRESS_OPTION_PE="--rf-stranded"
	EXPRESS_OPTION_SE="--r-stranded"
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
BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR
DIRECTMAPPING_DIR=direct_transcriptome_mapping && mkdir -p $DIRECTMAPPING_DIR
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

#############
# Variables #
#############
# step counter
STEP=1
# job uid
JOBUID=`echo ${INPUT_FASTQ}${LEFT_FASTQ}${RIGHT_FASTQ} | md5sum | cut -d" " -f1`
INPUT_FASTQ_NAME=`basename $INPUT_FASTQ`
LEFT_FASTQ_NAME=`basename $LEFT_FASTQ`
RIGHT_FASTQ_NAME=`basename $RIGHT_FASTQ`

if [[ -n $PE_MODE ]]; then
	PREFIX=`echo -e "${LEFT_FASTQ_NAME}\n${RIGHT_FASTQ_NAME}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
	[ -z "${PREFIX}" ] && export PREFIX=${LEFT_FASTQ_NAME%.f[aq]*} # if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
fi
[[ -n $SE_MODE ]] && PREFIX=${INPUT_FASTQ_NAME%.f[aq]*}

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
# bin size for graph
export BINSIZE=1000
##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] RNA-Seq pipeline version $RNASEQ_VERSION" 

if [[ -n $PE_MODE ]]; then
# Paired-End
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
		--limitOutSAMoneReadBytes 1000000 \
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
		samtools view -bS ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam && \
		rm -rf  ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam && \
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
		bedtools_piPipes bamtobed -bed12 -tag NH -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
		awk -v strand=$END_TO_REVERSE_STRAND 'BEGIN{FS=OFS="\t"}{e=substr($4,length($4)); if (e==strand) $6=($6=="+"?"-":"+"); if ($5==1) print $0; }' > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 && \
		para_file=${BW_OUTDIR}/${RANDOM}${RANDOM}.makeBigWigPE.para && \
		echo "bedtools_piPipes genomecov -scale $NormScale -split -bg -strand + -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bigWig"  >  $para_file && \
		echo "bedtools_piPipes genomecov -scale $NormScale -split -bg -strand - -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bigWig" >> $para_file && \
		ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
		rm -rf ${para_file} ${para_file}.completed ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph && \
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
	echo2 "Mapping to genes and transposon directly with Bowtie2"
	. $COMMON_FOLDER/genomic_features
	# for fly genome, the transcripts from piRNA cluster are usually undetectable. including them in eXpress will actually have negative influence.
	if [ "$GENOME" == "dm3" ]; then
		TRANSCRIPTOME_INDEX="gene+transposon"
	else
		TRANSCRIPTOME_INDEX="gene+cluster+repBase"
	fi
	TRANSCRIPTOME_SIZES=$BOWTIE2_INDEXES/${TRANSCRIPTOME_INDEX}.sizes

	[ ! -f .${JOBUID}.status.${STEP}.direct_mapping ] && \
	bowtie2 -x ${TRANSCRIPTOME_INDEX} \
		-1 ${xrRNA_LEFT_FQ} \
		-2 ${xrRNA_RIGHT_FQ} \
		-q \
		$bowtie2PhredOption \
		-a \
		-X 800 \
		--no-mixed \
		--quiet \
		-p $CPU \
		2> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.log | \
	samtools view -bS - > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam && \
	samtools sort -o -@ $CPU ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam ${DIRECTMAPPING_DIR}/foo | \
	bedtools_piPipes bamtobed -i - | \
	awk -v etr=$END_TO_REVERSE_STRAND -v MAPQ=10 'BEGIN{FS=OFS="\t"}{l=split($4,arr,""); if (arr[l]==etr) $6=($6=="+"?"-":"+"); if ($5 > MAPQ) print $0}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed && \
	touch .${JOBUID}.status.${STEP}.direct_mapping

	echo2 "Making summary graph"
	[ ! -f .${JOBUID}.status.${STEP}.make_direct_mapping_sum ] && \
	grep -v 'NM_' $TRANSCRIPTOME_SIZES | grep -v 'NR_' > ${DIRECTMAPPING_DIR}/transposon.sizes && \
	bedtools_piPipes genomecov -i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed -g $TRANSCRIPTOME_SIZES -strand + -bg \
		> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph && \
	bedtools_piPipes genomecov -i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed -g $TRANSCRIPTOME_SIZES -strand - -bg \
		| awk 'BEGIN{FS=OFS="\t"}{$4=-$4;print $0}' \
		> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph && \
	bedGraphToBigWig ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph  $TRANSCRIPTOME_SIZES ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig && \
	bedGraphToBigWig ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph $TRANSCRIPTOME_SIZES ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig && \
	rm -f ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph && \
	paraFile=${DIRECTMAPPING_DIR}/bigWigSummary.para && \
	bgP=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig && \
	bgM=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig && \
	awk -v bgP=${bgP} -v bgM=${bgM} -v binSize=${BINSIZE} '{print "bigWigSummary", bgP, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x27 >", bgP"."$1; print "bigWigSummary", bgM, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x27 >", bgM"."$1;}' ${DIRECTMAPPING_DIR}/transposon.sizes > $paraFile && \
	ParaFly -c $paraFile -CPU $CPU && \
	paraFile=${OUTDIR}/drawFigures && \
	rm -rf ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary && \
	for i in `cut -f1 ${DIRECTMAPPING_DIR}/transposon.sizes`; do \
		[ ! -s ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i  ] && awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i
		[ ! -s ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i ] && awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i
		awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i  > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t
		awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t
		paste ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t | cut -f1,2,3,6 >> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary
		rm -rf ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t
	done && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary $CPU $NormScale 1>&2 && \
	PDFs=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted*pdf && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.pdf ${PDFs} && \
	rm -rf ${PDFs} && \
	touch .${JOBUID}.status.${STEP}.make_direct_mapping_sum

	[ ! -f .${JOBUID}.status.${STEP}.eXpress_quantification ] && \
	express \
		$EXPRESS_OPTION_PE \
		-B $eXpressBATCH \
		-o $DIRECTMAPPING_DIR \
		--no-update-check \
		--library-size ${MapMass%.*} \
		$COMMON_FOLDER/${GENOME}.${TRANSCRIPTOME_INDEX}.fa \
		${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam \
		1>&2 2> $DIRECTMAPPING_DIR/${PREFIX}.${TRANSCRIPTOME_INDEX}.eXpress.log && \
	awk -v depth=$NormScale 'BEGIN{OFS="\t"; getline; print}{$8*=depth; print}' $DIRECTMAPPING_DIR/results.xprs > \
		$DIRECTMAPPING_DIR/results.xprs.normalized && \
	touch .${JOBUID}.status.${STEP}.eXpress_quantification
	STEP=$((STEP+1))

	#############
	# finishing #
	#############
	if [[ "$CLEAN" == 1 ]]; then
		rm -f $GENOMIC_MAPPING_DIR/*mate1 $GENOMIC_MAPPING_DIR/*mate2
		rm -f ${READS_DIR}/${PREFIX}.x_rRNA.1.fq ${READS_DIR}/${PREFIX}.x_rRNA.2.fq
		rm -f ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed
	fi
else # Single-End
	
	###########################
	# determine fastQ version #
	###########################
	echo2 "Determining the version of fastQ using SolexaQA"
	# determine version of fastq used, using a modified SolexaQA.pl
	PHRED_SCORE=`perl $PIPELINE_DIRECTORY/bin/SolexaQA_piPipes.pl ${INPUT_FASTQ}`
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
		-U $INPUT_FASTQ \
		-q \
		$bowtie2PhredOption \
		--very-fast \
		-k 1 \
		--un ${READS_DIR}/${PREFIX}.x_rRNA.fq \
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
	xrRNA_FQ=${READS_DIR}/${PREFIX}.x_rRNA.fq && \
	[ ! -f .${JOBUID}.status.${STEP}.genome_mapping ] && \
	STAR \
		--runMode alignReads \
		--limitOutSAMoneReadBytes 1000000 \
		--genomeDir $STARINDEX \
		--readFilesIn ${xrRNA_FQ} \
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
		samtools view -bS ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam && \
		samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted && \
		samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam && \
		rm -rf  ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.bam && \
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
	# in order to make bigWig for unique mappers
	echo2 "Making bigWig from sorted bam"
	[ ! -f .${JOBUID}.status.${STEP}.make_bigWig ] && \
	if [[ "${LIGATIONLIB}" == 1 ]]; then
		bedtools_piPipes bamtobed -bed12 -tag NH -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
		awk 'BEGIN{OFS="\t"}{if ($5==1) print $0; }' > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12
	else
		bedtools_piPipes bamtobed -bed12 -tag NH -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam | \
		awk 'BEGIN{OFS="\t"}{ $6=($6=="+"?"-":"+"); if ($5==1) print $0;}' > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12
	fi && \
	para_file=${BW_OUTDIR}/${RANDOM}${RANDOM}.makeBigWig.para && \
	echo "bedtools_piPipes genomecov -scale $NormScale -split -bg -strand + -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bigWig"  >  $para_file && \
	echo "bedtools_piPipes genomecov -scale $NormScale -split -bg -strand - -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph && bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bigWig" >> $para_file && \
	ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	rm -rf ${para_file} ${para_file}.completed  ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Watson.bedGraph ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.Crick.bedGraph && \
	touch .${JOBUID}.status.${STEP}.make_bigWig
	STEP=$((STEP+1))

	#############################################
	# genomic feature counting with htSeq-count #
	#############################################
	# echo2 "Quantifying genomic features from genomic mapping using HTSeq-count"
	# [ ! -f .${JOBUID}.status.${STEP}.htseq_count ] && \
	# 	. $COMMON_FOLDER/genomic_features && \
	# 	[ ! -z $HTSEQ_TARGETS ] && \
	# 	para_file=$HTSEQ_DIR/${RANDOM}${RANDOM}.para && \
	# 	for t in "${HTSEQ_TARGETS[@]}"; do \
	# 		echo "htseq-count -m intersection-strict -s $SENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.strict.S.out" >> $para_file
	# 		echo "htseq-count -m intersection-strict -s $ANTISENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.strict.AS.out" >> $para_file
	# 		echo "htseq-count -m union -s $SENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.union.S.out" >> $para_file
	# 		echo "htseq-count -m union -s $ANTISENSE_HTSEQ_OPT -t exon -i transcript_id -q ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam ${!t} | awk 'BEGIN{FS=OFS=\"\t\"}{n=split (\$1,a,\".\"); ct[a[1]]+=\$2; }END{for (x in ct) {print x, ct[x]}}' | sort -k1,1 > ${HTSEQ_DIR}/${PREFIX}.x_rRNA.${GENOME}.${t}.htseqcount.union.AS.out" >> $para_file
	# 	done && \
	# ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
	# rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.Aligned.out.sam && \
	# touch .${JOBUID}.status.${STEP}.htseq_count
	# STEP=$((STEP+1))

	##################################################
	# direct mapping and quantification with eXpress #
	##################################################
	echo2 "Mapping to genes and transposon directly with Bowtie2"
	. $COMMON_FOLDER/genomic_features
	# for fly genome, the transcripts from piRNA cluster are usually undetectable. including them in eXpress will actually have negative influence.
	if [ "$GENOME" == "dm3" ]; then
		TRANSCRIPTOME_INDEX="gene+transposon"
	else
		TRANSCRIPTOME_INDEX="gene+cluster+repBase"
	fi
	TRANSCRIPTOME_SIZES=$BOWTIE2_INDEXES/${TRANSCRIPTOME_INDEX}.sizes

	[ ! -f .${JOBUID}.status.${STEP}.direct_mapping ] && \
	bowtie2 -x ${TRANSCRIPTOME_INDEX} \
		-U ${xrRNA_FQ} \
		-q \
		$bowtie2PhredOption \
		-a \
		--quiet \
		-p $CPU \
		2> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.log | \
	samtools view -bS - > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam && \
	if [[ "${LIGATIONLIB}" == 1 ]]; then
		samtools sort -o -@ $CPU ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam ${DIRECTMAPPING_DIR}/foo | \
		bedtools_piPipes bamtobed -i - | \
		awk -v MAPQ=10 'BEGIN{FS=OFS="\t"}{ if ($5 > MAPQ) print $0 }' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed && \
		touch .${JOBUID}.status.${STEP}.direct_mapping
	else
		samtools sort -o -@ $CPU ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam ${DIRECTMAPPING_DIR}/foo | \
		bedtools_piPipes bamtobed -i - | \
		awk -v MAPQ=10 'BEGIN{FS=OFS="\t"}{ $6=($6=="+"?"-":"+"); if ($5 > MAPQ) print $0}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed && \
		touch .${JOBUID}.status.${STEP}.direct_mapping
	fi

	echo2 "Making summary graph"
	[ ! -f .${JOBUID}.status.${STEP}.make_direct_mapping_sum ] && \
	grep -v 'NM_' $TRANSCRIPTOME_SIZES | grep -v 'NR_' > ${DIRECTMAPPING_DIR}/transposon.sizes && \
	bedtools_piPipes genomecov -i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed -g $TRANSCRIPTOME_SIZES -strand + -bg \
		> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph && \
	bedtools_piPipes genomecov -i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed -g $TRANSCRIPTOME_SIZES -strand - -bg \
		| awk 'BEGIN{FS=OFS="\t"}{$4=-$4;print $0}' \
		> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph && \
	bedGraphToBigWig ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph  $TRANSCRIPTOME_SIZES ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig && \
	bedGraphToBigWig ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph $TRANSCRIPTOME_SIZES ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig && \
	rm -f ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bedGraph ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bedGraph && \
	paraFile=${DIRECTMAPPING_DIR}/bigWigSummary.para && \
	bgP=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig && \
	bgM=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig && \
	awk -v bgP=${bgP} -v bgM=${bgM} -v binSize=${BINSIZE} '{print "bigWigSummary", bgP, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x27 >", bgP"."$1; print "bigWigSummary", bgM, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x27 >", bgM"."$1;}' ${DIRECTMAPPING_DIR}/transposon.sizes > $paraFile && \
	ParaFly -c $paraFile -CPU $CPU && \
	paraFile=${OUTDIR}/drawFigures && \
	rm -rf ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary && \
	for i in `cut -f1 ${DIRECTMAPPING_DIR}/transposon.sizes`; do \
		[ ! -s ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i  ] && awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i
		[ ! -s ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i ] && awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i
		awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i  > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t
		awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t
		paste ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t | cut -f1,2,3,6 >> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary
		rm -rf ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.plus.bigWig.$i.t ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.minus.bigWig.$i.t
	done && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bigWig.summary $CPU $NormScale 1>&2 && \
	PDFs=${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted*pdf && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.pdf ${PDFs} && \
	rm -rf ${PDFs} && \
	touch .${JOBUID}.status.${STEP}.make_direct_mapping_sum

	# we currently don't specify the direction for eXpress to allow counting of antisense transcripts;
	[ ! -f .${JOBUID}.status.${STEP}.eXpress_quantification ] && \
	express \
		$EXPRESS_OPTION_SE \
		-B $eXpressBATCH \
		-o $DIRECTMAPPING_DIR \
		--no-update-check \
		--library-size ${MapMass%.*} \
		$COMMON_FOLDER/${GENOME}.${TRANSCRIPTOME_INDEX}.fa \
		${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.bam \
		1>&2 2> $DIRECTMAPPING_DIR/${PREFIX}.${TRANSCRIPTOME_INDEX}.eXpress.log && \
	awk -v depth=$NormScale 'BEGIN{OFS="\t"; getline; print}{$8*=depth; print}' $DIRECTMAPPING_DIR/results.xprs > \
		$DIRECTMAPPING_DIR/results.xprs.normalized && \
	touch .${JOBUID}.status.${STEP}.eXpress_quantification
	STEP=$((STEP+1))

	#############
	# finishing #
	#############
	if [[ "$CLEAN" == 1 ]]; then
		rm -f $BW_OUTDIR/*bedGraph
		rm -f $GENOMIC_MAPPING_DIR/*bedpe
		rm -f $GENOMIC_MAPPING_DIR/*mate1 $GENOMIC_MAPPING_DIR/*mate2
		rm -f ${READS_DIR}/${PREFIX}.x_rRNA.1.fq ${READS_DIR}/${PREFIX}.x_rRNA.2.fq
		rm -f ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_INDEX}.sorted.unique.bed
	fi
fi # SE or PE

echo2 "Finished running ${PACKAGE_NAME} RNA-Seq pipeline version $RNASEQ_VERSION"
echo2 "---------------------------------------------------------------------------------"
touch .${GENOME}.${LIGATIONLIB}.RNASEQ_VERSION.${RNASEQ_VERSION}
