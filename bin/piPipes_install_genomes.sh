
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
export GENOME_INSTALL_VERSION=1.0.0

#########
# USAGE #
#########
usage () {
cat << EOF

Genome install pipeline v$GENOME_INSTALL_VERSION from $BOLD$PACKAGE_NAME$RESET
$INSTALL_USAGE${RESET}
Please email $CONTACT_EMAILS for any questions or bugs. 
Thank you for using it. 

${UNDERLINE}usage${RESET}:
	piPipes install \ 
		-g [dm3|mm9|...] \  
		-l http://www.link.to.iGenome.tar.gz \ 
		-D 

OPTIONS:
	-h      Show this message
	-v      Print out version
${REQUIRED}[ required ]
	-g      Name of the genome to install
${OPTIONAL}[ optional ]
	-l      Link to the iGenome, this is required if your genome is not in piPipes
	-D      Only do downloading but not other computation
		 This is designed be used when the user wants to separate downloading and other works. For instance, only the
		 head node on a hpcc has internet access but it is not appropriate to be used to build index.
	-c      Number of CPU to use
EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hg:c:l:vD" OPTION; do
	case $OPTION in
		h)	usage && exit 1 ;;
		g)	export GENOME=$OPTARG  ;;
		c)  CPU=$OPTARG  ;;
		v)	echo2 "GENOME_INSTALL_VERSION: v$GENOME_INSTALL_VERSION" && exit 0 ;;
		l)	LINK=$OPTARG  ;;
		D)	DOWNLOAD_ONLY=1 ;;
		*)	usage && exit 1 ;;
	esac
done
[[ -z $GENOME ]] && usage && echo2 "Missing option -g for version of genome assembly to install" "error"
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
mkdir -p $PIPELINE_DIRECTORY/common/$GENOME || echo2 "Cannot create directory $PIPELINE_DIRECTORY/common/$GENOME... Exiting..." "error"
cd $PIPELINE_DIRECTORY/common/$GENOME || echo2 "Cannot access directory $PIPELINE_DIRECTORY/common/$GENOME... Exiting..." "error"

#################################
# creating output files/folders #
#################################
mkdir -p $PIPELINE_DIRECTORY/Rlib
mkdir -p STARIndex
mkdir -p mrFastIndex

########################
# running binary check #
########################
checkBin "wget"
checkBin "samtools"
checkBin "faSize"
checkBin "bowtie-build"
checkBin "bowtie2-build"
checkBin "bwa"
checkBin "bowtie"
checkBin "gtfToGenePred"
checkBin "genePredToBed"
checkBin "bedtools_piPipes"
checkBin "mrfast"
checkBin "Rscript"

######################################
# reading the genome record of piPipes #
######################################
echo2 "Reading iGenome URL"
. $PIPELINE_DIRECTORY/common/iGenome_URL.txt
[ -z "${!GENOME}" -a -z "$LINK" ] && echo2 "It appeared that $GENOME is not in piPipes's record. Please provide the link to download the iGenome file with -l options" "error"
echo2 "$GENOME is in the record!"

##############################################
# instructions on getting repBase annotation #
##############################################
[ ! -f "$PIPELINE_DIRECTORY/common/$GENOME/${GENOME}.repBase.fa" ] && \
	echo2 "It appears that piPipes does not have the transposon consensus sequence ready for this genome.\nPlease provide this information in fasta format and name it \n\n\t$PIPELINE_DIRECTORY/common/$GENOME/${GENOME}.repBase.fa\n\n*For your convinience, the fasta of all repBase record can be found under the \"common\" folder.\nDo you still wish to proceed? [y/n]" "warning" && \
	read PROCEED && \
	case $PROCEED in 
		N|n|no) exit ;;
		Y|y|yes) ;;
		*) echo2 "unreognized answer" "error" ;;
	esac

#################################################
# to define the length range of siRNA and piRNA #
#################################################
echo2 "Reading length definition from the user"
# reading variables from the user
echo2 "Please input parameters for small RNA seq \n( If you have installed this genome before and don't want to change the parameters, please press ENTER to skip )"
echo2 "How many mismatches should be allowed for rRNA mapping by bowtie?"
read rRNA_MM
[ ! -z $rRNA_MM ] && echo "export rRNA_MM=$rRNA_MM" > $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "How many mismatches should be allowed for microRNA hairping mapping by bowtie?"
read hairpin_MM 
[ ! -z $hairpin_MM ] && echo "export hairpin_MM=$hairpin_MM" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "How many mismatches should be allowed for genome mapping by bowtie?"
read genome_MM 
[ ! -z $genome_MM ] && echo "export genome_MM=$genome_MM" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "How many mismatches should be allowed for trasnposons/piRNAcluster mapping by bowtie?"
read transposon_MM 
[ ! -z $transposon_MM ] && echo "export transposon_MM=$transposon_MM" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "What is the shortest length for siRNA?"
read siRNA_bot 
[ ! -z $siRNA_bot ] && echo "export siRNA_bot=$siRNA_bot" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "What is the longest length for siRNA?"
read siRNA_top 
[ ! -z $siRNA_top ] && echo "export siRNA_top=$siRNA_top" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "What is the shortest length for piRNA?"
read piRNA_bot 
[ ! -z $piRNA_bot ] && echo "export piRNA_bot=$piRNA_bot" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "What is the longest length for piRNA?"
read piRNA_top 
[ ! -z $piRNA_top ] && echo "export piRNA_top=$piRNA_top" >> $PIPELINE_DIRECTORY/common/$GENOME/variables
echo2 "Done. If you would like to change the variables, please edit file: $PIPELINE_DIRECTORY/common/$GENOME/variables manually" "warning"
echo2 "------------------------------------------"

##############################
# beginning running pipeline #
##############################
echo2 "Begining installing the genome $GENOME"
eval LINK='$'`echo $GENOME`

#######################################
# install R packages if not available #
#######################################
echo2 "Testing/Installing missing R packages"
[ ! -f $PIPELINE_DIRECTORY/common/.R_pkg_installed ] && Rscript $PIPELINE_DIRECTORY/bin/piPipes_install_packages.R 1>&2 && touch $PIPELINE_DIRECTORY/common/.R_pkg_installed

#################
# Download only #
#################
if [ "$DOWNLOAD_ONLY" == "1" ] ; then \
	echo2 "Downloading iGenome in Download ONLY mode."
	rm -rf $IGENOME_TAR_NAME; $wget "$LINK" || echo2 "Failed to download the genome file, please check the internet." "error"
	[ "$GENOME" == "dm3" ] && [ ! -s chrU.fa ] && echo2 "Downloading chrU for dm3" && wget "ftp://hgdownload.cse.ucsc.edu/goldenPath/dm3/chromosomes/chrU.fa.gz"
	[ ! -s UCSC.RepeatMask.bed -a ! -s UCSC.RepeatMask.bed.gz ] && \
		echo2 "Downloading repeatMasker files from UCSC." && \
		mkdir -p rmsk && cd rmsk && \
		rsync -a -P "rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/*rmsk.txt.gz" . && \
		zcat *gz | \
		awk 'BEGIN{FS=OFS="\t"}{print $6, $7-1, $8, $11, $1, $10}' > ../UCSC.RepeatMask.bed && \
		cd ..
	echo2 "Downloading finished but installation is not. Now you can run \"install\" again to install the rest without internet."
	exit 0;
fi

######################################
# iGenome download and uncompressing #
######################################
echo2 "Downloading iGenome $GENOME"
IGENOME_TAR_NAME=`basename $LINK`
IGENOME_DIR_NAME=${IGENOME_TAR_NAME%_UCSC*} # this only works for UCSC version of the iGenome; if you modify the code to use other assembly, please change here as well
[ ! -s $IGENOME_TAR_NAME ] && ( rm -rf $IGENOME_TAR_NAME && wget "$LINK" || echo2 "Failed to download the genome file, please check the internet." "error" )
echo2 "Uncompressing genome $GENOME"
[ ! -d $IGENOME_DIR_NAME -a -s $IGENOME_TAR_NAME ] && ( tar -zxvf $IGENOME_TAR_NAME || echo2 "Failed to unarchiving iGenome" "error" )

#############################
# links files & build index #
#############################
echo2 "Preparing genomic sequence/annotation and making indexes"
	
# patch for dm3; iGenome does not have chrU and X-TAS, hence add it manually. Indexes have to be rebuild as well. 
case $GENOME in	
dm3)	
	mkdir -p BowtieIndex
	mkdir -p Bowtie2Index
	mkdir -p BWAIndex
	echo2 "Applying patch to dm3: incorporating chrU and X-TAS"
	TAS='X-TAS.fa'
	[ ! -s $TAS ] && echo2 "Cannot file X-TAS.fa file, please reclone your git" "error"
	[ ! -s chrU.fa -a ! -s chrU.fa.gz ] && ( wget "ftp://hgdownload.cse.ucsc.edu/goldenPath/dm3/chromosomes/chrU.fa.gz" || echo2 "Failed to download chrU for dm3." "error" )
	[ ! -s chrU.fa -a -s chrU.fa.gz ] && ( gunzip chrU.fa.gz || echo2 "Failed to uncompressing chrU.fa.gz. Please re-download it from ftp://hgdownload.cse.ucsc.edu/goldenPath/dm3/chromosomes/chrU.fa.gz" "error" )
	[ ! -s ${GENOME}.fa ] && cat $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/WholeGenomeFasta/genome.fa chrU.fa $TAS > ${GENOME}.fa
	[ ! -s ${GENOME}.fa.fai ] && samtools faidx ${GENOME}.fa
	[ ! -s ${GENOME}.ChromInfo.txt ] && faSize -tab -detailed ${GENOME}.fa > ${GENOME}.ChromInfo.txt
# bowtie
	[ ! -s BowtieIndex/genome.1.ebwt ] && bowtie-build ${GENOME}.fa BowtieIndex/genome
	[ ! -s BowtieIndex/transposon.sizes ] && bowtie-build ${GENOME}.transposon.fa BowtieIndex/transposon && faSize -tab -detailed ${GENOME}.transposon.fa > BowtieIndex/transposon.sizes
# bowtie2
	[ ! -s Bowtie2Index/genome.1.bt2 ] && bowtie2-build ${GENOME}.fa Bowtie2Index/genome
	[ ! -s Bowtie2Index/transposon.sizes ] && bowtie2-build ${GENOME}.transposon.fa Bowtie2Index/transposon && faSize -tab -detailed ${GENOME}.transposon.fa > Bowtie2Index/transposon.sizes
# BWA
	[ ! -s BWAIndex/genome.bwt ] && bwa index -p BWAIndex/genome ${GENOME}.fa

;;
*)
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/WholeGenomeFasta/genome.fa ${GENOME}.fa
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/WholeGenomeFasta/genome.fa.fai ${GENOME}.fa.fai
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Annotation/Genes/ChromInfo.txt ${GENOME}.ChromInfo.txt
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/BowtieIndex
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/Bowtie2Index
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/BWAIndex
	ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Annotation/Genes/cytoBand.txt
;;
esac

# converting gtf to bed and extract the fasta
echo2 "Extracting sequence from iGenome gtf file"
ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Annotation/Genes/genes.gtf ${GENOME}.genes.gtf
[ ! -s ${GENOME}.genes.bed12 ] && gtfToGenePred ${GENOME}.genes.gtf ${GENOME}.genes.gp && genePredToBed ${GENOME}.genes.gp ${GENOME}.genes.bed12
[ ! -s ${GENOME}.genes.fa ] && bedtools_piPipes getfasta -fi ${GENOME}.fa -bed ${GENOME}.genes.bed12 -fo ${GENOME}.genes.fa -name -split -s

# STAR index for the genome
echo2 "Building STAR index for genome"
[ ! -s STARIndex/SAindex ] && \
STAR --runMode genomeGenerate --runThreadN $CPU --genomeDir STARIndex --genomeFastaFiles ${GENOME}.fa --sjdbGTFfile ${GENOME}.genes.gtf --sjdbOverhang 99 # TODO: this 99 is not really optimized...

# mrFast index for the genome
echo2 "Building mrFast index for genome"
[ ! -s mrFastIndex/${GENOME}.fa.index ] && \
	ln -st mrFastIndex/ ../${GENOME}.fa && \
	ln -st mrFastIndex/ ../${GENOME}.fa.fai && \
	mrfast --index mrFastIndex/${GENOME}.fa

# rRNA index
echo2 "Building Bowtie/Bowtie2 index for rRNA"
ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/AbundantSequences/*ibosomal.fa rRNA.fa
[ ! -s BowtieIndex/rRNA.sizes ] && bowtie-build rRNA.fa BowtieIndex/rRNA && faSize -tab -detailed rRNA.fa > BowtieIndex/rRNA.sizes
[ ! -s Bowtie2Index/rRNA.sizes ] && bowtie2-build rRNA.fa Bowtie2Index/rRNA && faSize -tab -detailed rRNA.fa > Bowtie2Index/rRNA.sizes

# microRNA and hairpin index
echo2 "Building index for microRNA hairpin"
[ ! -s ${GENOME}.hairpin.fa ] && awk '{if ($1~/^>/) print $1; else {gsub ("U","T", $0); print}}' $IGENOME_DIR_NAME/UCSC/$GENOME/Annotation/SmallRNA/precursor.fa > ${GENOME}.hairpin.fa
[ ! -s ${GENOME}.mature.fa ] &&  awk '{if ($1~/^>/) print $1; else {gsub ("U","T", $0); print}}' $IGENOME_DIR_NAME/UCSC/$GENOME/Annotation/SmallRNA/mature.fa > ${GENOME}.mature.fa
[ ! -s BowtieIndex/hairpin.sizes ] && bowtie-build ${GENOME}.hairpin.fa BowtieIndex/hairpin && faSize -tab -detailed ${GENOME}.hairpin.fa > BowtieIndex/hairpin.sizes
[ ! -s mature2hairpin.uniq.bed ]  && bowtie -S -f -v 0 -m 1 --best --strata --max ${GENOME}.mature.multiMapper.fa BowtieIndex/hairpin ${GENOME}.mature.fa 1> /dev/stdout 2> /dev/null | samtools view -uS - | bedtools_piPipes bamtobed -i - | awk '$6=="+"' > mature2hairpin.uniq.bed
[ ! -s mature2hairpin.multi.bed ] && bowtie -S -f -v 0 -a   --best --strata BowtieIndex/hairpin ${GENOME}.mature.multiMapper.fa 1> /dev/stdout 2> /dev/null | samtools view -uS - | bedtools_piPipes bamtobed -i - | awk '$6=="+"' > mature2hairpin.multi.bed
[ ! -s mature2hairpin.allMapper.bed ] && cat mature2hairpin.uniq.bed mature2hairpin.multi.bed > mature2hairpin.allMapper.bed

# repBase | transposon indexes # the pipiline should include the repBase.fa
echo2 "Building Bowtie/BWA index for repBase transposon annotation" 
[ ! -s ${GENOME}.repBase.fa ] && echo2 "Missing ${GENOME}.repBase.fa, if you are installing genomes other than dm3 or mm9, please retrieve that file from repBase" "error"
[ ! -s BowtieIndex/repBase.sizes ] && bowtie-build ${GENOME}.repBase.fa BowtieIndex/repBase && faSize -tab -detailed ${GENOME}.repBase.fa > BowtieIndex/repBase.sizes
[ ! -s ${GENOME}.repBase.eref ] && mkdir -p repBase && faSplit byname ${GENOME}.repBase.fa repBase/ && for i in repBase/*; do echo -e "`basename ${i%.fa}`\t`readlink -f $i`"; done > ${GENOME}.repBase.eref

# piRNA cluster indexes
echo2 "Building Bowtie/BWA index for piRNA cluster"
[ ! -s ${GENOME}.piRNAcluster.bed.gz ] && echo2 "Missing ${GENOME}.piRNAcluster.bed.gz, you are using a genome that is not optimized, some functions of the pipeline won't work." "warning"
[ ! -s ${GENOME}.piRNAcluster.fa ] && bedtools_piPipes getfasta -fi ${GENOME}.fa -bed ${GENOME}.piRNAcluster.bed.gz -fo ${GENOME}.piRNAcluster.fa -name -split -s
[ ! -s BowtieIndex/piRNAcluster.sizes ] && bowtie-build ${GENOME}.piRNAcluster.fa BowtieIndex/piRNAcluster && faSize -tab -detailed ${GENOME}.piRNAcluster.fa > BowtieIndex/piRNAcluster.sizes

# genes + repBase + cluster indexes
echo2 "Building Bowtie/BWA index for repBase + piRNA cluster + genes"
[ ! -s ${GENOME}.gene+cluster+repBase.fa ] && cat ${GENOME}.genes.fa  ${GENOME}.piRNAcluster.fa  ${GENOME}.repBase.fa  >  ${GENOME}.gene+cluster+repBase.fa
[ ! -s BowtieIndex/gene+cluster+repBase.sizes ] && bowtie-build ${GENOME}.gene+cluster+repBase.fa BowtieIndex/gene+cluster+repBase && faSize -tab -detailed ${GENOME}.gene+cluster+repBase.fa > BowtieIndex/gene+cluster+repBase.sizes
[ ! -s Bowtie2Index/gene+cluster+repBase.sizes ] && bowtie2-build ${GENOME}.gene+cluster+repBase.fa Bowtie2Index/gene+cluster+repBase && faSize -tab -detailed ${GENOME}.gene+cluster+repBase.fa > Bowtie2Index/gene+cluster+repBase.sizes

# unzipping the UCSC.RepeatMask.bed.gz shipped with the pipeline
# if the piPipes already has this in the github, just unzip it
[ ! -s UCSC.RepeatMask.bed -a -s UCSC.RepeatMask.bed.gz ] && gunzip UCSC.RepeatMask.bed.gz 
# if the piPipes does not have it, download it from UCSC
[ ! -s UCSC.RepeatMask.bed -a ! -s UCSC.RepeatMask.bed.gz ] && \
	mkdir -p rmsk && \
	cd rmsk && \
	rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/*rmsk.txt.gz . && \
	zcat *gz | \
	awk 'BEGIN{FS=OFS="\t"}{print $6, $7-1, $8, $11, $1, $10}' > ../UCSC.RepeatMask.bed && \
	cd ..

# making gtf files for htseq-count
echo2 "Making GTF file for HTSeq-count"
case $GENOME in	
dm3)
	echo2 "Nothing needs to be done for dm3: the gtf file for htseq-count is stored on github"
;;
mm9)
	echo2 "Merging repeat masker with prepachytene and pachytene clusters defined in the Zamore and ZLab"
	[ ! -s UCSC.RepeatMask.gtf ] && awk 'BEGIN{FS=OFS="\t"}{ if (!c[$4]) c[$4]=0; ++c[$4]; $4=$4"."c[$4]; print $0}' UCSC.RepeatMask.bed | bedToGenePred stdin /dev/stdout | genePredToGtf file stdin /dev/stdout | awk '$3=="exon"' > UCSC.RepeatMask.gtf
	[ ! -s Zamore.NM.gtf ] && zcat Zamore.NM.bed12.gz | bedToGenePred stdin /dev/stdout | genePredToGtf file stdin /dev/stdout | awk '$3=="exon"' > Zamore.NM.gtf
	[ ! -s Zamore.NR.gtf ] && zcat Zamore.NR.bed12.gz | bedToGenePred stdin /dev/stdout | genePredToGtf file stdin /dev/stdout | awk '$3=="exon"' > Zamore.NR.gtf
	[ ! -s ${GENOME}.genes+repBase+cluster.gtf ] && cat Zamore.NM.gtf Zamore.NR.gtf UCSC.RepeatMask.gtf > ${GENOME}.genes+repBase+cluster.gtf && rm -rf Zamore.NM.gtf Zamore.NR.gtf UCSC.RepeatMask.gtf
;;
*)	
	[ ! -s UCSC.RepeatMask.gtf ] && awk 'BEGIN{FS=OFS="\t"}{ if (!c[$4]) c[$4]=0; ++c[$4]; $4=$4"."c[$4]; print $0}' UCSC.RepeatMask.bed | bedToGenePred stdin /dev/stdout | genePredToGtf file stdin /dev/stdout | awk '$3=="exon"' > UCSC.RepeatMask.gtf
	[ ! -s ${GENOME}.piRNAcluster.gtf ] && zcat ${GENOME}.piRNAcluster.bed.gz | bedToGenePred stdin /dev/stdout | genePredToGtf file stdin /dev/stdout | awk '$3=="exon"' > ${GENOME}.piRNAcluster.gtf
	[ ! -s ${GENOME}.genes+repBase+cluster.htseq.gtf ] && cat ${GENOME}.genes.gtf  UCSC.RepeatMask.gtf  ${GENOME}.piRNAcluster.gtf  >  ${GENOME}.genes+repBase+cluster.htseq.gtf && rm -rf UCSC.RepeatMask.gtf  ${GENOME}.piRNAcluster.gtf
;;
esac

echo $GENOME >> $PIPELINE_DIRECTORY/common/genome_supported.txt





