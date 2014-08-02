
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
This script is only used to install the D. melanogaster release 6 (BDGP6), which has
not been included in iGenome yet.
EOF
echo -e "${COLOR_END}"
}

export GENOME=dm6
CPU=$1
DOWNLOAD_ONLY=$2
	
mkdir -p $PIPELINE_DIRECTORY/common/$GENOME || echo2 "Cannot create directory $PIPELINE_DIRECTORY/common/$GENOME... Exiting..." "error"
cd $PIPELINE_DIRECTORY/common/$GENOME || echo2 "Cannot access directory $PIPELINE_DIRECTORY/common/$GENOME... Exiting..." "error"

#################################
# creating output files/folders #
#################################
mkdir -p $PIPELINE_DIRECTORY/Rlib
mkdir -p STARIndex
mkdir -p mrFastIndex

# download function for dm6, 6.01
function download_dm6 {
	DM6_GENOME_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-chromosome-r6.01.fasta.gz'
	DM6_GENE_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-gene-r6.01.fasta.gz'
	DM6_TRANSCRIPTS_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-transcript-r6.01.fasta.gz'
	# DM6_MIRNA_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-miRNA-r6.01.fasta.gz'
	DM6_TRNA_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-tRNA-r6.01.fasta.gz'
	DM6_NCRNA_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-ncRNA-r6.01.fasta.gz'
	DM6_TRANSPOSON_FA='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/fasta/dmel-all-transposon-r6.01.fasta.gz'
	DM6_GTF='ftp://ftp.flybase.net/releases/FB2014_04/dmel_r6.01/gtf/dmel-all-r6.01.gtf.gz'
	declare -a DOWNLOAD_=("DM6_GENOME_FA" "DM6_GENE_FA" "DM6_TRANSCRIPTS_FA" "DM6_TRNA_FA" "DM6_NCRNA_FA" "DM6_TRANSPOSON_FA" "DM6_GTF")
	for i in "${DOWNLOAD_[@]}"; do
		[ ! -f $i ] && wget --continue -O - ${!i} | gunzip > $i 
	done
}
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
echo2 "Downloading dm6 files from flyBase"
download_dm6
if [ "$DOWNLOAD_ONLY" == "1" ] ; then 
	echo2 "Finish downloading dm6 files from flyBase. Exiting..."
	exit; 
fi

#############################
# links files & build index #
#############################
echo2 "Preparing genomic sequence/annotation and making indexes"

mkdir -p BowtieIndex
mkdir -p Bowtie2Index
mkdir -p BWAIndex
[ ! -s ${GENOME}.fa ] && ln -s DM6_GENOME_FA ${GENOME}.fa
[ ! -s ${GENOME}.fa.fai ] && samtools faidx ${GENOME}.fa
[ ! -s ${GENOME}.ChromInfo.txt ] && faSize -tab -detailed ${GENOME}.fa > ${GENOME}.ChromInfo.txt
[ ! -s ${GENOME}.transposon.fa ] && ln -s DM6_TRANSPOSON_FA ${GENOME}.transposon.fa

# bowtie
[ ! -s BowtieIndex/genome.1.ebwt ] && bowtie-build ${GENOME}.fa BowtieIndex/genome
[ ! -s BowtieIndex/transposon.sizes ] && bowtie-build ${GENOME}.transposon.fa BowtieIndex/transposon && faSize -tab -detailed ${GENOME}.transposon.fa > BowtieIndex/transposon.sizes
# bowtie2
[ ! -s Bowtie2Index/genome.1.bt2 ] && bowtie2-build ${GENOME}.fa Bowtie2Index/genome
[ ! -s Bowtie2Index/transposon.sizes ] && bowtie2-build ${GENOME}.transposon.fa Bowtie2Index/transposon && faSize -tab -detailed ${GENOME}.transposon.fa > Bowtie2Index/transposon.sizes
# BWA
[ ! -s BWAIndex/genome.bwt ] && bwa index -p BWAIndex/genome ${GENOME}.fa

# converting gtf to bed and extract the fasta
# TODO: correct this after flybase correct the mdg4 error
[ ! -s ${GENOME}.genes.gtf ] && grep -v mdg4 DM6_GTF > ${GENOME}.genes.gtf # ln -s DM6_GTF ${GENOME}.genes.gtf
[ ! -s ${GENOME}.genes.bed12 ] && gtfToGenePred ${GENOME}.genes.gtf ${GENOME}.genes.gp && genePredToBed ${GENOME}.genes.gp ${GENOME}.genes.bed12
[ ! -s ${GENOME}.genes.fa ] && bedtools_piPipes getfasta -fi ${GENOME}.fa -bed ${GENOME}.genes.bed12 -fo ${GENOME}.genes.fa -name -split -s

# STAR index for the genome
echo2 "Building STAR index for genome"
[ ! -s STARIndex/SAindex ] && \
STAR --runMode genomeGenerate --runThreadN $CPU --genomeDir STARIndex --genomeFastaFiles ${GENOME}.fa --sjdbGTFfile ${GENOME}.genes.gtf --sjdbOverhang 99 # TODO: this 99 is not really optimized...

# mrFast index for the genome
echo2 "Building mrFast index for genome"
# TODO: seems that mrfast does not like the fasta...
if [ ! -s mrFastIndex/${GENOME}.fa.index ]; then
	[ ! -s ${GENOME}.fa ] && ln -st mrFastIndex/ ../${GENOME}.fa 
	[ ! -s ${GENOME}.fa.fai ] && ln -st mrFastIndex/ ../${GENOME}.fa.fai
	mrfast --index mrFastIndex/${GENOME}.fa
fi

# rRNA index
echo2 "Building Bowtie/Bowtie2 index for rRNA"
# ln -s $IGENOME_DIR_NAME/UCSC/$GENOME/Sequence/AbundantSequences/*ibosomal.fa rRNA.fa 2>/dev/null
# rRNA.fa has been included in the Github
[ ! -s BowtieIndex/rRNA.sizes ] && bowtie-build rRNA.fa BowtieIndex/rRNA && faSize -tab -detailed rRNA.fa > BowtieIndex/rRNA.sizes
[ ! -s Bowtie2Index/rRNA.sizes ] && bowtie2-build rRNA.fa Bowtie2Index/rRNA && faSize -tab -detailed rRNA.fa > Bowtie2Index/rRNA.sizes

# microRNA and hairpin index
echo2 "Building index for microRNA hairpin, still using fasta from mirBase"
[ ! -s BowtieIndex/hairpin.sizes ] && bowtie-build ${GENOME}.hairpin.fa BowtieIndex/hairpin && faSize -tab -detailed ${GENOME}.hairpin.fa > BowtieIndex/hairpin.sizes
[ ! -s mature2hairpin.uniq.bed ]  && bowtie -S -f -v 0 -m 1 --best --strata --max ${GENOME}.mature.multiMapper.fa BowtieIndex/hairpin ${GENOME}.mature.fa 1> /dev/stdout 2> /dev/null | samtools view -uS - | bedtools_piPipes bamtobed -i - | awk '$6=="+"' > mature2hairpin.uniq.bed
[ ! -s mature2hairpin.multi.bed ] && bowtie -S -f -v 0 -a   --best --strata BowtieIndex/hairpin ${GENOME}.mature.multiMapper.fa 1> /dev/stdout 2> /dev/null | samtools view -uS - | bedtools_piPipes bamtobed -i - | awk '$6=="+"' > mature2hairpin.multi.bed
[ ! -s mature2hairpin.allMapper.bed ] && cat mature2hairpin.uniq.bed mature2hairpin.multi.bed > mature2hairpin.allMapper.bed

# repBase | transposon indexes # the pipiline should include the repBase.fa
echo2 "Building Bowtie/BWA index for repBase transposon annotation" 
[ ! -s BowtieIndex/repBase.sizes ] && bowtie-build ${GENOME}.repBase.fa BowtieIndex/repBase && faSize -tab -detailed ${GENOME}.repBase.fa > BowtieIndex/repBase.sizes
[ ! -s ${GENOME}.repBase.eref ] && mkdir -p repBase && faSplit byname ${GENOME}.repBase.fa repBase/ && for i in repBase/*; do echo -e "`basename ${i%.fa}`\t`readlink -f $i`"; done > ${GENOME}.repBase.eref

# piRNA cluster indexes
echo2 "Building Bowtie/BWA index for piRNA cluster"
# TODO: annotate new piRNA cluster
[ ! -s ${GENOME}.piRNAcluster.bed ] && ln -s Brennecke.piRNAcluster.bed6 ${GENOME}.piRNAcluster.bed
[ ! -s ${GENOME}.piRNAcluster.fa ] && bedtools_piPipes getfasta -fi ${GENOME}.fa -bed ${GENOME}.piRNAcluster.bed -fo ${GENOME}.piRNAcluster.fa -name -split -s
[ ! -s BowtieIndex/piRNAcluster.sizes ] && bowtie-build ${GENOME}.piRNAcluster.fa BowtieIndex/piRNAcluster && faSize -tab -detailed ${GENOME}.piRNAcluster.fa > BowtieIndex/piRNAcluster.sizes

[ ! -s ${GENOME}.transposon.fa ] && ln -s DM6_TRANSPOSON_FA ${GENOME}.transposon.fa
[ ! -s ${GENOME}.gene+transposon.fa ] && cat ${GENOME}.genes.fa  ${GENOME}.transposon.fa  >  ${GENOME}.gene+transposon.fa
[ ! -s BowtieIndex/gene+transposon.sizes ] && bowtie-build		${GENOME}.gene+transposon.fa BowtieIndex/gene+transposon  && faSize -tab -detailed ${GENOME}.gene+transposon.fa > BowtieIndex/gene+transposon.sizes
[ ! -s Bowtie2Index/gene+transposon.sizes ] && bowtie2-build	${GENOME}.gene+transposon.fa Bowtie2Index/gene+transposon && faSize -tab -detailed ${GENOME}.gene+transposon.fa > Bowtie2Index/gene+transposon.sizes

# genes + repBase + cluster indexes
echo2 "Building Bowtie/BWA index for repBase + piRNA cluster + genes"
[ ! -s ${GENOME}.gene+cluster+repBase.fa ] && cat ${GENOME}.genes.fa  ${GENOME}.piRNAcluster.fa  ${GENOME}.repBase.fa  >  ${GENOME}.gene+cluster+repBase.fa
[ ! -s BowtieIndex/gene+cluster+repBase.sizes ] && bowtie-build ${GENOME}.gene+cluster+repBase.fa BowtieIndex/gene+cluster+repBase && faSize -tab -detailed ${GENOME}.gene+cluster+repBase.fa > BowtieIndex/gene+cluster+repBase.sizes
[ ! -s Bowtie2Index/gene+cluster+repBase.sizes ] && bowtie2-build ${GENOME}.gene+cluster+repBase.fa Bowtie2Index/gene+cluster+repBase && faSize -tab -detailed ${GENOME}.gene+cluster+repBase.fa > Bowtie2Index/gene+cluster+repBase.sizes

echo $GENOME >> $PIPELINE_DIRECTORY/common/genome_supported.txt





