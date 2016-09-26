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
declare -x  GenomeInstallModuleVesion=2.0.0
declare -xi Threads=8
declare -x  Genome=""
declare -a RequiredPrgrams=( \
	wget rsync mysql \
	faSize twoBitToFa gtfToGenePred genePredToBed \
	bowtie-build bowtie2-build bwa \
	samtools bedtools_piPipes Rscript \
	)
declare -a RequiredAnnotationFiles=( \
	rRNA.fa \
	hairpin.fa \
	mature.fa \
	transposon.fa \
	piRNAcluster.bed \
	)

#########
# USAGE #
#########
usage () {
cat << EOF

Genome install pipeline v$GenomeInstallModuleVesion from $BOLD$PACKAGE_NAME$RESET
$INSTALL_USAGE${RESET}
Please email $CONTACT_EMAILS for any questions or bugs.
Thank you for using it.

From version 2.0, ${PACKAGE_NAME} no longer obtain sequence and annotation from iGenome.
Instead, genome sequence and gene annotations are obtained from UCSC using either rsync
or mysql. However, a few annotation files need to be provided manually, including:
rRNA.fa		rRNA sequences in fasta, you will have to obtain the rRNA sequences.
mature.fa and hairpin.fa	mature and hairpin miRNA sequences in fasta, you can extract their sequences using "${MYBIN}/piPipes_extract_organiam_from_fa.py" from mature.fa and hairpin.fa files from mirBase.

Some of those files are included in the iGenome annotation, you can download them
manually and rename the file.

${UNDERLINE}usage${RESET}:
	piPipes install \
		-g [dm6|mm10|...] \
		-D

OPTIONS:
	-h      Show this message
	-v      Print out version
${REQUIRED}[ required ]
	-g      Name of the Genome assembly to install, choose from ftp://hgdownload.cse.ucsc.edu/goldenPath/
${OPTIONAL}[ optional ]
	-c      Number of Threads to use

EOF
echo -e "${COLOR_END}"
}

#############################
# ARGS reading and checking #
#############################
while getopts "hg:c:vDC" OPTION; do
	case $OPTION in
		h)	usage && exit 1 ;;
		v)	echo2 "GenomeInstallModuleVesion: v$GenomeInstallModuleVesion" && exit 0 ;;
		g)	Genome=$OPTARG  ;;
		c)  Threads=$OPTARG ;;
		*)	usage && exit 1 ;;
	esac
done

if [[ -z $Genome ]]; then 
	usage
	echo2 "Missing option -g for version of Genome assembly to install" error
fi 

########################
# running binary check #
########################
for program in "${RequiredPrgrams[@]}"; do assertBinExists $program; done

#################################
# creating output files/folders #
#################################
mkdir -p $PIPELINE_DIRECTORY/common/$Genome || echo2 "Cannot create directory $PIPELINE_DIRECTORY/common/$Genome... Exiting..." error
cd $PIPELINE_DIRECTORY/common/$Genome || echo2 "Cannot access directory $PIPELINE_DIRECTORY/common/$Genome... Exiting..." error
mkdir -p $PIPELINE_DIRECTORY/Rlib STARIndex BowtieIndex Bowtie2Index BWAIndex log

#######################################
# install R packages if not available #
#######################################
echo2 "Testing/Installing missing R packages"
Rscript $PIPELINE_DIRECTORY/bin/piPipes_install_packages.R &> log/install_R.log

######################################
# checking required annotation files #
######################################
for f in "${RequiredAnnotationFiles[@]}"; do
	if [[ ! -f ${f} ]]; then
		echo2 "Missing required file ${f}; Please either provide the correct file under ${PWD}/${f} or generate a fake ${f} file with a dummy sequence such as AAAAAAAAAAAAAAA" error 
	fi
done

############################
# read variables from user #
############################
if [[ ! -s $PIPELINE_DIRECTORY/common/$Genome/variables ]]; then
	echo2 "Reading length definition from the user"
	# reading variables from the user
	ReadUserParameter() {
		local message="${1}"
		local varname=${2}
		echo2 ${message}
		read ${varname}
		[[ ! -z ${!varname} ]] \
		&& echo "export ${varname}=${!varname}" >> $PIPELINE_DIRECTORY/common/$Genome/variables \
		|| echo2 "Invalid input" error
	}

	ReadUserParameter "How many mismatches should be allowed for rRNA mapping by bowtie?" rRNA_MM
	ReadUserParameter "How many mismatches should be allowed for microRNA hairping mapping by bowtie?" hairpin_MM
	ReadUserParameter "How many mismatches should be allowed for Genome mapping by bowtie?" Genome_MM
	ReadUserParameter "How many mismatches should be allowed for trasnposons/piRNAcluster mapping by bowtie?" transposon_MM
	ReadUserParameter "What is the shortest length for siRNA?" siRNA_bot
	ReadUserParameter "What is the longest length for siRNA?" siRNA_top
	ReadUserParameter "What is the shortest length for piRNA?" piRNA_bot
	ReadUserParameter "What is the longest length for piRNA?" piRNA_top

	echo2 "Done. If you would like to change the variables, please edit file: $PIPELINE_DIRECTORY/common/$Genome/variables manually"
fi

##############################
# beginning running pipeline #
##############################
echo2 "Beginning to install the Genome $Genome"

###################
# Genome sequence #
###################
if [[ ! -s ${Genome}.fa ]]; then
	echo2 "Downloading 2bit file and convert it to fasta"
	rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/${Genome}/bigZips/${Genome}.2bit ./ \
    && twoBitToFa ${Genome}.2bit ${Genome}.fa \
    && rm ${Genome}.2bit
fi

# adding TAS to Drosophila melanogaster genomes
if [[ $Genome = dm* ]] ; then
	if [[ ! -f .TAS.Done ]]; then
		declare TAS="$PIPELINE_DIRECTORY/common/supplementary/X-TAS.fa"
		if [[ ! -s $TAS ]]; then echo2 "Cannot file X-TAS.fa file, please reclone your git" error; fi
		cat ${Genome}.fa ${TAS} > ${Genome}.fa1 \
		&& mv ${Genome}.fa1 ${Genome}.fa \
		&& touch .TAS.Done
	fi
fi

if [[ ! -s ${Genome}.fa.fai ]]; then samtools faidx ${Genome}.fa; fi
if [[ ! -s ${Genome}.ChromInfo.txt ]]; then faSize -tab -detailed ${Genome}.fa > ${Genome}.ChromInfo.txt; fi

echo2 "Building genome indices"
if ! bowtieCheck BowtieIndex/genome; then 
	echo2 "Building bowtie index for genome"
	bowtie-build  ${Genome}.fa BowtieIndex/genome &> log/bowtie_genome.log \
	|| echo2 "failed to generate Bowtie index, exiting..." error; 
fi

if ! bowtie2Check Bowtie2Index/genome; then 
	echo2 "Building bowtie2 index for genome"
	bowtie2-build ${Genome}.fa Bowtie2Index/genome &> log/bowtie2_genome.log \
	|| echo2 "failed to generate Bowtie2 index, exiting..." error; 
fi

if ! bwaCheck BWAIndex/genome; then 
	echo2 "Building BWA index for genome"
	bwa index -p BWAIndex/genome ${Genome}.fa &> log/bwa_genome.log \
	|| echo2 "failed to generate BWA index, exiting..." error; 
fi

############################
# Transcriptome annotation #
############################
echo2 "Obtaining GTF file from UCSC"
if [[ ! -s gene.gtf ]]; then
	echo2 "Obtaining gene gtf file from UCSC"
	# make config for kent tools  
    if [[ ! -s ${HOME}/.hg.conf ]]; then
        cat > ${HOME}/.hg.conf << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
db.host=genome-mysql.cse.ucsc.edu
db.user=genomep
db.password=password
central.db=hgcentral
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    	chmod 600 ${HOME}/.hg.conf
    fi # end of [[ ! -s ${HOME}/.hg.conf ]]    
    genePredToGtf ${Genome} refFlat gene.gtf &> log/genePredToGtf.genome.log
fi # end of [[ -f gene.gtf ]]

if ! StarCheck STARIndex/; then
	echo2 "Building STAR index"
	STAR \
		--runMode genomeGenerate \
		--runThreadN $Threads \
		--genomeDir STARIndex \
		--genomeFastaFiles ${Genome}.fa \
		--sjdbGTFfile gene.gtf \
		--sjdbOverhang 99 \
		&> log/STAR_genome.log \
	|| echo2 "failed to generate STAR index, exiting..." error; 
fi

echo2 "Obtaining refFlat file from UCSC"
if [[ ! -s refFlat.txt ]]; then
	mysql \
		-h genome-mysql.cse.ucsc.edu \
		-u genome \
		-D ${Genome} \
		-N \
		-A \
		-e 'SELECT * FROM refFlat' \
		> refFlat.txt \
	|| echo2 "Cannot obtain refFlat file from UCSC" error
fi

# converting gtf to bed and extract the fasta
if [[ ! -s genes.bed12 ]]; then 
	echo2 "Converting GTF file to BED12"
	gtfToGenePred -ignoreGroupsWithoutExons gene.gtf genes.gp \
	&& genePredToBed genes.gp genes.bed12 \
	|| echo2 "cannot convert gene.gtf to genes.bed12" error
fi

if [[ ! -s gene.fa ]]; then
	echo2 "Extracting transcriptome sequences from genome"
	bedtools_piPipes getfasta \
		-fi ${Genome}.fa \
		-bed genes.bed12 \
		-fo gene.fa \
		-name -split -s \
	|| echo2 "cannot extract transcriptome sequences" error
fi

##################
# rRNA and miRNA #
##################
for index in rRNA hairpin mature; do
	if ! bowtieCheck BowtieIndex/${index}; then 
		echo2 "Building Bowtie Index for ${index}"; 
		bowtie-build ${index}.fa BowtieIndex/${index} &> log/bowtie_${index}.log \
		|| echo2 "failed to generate Bowtie index, exiting..." error; 
	fi
	
	if ! bowtie2Check Bowtie2Index/${index}; then 
		echo2 "Building Bowtie2 Index for ${index}";
		bowtie2-build ${index}.fa Bowtie2Index/${index} &> log/bowtie2_${index}.log \
		|| echo2 "failed to generate Bowtie2 index, exiting..." error; 
	fi
done

echo2 "align mature miRNA sequences to the hairpins to determine their positions"
if [[ ! -s mature2hairpin.uniq.bed ]]; then
	bowtie --norc -S -f -v 0 -m 1 --best --strata \
		--max mature.multiMapper.fa \
		BowtieIndex/hairpin \
		mature.fa \
		2> log/mature2hairpin.bowtie.log \
	| samtools view -uS - \
	| bedtools_piPipes bamtobed -i - \
	> mature2hairpin.uniq.bed \
	|| echo2 "cannot align mature.fa to hairpin.fa by bowtie and generate bed file" error
fi 
if [[ ! -s mature2hairpin.multi.bed ]]; then 
	bowtie --norc -S -f -v 0 -a   --best --strata \
		BowtieIndex/hairpin \
		mature.multiMapper.fa \
		2> log/multimapper_mature2hairpin.bowtie.log \
	| samtools view -uS - \
	| bedtools_piPipes bamtobed -i - > mature2hairpin.multi.bed \
	|| echo2 "cannot align mature.multimapper.fa to hairpin.fa by bowtie and generate bed file" error
fi 
if [[ ! -s mature2hairpin.allMapper.bed ]]; then 
	cat \
		mature2hairpin.uniq.bed \
		mature2hairpin.multi.bed \
	> mature2hairpin.allMapper.bed
fi

##############
# Transposon #
##############
echo2 "Building Bowtie/BWA index for transposon annotation"
if ! bowtieCheck BowtieIndex/transposon; then bowtie-build transposon.fa BowtieIndex/transposon &> log/bowtie_transposon.log ; fi
if ! bowtie2Check Bowtie2Index/transposon; then bowtie2-build transposon.fa Bowtie2Index/transposon &> log/bowtie2_transposon.log; fi
if [[ ! -f transposon.sizes ]]; then faSize -tab -detailed transposon.fa > transposon.sizes; fi

if [[ ! -s transposon.eref ]]; then
	mkdir -p transposon \
	&& faSplit byname transposon.fa transposon/ \
	&& for i in transposon/*; do echo -e "`basename ${i%.fa}`\t`readlink -f $i`"; done \
	> ${Genome}.transposon.eref
fi

#################
# piRNA cluster #
#################
echo2 "Building Bowtie/BWA index for piRNA cluster"
if [[ ! -s piRNAcluster.fa ]]; then 
	bedtools_piPipes getfasta \
		-fi  ${Genome}.fa \
		-bed piRNAcluster.bed \
		-fo piRNAcluster.fa \
		-name -split -s \
	|| echo2 "cannot extract piRNA cluster sequences" error
fi

if ! bowtieCheck BowtieIndex/piRNAcluster; then bowtie-build piRNAcluster.fa BowtieIndex/piRNAcluster  &> log/bowtie_piRNAcluster.log || echo2 "failed to generate bowtie index for piRNAcluster.fa" error; fi 
if [[ ! -f piRNAcluster.sizes ]]; then faSize -tab -detailed piRNAcluster.fa > piRNAcluster.sizes || echo2 "failed to generate piRNAcluster.sizes file" error; fi

#################
# Transcriptome #
#################
case $Genome in
dm3|dm6)
	echo2 "Building Bowtie2 index for transposon + genes"
	if [[ ! -s transcriptome.fa ]]; then cat gene.fa transposon.fa > transcriptome.fa; fi
;;
*)
	echo2 "Building Bowtie2 index for transposon + piRNA cluster + genes"
	if [[ ! -s transcriptome.fa ]]; then cat gene.fa piRNAcluster.fa transposon.fa > transcriptome.fa; fi
;;
esac

if ! bowtie2Check Bowtie2Index/transcriptome; then bowtie2-build transcriptome.fa Bowtie2Index/transcriptome &> log/bowtie2_transcriptome.log; fi
if [[ ! -s transcriptome.sizes ]]; then faSize -tab -detailed transcriptome.fa > transcriptome.sizes; fi

########
# Done #
########
# TODO: lock file
sort -u <(cat $PIPELINE_DIRECTORY/common/genome_supported.txt) <(echo $Genome) > $PIPELINE_DIRECTORY/common/genome_supported.txt.${Genome} \
&& mv $PIPELINE_DIRECTORY/common/genome_supported.txt.${Genome} $PIPELINE_DIRECTORY/common/genome_supported.txt  
