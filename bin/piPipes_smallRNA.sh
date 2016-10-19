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
declare -x  SmallRnaModuleVersion=${PROG_VERSION}
declare -xi Threads=${DEFAULT_THREADS}
declare -x  Genome
declare -x  AnnotationDir
declare -x  InputFastq
declare -x  Outdir
declare -x  Normmethod="unique"
declare -a  FilterMappingFileList=()
declare -a  PreGenomeMappingFileList=()
declare -a  PostGenomeMappingFileList=()

declare -a  RequiredPrgrams=( \
    sort md5sum gawk grep python gs Rscript \
    samtools bowtie ParaFly bedGraphToBigWig \
    bedtools_piPipes piPipes_bed2Summary piPipes_fastq_to_insert piPipes_insertBed_to_bed2 \
    )

#########
# const #
#########
declare -ix NucCompExtLen=30

#########
# USAGE #
#########
usage () {
cat << EOF

small RNA Seq pipeline single library mode v$SmallRnaModuleVersion from the $BOLD$PACKAGE_NAME$RESET
$SMALLRNA_INTRO${RESET}
Please email $CONTACT_EMAILS for any questions or bugs.
Thank you for using it.

${UNDERLINE}usage${RESET}:
    piPipes small \ 
        -i input.fq[.gz] \ 
        -g dm3 \ 
        -N uniqueXmiRNA [unique] \ 
        -o output_directory [current working directory] \ 
        -c cpu [8] \  
        -P miniwhite.fa \ 
        -O gfp.fa,luciferase.fa 

OPTIONS:
    -h      Show this message

    -v      Print out the version

${REQUIRED}[ required ]
    -i      Input file in fastq or gzipped fastq format; Needs adaptor and barcode removed
            Since this small RNA pipeline does not consider quality, we strongly recommend a quality filtering step.
    
    -g      Genome assembly name, like mm9 or dm3
            Check "$PIPELINE_DIRECTORY/common/genome_supported.txt" for genome assemblies currently installed.
            Use "install" to install new genome.

${OPTIONAL}[ optional ]
    -N      Normalization method, choose from " input | rRNA | unique | uniqueXmiRNA | all | allXmiRNA | miRNA "
                ${UNDERLINE}unique${RESET}${OPTIONAL}: use non-rRNA genomic unique mappers <default>.
                ${UNDERLINE}input${RESET}${OPTIONAL}: use the number of reads input to the pipeline, this will include genome 
            unmappers. But might be useful when you have additional sequence in the genome, like a transgene.
                ${UNDERLINE}rRNA${RESET}${OPTIONAL}: use the number of reads mapped to rRNA.
                ${UNDERLINE}uniqueXmiRNA${RESET}${OPTIONAL}: use non-rRNA genomic unique mappers excluding microRNAs <for oxidized library from piRNA mutant>.
                ${UNDERLINE}all${RESET}${OPTIONAL}: use non-rRNA genomic all mappers including microRNAs.
                ${UNDERLINE}allXmiRNA${RESET}${OPTIONAL}: use non-rRNA genomic all mappers excluding microRNAs.
                ${UNDERLINE}miRNA${RESET}${OPTIONAL}: use microRNAs. normalized to: reads per millions of miRNA <for unoxidized library from piRNA mutant>.
            * Different normalization methods, including "siRNA", are available in the dual sample mode.
            * You are able to run the same library multiple times with different normalization method. They will not collapse.

    -c      Number of Threadss to use, default: $Threads
    
    -o      Output directory, default: current directory $PWD
    
    -F      A Fasta file used for filtering (other than rRNA precursor sequence provide).
    
    -P      A Fasta files used to do pre-genomic mapping and analysis. For example, given "-P miniwhite.fa -P virus.fa", 
            after removing reads mappable to rRNA and miRNA hairpin, reads are mapped to miniWhite sequence first. 
            Only the non-miniWhite-mappers are mapped to virus sequence. And only the non-miniWhite, non-virus mappers will be 
            used in the genome mapping and further analysis.
    
    -O      A Fasta files used to do post-genomic mapping and analysis. For example, given "-O gfp.fa,luciferase.fa", 
            after removing reads mappable rRNA, miRNA hairpin and genome, reads are mapped to gfp sequence first. 
            Only the non-genome non-gfp mappers are mapped to luciferase sequence. If more than one sequences are put in 
            one Fasta file, they will be treated equally. 
    
    For -F, -P and -O options, you can provide more than one fasta files, each with their own tag, such as "-P white.fa -P virus.fa". 
    They will be mapped and analyzed sequentially.
    

${RESET}

EOF
echo -e "${COLOR_END}"
}

########
# ARGS #
########
while getopts "hi:c:o:g:vN:F:P:O:" OPTION; do
    case $OPTION in
        h)  usage && exit 0;;
        
        v)  echo2 "SmallRnaModuleVersion: v$SmallRnaModuleVersion" && exit 0;;
        
        i)  InputFastq=$(readlink -f "${OPTARG}"); 
            assertFileExists "${InputFastq}";;
        
        o)  Outdir=$(readlink -f "${OPTARG}");;
        
        c)  Threads=$OPTARG;;
        
        g)  Genome=${OPTARG}
            check_genome $Genome
            AnnotationDir=${PIPELINE_DIRECTORY}/common/${Genome}/
            declare -x COMMON_FOLDER=${AnnotationDir} # backward compatibility
            ;;
        
        N)  Normmethod=$(echo "${OPTARG}" | tr '[A-Z]' '[a-z]');;
        
        F)  declare optarg=$(readlink -f "${OPTARG}"); 
            assertFileExists $optarg;
            FilterMappingFileList+=( "${optarg}" );;
        
        P)  declare optarg=$(readlink -f "${OPTARG}"); 
            assertFileExists $optarg;
            PreGenomeMappingFileList+=( "${optarg}" );;
        
        O)  declare optarg=$(readlink -f "${OPTARG}");
            assertFileExists $optarg;
            PostGenomeMappingFileList+=( "${optarg}" );;
        
        *)  usage && exit 1;;
    esac
done

# if InputFastq or GENOME is undefined, print out usage and exit
if [[ -z ${InputFastq} ]]; then echo2 "Missing option -i for input fastq file" error; fi
if [[ -z ${Genome} ]]; then echo2 "Missing option -g for specifying which genome assembly to use" error; fi
declare FqName=$(basename "${InputFastq}")
declare -x Prefix=${FqName%.f[qa]*}
if [[ -z $Outdir ]]; then Outdir=$PWD; fi
if [[ ! -d $Outdir ]] ; then mkdir -p $Outdir; fi
assertDirExists $Outdir
assertDirWritable $Outdir
cd ${Outdir} || echo2 "Cannot access directory ${Outdir}... Exiting..." error

########################
# running binary check #
########################
for program in "${RequiredPrgrams[@]}"; do assertBinExists $program; done

#################################
# creating output files/folders #
#################################
declare -x TableDir=tables && mkdir -p $TableDir
declare -x StatsTable=$TableDir/${Prefix}.basic_stats
declare -x PdfDir=pdfs && mkdir -p $PdfDir
declare -x ReadsDir=input_read_files && mkdir -p $ReadsDir
declare -x JobDir=jobs && mkdir -p $JobDir
declare -x rRnaDir=$JobDir/rRNA_mapping && mkdir -p $rRnaDir
declare -x miRnaDir=$JobDir/hairpins_mapping && mkdir -p $miRnaDir
declare -x FilterDir=$JobDir/custom_filter && mkdir -p $FilterDir
declare -x PreGenomeDir=$JobDir/pre_genome_mapping && mkdir -p $PreGenomeDir
declare -x PostGenomeDir=$JobDir/post_genome_mapping && mkdir -p $PostGenomeDir
declare -x GenomeMappingDir=$JobDir/genome_mapping && mkdir -p $GenomeMappingDir
declare -x IntersectDir=$JobDir/intersect_genomic_features && mkdir -p $IntersectDir
declare -x SummaryDir=summaries && mkdir -p $SummaryDir
declare -x BigwigDir=$JobDir/bigWig_normalized_by_$Normmethod && mkdir -p $BigwigDir
declare -x TransposonDir=$JobDir/transposon_piRNAcluster_mapping_normalized_by_$Normmethod && mkdir -p $TransposonDir
declare -x LogDir=log && mkdir -p ${LogDir}
declare -x SentinelDir=sentinels && mkdir -p ${SentinelDir}
declare -x NormScale=1
declare -xi Step=1
declare -x RunUid=$(echo "${InputFastq}" | md5sum | cut -d" " -f1)
declare -x CommonFolder=$PIPELINE_DIRECTORY/common/$Genome
. $CommonFolder/variables
for var in rRNA_MM hairpin_MM Genome_MM transposon_MM siRNA_bot siRNA_top piRNA_bot piRNA_top; do
    if [[ -z ${!var} ]]; then echo2 "value of $var is missing, please check your $CommonFolder/variables file" error; fi
done
declare -x GenomeFa=$CommonFolder/${Genome}.fa && assertFileExists ${GenomeFa}
declare -x ChromSize=$CommonFolder/${Genome}.ChromInfo.txt && assertFileExists ${ChromSize}
declare -x BOWTIE_INDEXES=$CommonFolder/BowtieIndex # for Bowtie

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] small RNA pipeline single library mode version $SmallRnaModuleVersion"

########################################
## Pre Processing before any Mapping ###
########################################
declare Insert=$ReadsDir/${Prefix}.insert
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.Done ]]; then
    echo2 "Converting fastq format into insert format by dumping quality strings"
    piPipes_fastq_to_insert "${InputFastq}" ${Insert} \
    && touch $SentinelDir/${RunUid}.status.${Step}.Done \
    || echo2 "fq2insert failed" error
else 
    echo2 "Conversion from fastq to insert has been done previously" warning
fi
let Step+=1

####################
## pre filtering ###
####################
declare -xi CurrentMM=0

declare -x CurrentInput=$Insert
declare -x CurrentInputPrefix=${CurrentInput%.insert}
declare -x CurrentInputName=$(basename $CurrentInput) 
declare -x CurrentInputNamePrefix=${CurrentInputName%.insert}
declare -x NextInput
declare -x CurrentTarget
declare -x CurrentTargetName
declare -x CurrentTargetNamePrefix
declare -x CurrentOutdir
declare -x BowtieIndexName
declare -x LogFile

# parsing customer defined pre-genomic mapping variables
if [[ ! -z ${FilterMappingFileList[@]+"${FilterMappingFileList[@]}"} ]]; then
    echo2 "Mapping to customer defined filtering indexes"
    # eval $(echo $FilterMappingFileList | awk 'BEGIN{FS=","}{printf "export FilterMappingFiles=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}')
    for target in "${FilterMappingFileList[@]}"; do
        CurrentTarget=$target # already full path from argparse
        CurrentTargetName=$(basename $CurrentTarget) 
        CurrentTargetNamePrefix=${CurrentTargetName%.fa}
        NextInput=${CurrentInputPrefix}.x_${CurrentTargetNamePrefix}.insert
        CurrentOutdir=$FilterDir/${CurrentTargetNamePrefix}
        BowtieIndexName=$CurrentOutdir/$CurrentTargetNamePrefix
        LogFile=${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log
        if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_filtering_mapping ]]; then
            echo2 "Mapping to ${CurrentTargetNamePrefix}"
            mkdir -p $CurrentOutdir \
            || echo2 "Cannot create directory $CurrentOutdir, please check the permission. And try to only use letter and number to name the Fasta file" error
            
            bowtie-build \
                $CurrentTarget \
                $BowtieIndexName \
                &> ${LogDir}/bowtie_build.${CurrentTargetNamePrefix}.log \
            || echo2 "Failed to build the bowtie index for $targetFa" error
            
            bowtie -r -v 0 -a --best --strata -p $Threads -S \
                --un ${NextInput} \
                $BowtieIndexName \
                $CurrentInput \
                1> /dev/null \
                2> ${LogFile}  \
            && rm -rf ${BowtieIndexName}*ebwt \
            && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_filtering_mapping
        else 
            echo2 "Mapping to ${CurrentTargetNamePrefix} has been done previously" warning
        fi
        CurrentInput=${NextInput}
        CurrentInputPrefix=${CurrentInput%.insert}
        CurrentInputName=$(basename $CurrentInput) 
        CurrentInputNamePrefix=${CurrentInputName%.insert}
    done
else 
    echo2 "User has not specified any filtering sequenceing" warning
fi

################
# rRNA removal #
################
# getting rid of sequences mapping to rRNA, we use -k 1 option for speed purpose
NextInput=${CurrentInputPrefix}.x_rRNA.insert
CurrentTargetName=rRNA
CurrentTargetNamePrefix=rRNA
BowtieIndexName=rRNA
CurrentMM=$rRNA_MM
LogFile=${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix} ]]; then
    echo2 "Mapping to ${CurrentTargetNamePrefix}, with ${CurrentMM} mismatch(es) allowed"
    declare -i totalReads=$(awk '{a+=$2}END{printf "%d", a}' ${CurrentInput}) \
    && echo $totalReads > ${TableDir}/${RunUid}.totalReads \
    && bowtie -r -S -v $CurrentMM -k 1 -p $Threads \
        --un $NextInput \
        ${BowtieIndexName} \
        ${CurrentInput} \
        1> /dev/null \
        2> ${LogDir}/${Prefix}.bowtie2rRNA.log \
    && declare -i nonrRNAReads=$(awk '{a+=$2}END{printf "%d", a}' ${NextInput}) \
    && echo $nonrRNAReads > $TableDir/${RunUid}.nonrRNAReads \
    && declare -i rRNAReads=$((totalReads-nonrRNAReads)) \
    && echo $rRNAReads > $TableDir/${RunUid}.rRNAReads \
    && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix} \
    || echo2 "mapping to rRNA failed" error
else 
    echo2 "Mapping to ${CurrentTargetNamePrefix} has been done previously" warning
fi
let Step+=1
CurrentInput=${NextInput}
CurrentInputPrefix=${CurrentInput%.insert}
CurrentInputName=$(basename $CurrentInput) 
CurrentInputNamePrefix=${CurrentInputName%.*}
# reading values from file, this is for resuming the job, which won't run the previous step
totalReads=$(cat $TableDir/${RunUid}.totalReads)
rRNAReads=$(cat $TableDir/${RunUid}.rRNAReads)
nonrRNAReads=$(cat $TableDir/${RunUid}.nonrRNAReads)

#########################
# miRNA hairpin Mapping #
#########################
declare HairpinInsert=${CurrentInputPrefix}.hairpin.insert 
declare NextInput=${CurrentInputPrefix}.x_hairpin.insert 
declare HairpinBed2=$miRnaDir/${CurrentInputNamePrefix}.hairpin_v${hairpin_MM}m1.bed2
declare HairpinLendis=$miRnaDir/${CurrentInputNamePrefix}.hairpin_v${hairpin_MM}m1.lendis 
declare HairpinGenomeBed2=$GenomeMappingDir/${CurrentInputNamePrefix}.hairpin.${Genome}v${Genome_MM}a.bed2
CurrentTargetNamePrefix=hairpin
CurrentMM=$hairpin_MM
CurrentOutdir=$miRnaDir
BowtieIndexName=hairpin
LogFile=$LogDir/${CurrentInputNamePrefix}.hairpin.${Genome}v${Genome_MM}a.log # log file for hairpin mapping
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.hairpin_mapping ]]; then
    echo2 "Mapping to microRNA Hairpin, with $hairpin_MM mismatch(es) allowed; only keep unique mappers"
    bowtie -r --norc -v $CurrentMM -m 1 \
        --best --strata -p $Threads -S \
        --al $HairpinInsert \
        --un $NextInput \
        $BowtieIndexName \
        $CurrentInput \
        2> $LogDir/${Prefix}.hairpin_unique.log \
    | samtools view -bSF 0x4 - 2>/dev/null \
    | bedtools_piPipes bamtobed -i - \
        > ${CurrentInputPrefix}.${BowtieIndexName}m1.bed \
    && piPipes_insertBed_to_bed2 \
        $CurrentInput \
        ${CurrentInputPrefix}.${BowtieIndexName}m1.bed \
        > $HairpinBed2 \
    && rm -r ${CurrentInputPrefix}.${BowtieIndexName}m1.bed  \
    && bed2lendis $HairpinBed2 > $HairpinLendis \
    && bowtie -r -v $Genome_MM -a --best --strata -p $Threads \
        -S \
        genome \
        $HairpinInsert \
        2> $LogFile \
    | samtools view -uS -F0x4 - 2>/dev/null \
    | bedtools_piPipes bamtobed -i - \
        > ${miRnaDir}/${CurrentInputNamePrefix}.hairpin.${Genome}v${Genome_MM}a.bed \
    && piPipes_insertBed_to_bed2 \
        $HairpinInsert \
        ${miRnaDir}/${CurrentInputNamePrefix}.hairpin.${Genome}v${Genome_MM}a.bed \
        > $HairpinGenomeBed2 \
    && rm -f ${miRnaDir}/${CurrentInputNamePrefix}.hairpin.${Genome}v${Genome_MM}a.bed \
    && declare -i hairpinReads=$(bedwc $HairpinBed2) \
    && echo $hairpinReads > $TableDir/${RunUid}.hairpinReads \
    && touch $SentinelDir/${RunUid}.status.${Step}.hairpin_mapping
else 
    echo2 "Mapping to miRNA hairpin has been down previously" warning
fi
let Step+=1 
declare -i hairpinReads=$(cat $TableDir/${RunUid}.hairpinReads)
CurrentInput=${NextInput}
CurrentInputPrefix=${CurrentInput%.insert}
CurrentInputName=$(basename $CurrentInput) 
CurrentInputNamePrefix=${CurrentInputName%.*}

# run miRNA heterogeneity analysis
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.miRNA_heterogeneity ]]; then
    echo2 "Calculate microRNA heterogeneity"
    piPipes_calculate_miRNA_heterogeneity \
        $CommonFolder/mature2hairpin.uniq.bed \
        ${HairpinBed2} \
        1> ${HairpinBed2%.bed*}.sum \
    && touch $SentinelDir/${RunUid}.status.${Step}.miRNA_heterogeneity
else 
    echo2 "Calculate microRNA heterogeneity has been done previously" warning
fi
let Step+=1 
declare PostFilterInput=$CurrentInput # for transposon direct mapping

#############################
# custom pre-genome mapping #
#############################
CurrentMM=0 # haven't implement method to take mismatch # from user
# parsing customer defined pre-genomic mapping variables
if [[ ! -z ${PreGenomeMappingFileList[@]+"${PreGenomeMappingFileList[@]}"} ]]; then
    echo2 "Pre-genome mapping to customer defined sequences"
    for CurrentTarget in "${PreGenomeMappingFileList[@]}"; do
        CurrentTargetName=$(basename $CurrentTarget) 
        CurrentTargetNamePrefix=${CurrentTargetName%.fa}
        NextInput=${CurrentInputPrefix}.x_${CurrentTargetNamePrefix}.insert
        CurrentOutdir=$PreGenomeDir/${CurrentTargetNamePrefix}
        mkdir -p $CurrentOutdir || echo2 "Cannot create directory $CurrentOutdir, please check the permission. And try to only use letter and number to name the Fasta file" error
        BowtieIndexName=$CurrentOutdir/$CurrentTargetNamePrefix
        LogFile=${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log

        if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_pregenome_mapping ]]; then
            echo2 "Mapping to ${CurrentTargetName}"
            bash piPipes_map_smallRNA_to_target.sh \
            && PDFs=$PdfDir/${CurrentTargetNamePrefix}/*pdf \
            && gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
                -sOutputFile=$PdfDir/PreGenome.${CurrentTargetNamePrefix}.pdf \
                ${PDFs} \
            && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_pregenome_mapping \
            || echo2 "failed to run pre-genome mapping for ${CurrentTargetNamePrefix}" error
        else 
            echo2 "Mapping to ${CurrentTargetName} has been done previously" warning
        fi
        CurrentInput=${NextInput}
        CurrentInputPrefix=${CurrentInput%.insert}
        CurrentInputName=$(basename $CurrentInput) 
        CurrentInputNamePrefix=${CurrentInputName%.*}
    done
else 
    echo2 "User has not specified any pre-genomic mapping" warning
fi
let Step+=1 
            
##################
# GENOME Mapping #
##################
CurrentMM=$Genome_MM
CurrentTargetNamePrefix=${Genome}
BowtieIndexName=genome
declare -x GenomeAllmapBed2=$GenomeMappingDir/${CurrentInputNamePrefix}.genome_v${CurrentMM}a.bed2
declare -x GenomeUniquemapBed2=$GenomeMappingDir/${CurrentInputNamePrefix}.genome_v${CurrentMM}a.unique.bed2
GenomeAllmapLogFile=$LogDir/${CurrentInputNamePrefix}.genome_v${CurrentMM}a.log 
NextInput=${CurrentInputPrefix}.x_${CurrentTargetNamePrefix}.insert

declare -i totalGenomicMapCount
declare -i uniqueGenomicMapCount
declare -i multipGenomicMapCount

if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_mapping ]]; then
    echo2 "Mapping to ${CurrentTargetNamePrefix}, with ${CurrentMM} mismatch(es) allowed"
    bowtie -r -v $CurrentMM -a --best --strata -p $Threads \
        --al ${CurrentInputPrefix}.${CurrentTargetNamePrefix}.al.insert \
        --un $NextInput \
        -S \
        $BowtieIndexName \
        ${CurrentInput} \
        2> $GenomeAllmapLogFile \
    | samtools view -uS -F0x4 - 2>/dev/null \
    | bedtools_piPipes bamtobed -i - \
        > ${GenomeAllmapBed2}.temp \
    && piPipes_insertBed_to_bed2 \
        $CurrentInput \
        ${GenomeAllmapBed2}.temp \
        > ${GenomeAllmapBed2} \
    && rm -rf ${GenomeAllmapBed2}.temp \
    && awk 'BEGIN{OFS="\t"}{if ($5==1) print $0}' \
        ${GenomeAllmapBed2} \
        > ${GenomeUniquemapBed2} \
    && totalGenomicMapCount=$(bedwc ${GenomeAllmapBed2}) \
    && echo $totalGenomicMapCount > $TableDir/${RunUid}.totalGenomicMapCount \
    && uniqueGenomicMapCount=$(bedwc ${GenomeUniquemapBed2}) \
    && echo $uniqueGenomicMapCount > $TableDir/${RunUid}.uniqueGenomicMapCount \
    && multipGenomicMapCount=$((totalGenomicMapCount-uniqueGenomicMapCount)) \
    && echo $multipGenomicMapCount > $TableDir/${RunUid}.multipGenomicMapCount \
    && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_mapping
else
    echo2 "Mapping to ${CurrentTargetNamePrefix} has already been down" warning
fi
let Step+=1

totalGenomicMapCount=$(cat $TableDir/${RunUid}.totalGenomicMapCount)
uniqueGenomicMapCount=$(cat $TableDir/${RunUid}.uniqueGenomicMapCount)
multipGenomicMapCount=$(cat $TableDir/${RunUid}.multipGenomicMapCount)

CurrentInput=${NextInput}
CurrentInputPrefix=${CurrentInput%.insert}
CurrentInputName=$(basename $CurrentInput) 
CurrentInputNamePrefix=${CurrentInputName%.*}

##############################
# custom post-genome mapping #
##############################
CurrentMM=0 # haven't implement method to take mismatch # from user
if [[ ! -z ${PostGenomeMappingFileList[@]+"${PostGenomeMappingFileList[@]}"} ]]; then
    echo2 "Post-genome mapping to customer defined sequences"
    for CurrentTarget in "${PostGenomeMappingFileList[@]}"; do
        CurrentTargetName=$(basename $CurrentTarget) 
        CurrentTargetNamePrefix=${CurrentTargetName%.fa}
        NextInput=${CurrentInputPrefix}.x_${CurrentTargetNamePrefix}.insert
        CurrentOutdir=$PostGenomeDir/${CurrentTargetNamePrefix}
        mkdir -p $CurrentOutdir || echo2 "Cannot create directory $CurrentOutdir, please check the permission. And try to only use letter and number to name the Fasta file" error
        BowtieIndexName=$CurrentOutdir/$CurrentTargetNamePrefix
        LogFile=${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log

        if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_postgenome_mapping ]]; then
            echo2 "Mapping to ${CurrentTargetName}"
            bash piPipes_map_smallRNA_to_target.sh \
            && PDFs=$PdfDir/${CurrentTargetNamePrefix}/*pdf \
            && gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
                -sOutputFile=$PdfDir/PostGenome.${CurrentTargetNamePrefix}.pdf \
                ${PDFs} \
            && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetNamePrefix}_postgenome_mapping \
            || echo2 "failed to run post-genome mapping for ${CurrentTargetNamePrefix}" error
        else 
            echo2 "Mapping to ${CurrentTargetName} has been done previously" warning
        fi
        CurrentInput=${NextInput}
        CurrentInputPrefix=${CurrentInput%.insert}
        CurrentInputName=$(basename $CurrentInput) 
        CurrentInputNamePrefix=${CurrentInputName%.*}
    done
fi
let Step+=1 

#####################
# Length Separation #
#####################
declare GenomeAllmapSiRNABed=${GenomeAllmapBed2%bed2}siRNA.bed2
declare GenomeAllmapPiRNABed=${GenomeAllmapBed2%bed2}piRNA.bed2
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.sep_length ]]; then
    para_file=${RANDOM}${RANDOM}.para && \
    echo "awk '\$3-\$2>=$siRNA_bot && \$3-\$2<=$siRNA_top' ${GenomeAllmapBed2} > ${GenomeAllmapSiRNABed}" > $para_file \
    && echo "awk '\$3-\$2>=$piRNA_bot && \$3-\$2<=$piRNA_top' ${GenomeAllmapBed2} > ${GenomeAllmapPiRNABed}" >> $para_file \
    && ParaFly -c $para_file -CPU $Threads -failed_cmds ${para_file}.failedCommands &>/dev/null \
    || echo2 "Failed to separate siRNA and piRNA from genomic bed file" error
    rm -rf ${para_file}* \
    && touch $SentinelDir/${RunUid}.status.${Step}.sep_length
fi
let Step+=1 

# plotting length distribution
echo2 "Plotting length distribution"
declare GenomeAllmapBase=$(basename $GenomeAllmapBed2)
declare GenomeUniquemapBase=$(basename $GenomeUniquemapBed2)
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.plotting_length_dis ]]; then
    awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' \
        ${GenomeAllmapBed2} \
    | sort -k1,1n > $TableDir/${GenomeAllmapBase}.lendis \
    && awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' \
        ${GenomeUniquemapBed2} \
    | sort -k1,1n > $TableDir/${GenomeUniquemapBase}.lendis \
    && Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R $TableDir/${GenomeAllmapBase}.lendis    $PdfDir/${Step}.$(basename ${GenomeAllmapBed2%bed2})-hairpin &> /dev/null \
    && Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R $TableDir/${GenomeUniquemapBase}.lendis $PdfDir/${Step}.$(basename ${GenomeUniquemapBed2%bed2})-hairpin &> /dev/null \
    && awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' $TableDir/${GenomeAllmapBase}.lendis    $HairpinLendis | sort -k1,1n > $TableDir/${GenomeAllmapBase}.+hairpin.lendis \
    && awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' $TableDir/${GenomeUniquemapBase}.lendis $HairpinLendis | sort -k1,1n > $TableDir/${GenomeUniquemapBase}.+hairpin.lendis \
    && Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R $TableDir/${GenomeAllmapBase}.+hairpin.lendis    $PdfDir/${Step}.$(basename ${GenomeAllmapBed2%bed2})+hairpin &> /dev/null \
    && Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R $TableDir/${GenomeUniquemapBase}.+hairpin.lendis $PdfDir/${Step}.$(basename ${GenomeUniquemapBed2%bed2})+hairpin &> /dev/null \
    && touch $SentinelDir/${RunUid}.status.${Step}.plotting_length_dis
fi
let Step+=1

##################
# Print to table #
##################
# change dual library mode normalization method if change here
   echo -e "$(basename ${InputFastq})\ttotal_input_reads\t${totalReads}" > $StatsTable \
&& echo -e "$(basename ${InputFastq})\trRNA_reads\t${rRNAReads}" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tmiRNA_hairpin_reads\t${hairpinReads}" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tgenome_mapping\t$((totalGenomicMapCount+hairpinReads))" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tgenome_mapping_without_hairpin\t${totalGenomicMapCount}" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tunique_genome_mapping\t$((uniqueGenomicMapCount+hairpinReads))" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tunique_genome_mapping_without_hairpin\t${uniqueGenomicMapCount}" >> $StatsTable \
&& echo -e "$(basename ${InputFastq})\tmulti_genome_mapping_without_hairpin\t${multipGenomicMapCount}" >> $StatsTable
declare -x GenomicHairpinReads=$((totalGenomicMapCount+hairpinReads))

# normalization method
# input | rRNA | unique | uniqueXmiRNA | all | allXmiRNA | miRNA
case "$Normmethod" in
input)
    NormScale=$(head -1 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
rrna)
    NormScale=$(head -2 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
mirna)
    NormScale=$(head -3 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
all)
    NormScale=$(head -4 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
allxmirna)
    NormScale=$(head -5 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
unique)
    NormScale=$(head -6 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
uniquexmirna)
    NormScale=$(head -7 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
*)
    echo2 "unrecognized normalization option: $Normmethod; using the default method" warning
    NormScale=$(head -6 $StatsTable | tail -1 | awk '{print 1000000/$NF}')
;;
esac
echo $NormScale > $TableDir/${RunUid}.depth

####################################
# Intersecting with GENOME Feature #
####################################
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.intersect_with_genomic_features ]]; then
    echo2 "Intersecting with genomic features, make length distribution, nucleotide fraction for siRNA/piRNA assigned to each feature"
    bash piPipes_intersect_smallRNA_with_genomic_features.sh \
        $TableDir/${GenomeAllmapBase%.bed2} \
        $IntersectDir \
        1>&2 \
        && touch $SentinelDir/${RunUid}.status.${Step}.intersect_with_genomic_features
else 
    echo2 "Intersecting with genomic features has been done previously" warning
fi
let Step+=1

#######################
# Making BigWig Files #
#######################
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.bigWig_$Normmethod ]]; then
    echo2 "Making bigWig files for genome browser"
    bash piPipes_smallRNA_bed2_to_bw.sh \
        ${GenomeAllmapBed2} \
        ${ChromSize} \
        ${NormScale} \
        $BigwigDir \
    && bash piPipes_smallRNA_bed2_to_bw.sh \
        ${GenomeAllmapBed2%bed2}piRNA.bed2 \
        ${ChromSize} \
        ${NormScale} \
        $BigwigDir \
    && touch $SentinelDir/${RunUid}.status.${Step}.bigWig_$Normmethod
else 
    echo2 "bigWig files have already been made previously" warning
fi
let Step+=1

##############################################
# Direct mapping to transposon/piRNA cluster #
##############################################
CurrentInput=$PostFilterInput # use post-filtering/post-miRNA reads
CurrentInputPrefix=${CurrentInput%.insert}
CurrentInputName=$(basename $CurrentInput) 
CurrentInputNamePrefix=${CurrentInputName%.*}
CurrentMM=transposon_MM
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.direct_transposon_mapping ]]; then
    echo2 "Direct mapping to transposon and piRNA cluster and make distribution plot"
    . $CommonFolder/genomic_features
    if [[ ! -z ${PreGenomeMappingFileList[@]+"${PreGenomeMappingFileList[@]}"} ]]; then 
        for BowtieIndexName in "${DIRECT_MAPPING[@]}"; do \
            echo2 "Direct mapping to $BowtieIndexName"
            CurrentTargetNamePrefix=${BowtieIndexName}
            CurrentOutdir=$TransposonDir/${CurrentTargetNamePrefix} \
            && mkdir -p $CurrentOutdir || echo2 "Cannot create directory $CurrentOutdir, please check the permission. And try to only use letter and number to name the Fasta file" error
            CurrentTargetSize=$CommonFolder/${CurrentTargetNamePrefix}.sizes 
            if [[ ! -s $CurrentTargetSize ]]; then echo2 "Cannot find size file for ${BowtieIndexName}" error; fi
            LogFile=${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log

            bowtie -r -a --best --strata \
                -v ${CurrentMM} \
                -p $Threads \
                -S \
                ${BowtieIndexName} \
                ${CurrentInput} \
                2> ${LogFile} \
            | samtools view -uS -F0x4 - 2>/dev/null \
            | samtools sort -o -@ $Threads - foo 2>/dev/null \
            | bedtools_piPipes bamtobed -i - \
                > $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.bed \
            && piPipes_insertBed_to_bed2 \
                $CurrentInput \
                $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.bed \
                > $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.bed2 \
            && rm $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.bed \
            && piPipes_bed2Summary \
                -5 \
                -i $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.bed2 \
                -c $CurrentTargetSize \
                -o $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.summary \
            && Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R \
                $CurrentOutdir/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.summary \
                $CurrentOutdir/${CurrentInputNamePrefix} \
                $Threads \
                $NormScale \
                &> /dev/null \
            && PDFs=$CurrentOutdir/${CurrentInputNamePrefix}*pdf \
            && gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/${Step}.${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}.pdf ${PDFs} \
            || echo2 "Mapping to ${BowtieIndexName} failed" error
        done
    else 
        echo2 "DIRECT_MAPPING is empty, please edit your $CommonFolder/genomic_features file" warning
    fi
    touch $SentinelDir/${RunUid}.status.${Step}.direct_transposon_mapping
else 
    echo2 "Direct mapping has been done previously" warning
fi 
let Step+=1

################
# Joining Pdfs #
################
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.merge_pdfs ]]; then
    echo2 "Merging pdfs"
    Pdfs=$PdfDir/*pdf \
    && gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
        -sOutputFile=$PdfDir/${Prefix}.${Genome}.small_RNA_pipeline.v${SmallRnaModuleVersion}.pdf \
        ${Pdfs} \
    && touch $SentinelDir/${RunUid}.status.${Step}.merge_pdfs
else 
    echo2 "Pdf merging has been done before" warning
fi

#############
# finishing #
#############
echo2 "Finished running ${PACKAGE_NAME} small RNA pipeline version $SmallRnaModuleVersion"
echo2 "---------------------------------------------------------------------------------"
touch ${SentinelDir}/${FqName}.${Genome}.${SmallRnaModuleVersion}
