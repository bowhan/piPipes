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
declare -x  SmallRnaModuleVersion=2.0.0
declare -xi Threads=8
declare -x  Genome=""
declare -x  InputFastq=""
declare -x  Outdir=""
declare -x  Normmethod="unique"
declare     FilterMappingFileList=""
declare     PreGenomeMappingFileList=""
declare     PostGenomeMappingFileList=""

declare -a  RequiredPrgrams=( \
    sort md5sum awk grep python gs Rscript \
    samtools bowtie ParaFly bedGraphToBigWig \
    bedtools_piPipes piPipes_bed2Summary piPipes_fastq_to_insert piPipes_insertBed_to_bed2 \
    )

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
            Check "$PIPELINE_DIRECTORY/common/genome_supported.txt" for genome assemblies currently installed;
            Use "install" to install new genome
${OPTIONAL}[ optional ]
    -N      Normalization method, choose from " input | rRNA | unique | uniqueXmiRNA | all | allXmiRNA | miRNA "
            unique: use non-rRNA genomic unique mappers <default>.
            input: use the number of reads input to the pipeline, this will include genome unmappers. But might be useful when you have additional sequence in the genome, like a transgene.
            rRNA: use the number of reads mapped to rRNA.
            uniqueXmiRNA:    use non-rRNA genomic unique mappers excluding microRNAs <for oxidized library from piRNA mutant>.
            all: use non-rRNA genomic all mappers including microRNAs.
            allXmiRNA: use non-rRNA genomic all mappers excluding microRNAs.
            miRNA: use microRNAs. normalized to: reads per millions of miRNA <for unoxidized library from piRNA mutant>.
            * Different normalization methods, including "siRNA", are available in the dual sample mode.
            * You are able to run the same library multiple times with different normalization method. They will not collapse.
    -c      Number of Threadss to use, default: $Threads
    -o      Output directory, default: current directory $PWD
    -F      A list of Fasta files, delimited by comma, used to do filtering (other than rRNA precursor sequence provide).
    -P      A list of Fasta files, delimited by comma, used to do pre-genomic mapping and analysis. For example, given "-P miniwhite.fa,virus.fa", after removing reads mappable to rRNA and miRNA hairpin, reads are mapped to miniWhite sequence first. Only the non-miniWhite-mappers are mapped to virus sequence. And only the non-miniWhite, non-virus mappers will be used in the genome mapping and further analysis.
    -O      A list of Fasta files, delimited by comma, used to do post-genomic mapping and analysis. For example, given "-O gfp.fa,luciferase.fa", after removing reads mappable rRNA, miRNA hairpin and genome, reads are mapped to gfp sequence first. Only the non-genome non-gfp mappers are mapped to luciferase sequence. If more than one sequences are put in one Fasta file, they will be treated equally. ${UNDERLINE}Please only use letters and numbers as filename and USE \$HOME instead of ~ to indicate the home directory.${RESET}${OPTIONAL}

EOF
echo -e "${COLOR_END}"
}

########
# ARGS #
########
while getopts "hi:c:o:g:vN:F:P:O:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0;;
        v)  echo2 "SmallRnaModuleVersion: v$SmallRnaModuleVersion" && exit 0;;
        i)  InputFastq=$(readlink -f "${OPTARG}"); assertFileExists "${InputFastq}";;
        o)  Outdir=$(readlink -f "${OPTARG}");;
        c)  Threads=$OPTARG;;
        g)  Genome=${OPTARG}; check_genome $Genome;;
        N)  Normmethod=$(echo "${OPTARG}" | tr '[A-Z]' '[a-z]');;
        F)  FilterMappingFileList="${OPTARG}";;
        P)  PreGenomeMappingFileList="${OPTARG}";;
        O)  PostGenomeMappingFileList="${OPTARG}";;
        *)  usage && exit 1;;
    esac
done

# if InputFastq or GENOME is undefined, print out usage and exit
if [[ -z ${InputFastq} ]]; then echo2 "Missing option -i for input fastq file" error; fi
if [[ -z ${Genome} ]]; then echo2 "Missing option -g for specifying which genome assembly to use" error; fi
declare FqName=$(basename "${InputFastq}")
declare -x Prefix=${FqName%.f[qa]*}
if [[ -z $Outdir ]]; then Outdir=$PWD; fi
if [[ -d $Outdir ]] ; then mkdir -p $Outdir; fi
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
declare -x Table=${Prefix}.basic_stats
declare PdfDir=$Outdir/pdfs && mkdir -p $PdfDir
declare ReadsDir=input_read_files && mkdir -p $ReadsDir
declare rRnaDir=rRNA_mapping && mkdir -p $rRnaDir
declare miRnaDir=hairpins_mapping && mkdir -p $miRnaDir
declare FilterDir=custom_filter && mkdir -p $FilterDir
declare PreGenomeDir=pre_genome_mapping && mkdir -p $PreGenomeDir
declare PostGenomeDir=post_genome_mapping && mkdir -p $PostGenomeDir
declare GenomeMappingDir=genome_mapping && mkdir -p $GenomeMappingDir
declare IntersectDir=intersect_genomic_features && mkdir -p $IntersectDir
declare SummaryDir=summaries && mkdir -p $SummaryDir
declare BigwigDir=bigWig_normalized_by_$Normmethod && mkdir -p $BigwigDir
declare TransposonDir=transposon_piRNAcluster_mapping_normalized_by_$Normmethod && mkdir -p $TransposonDir
declare LogDir=log && mkdir -p ${LogDir}
declare SentinelDir=sentinels && mkdir -p ${SentinelDir}

declare -i Step=1
declare RunUid=$(echo "${InputFastq}" | md5sum | cut -d" " -f1)
declare -x CommonFolder=$PIPELINE_DIRECTORY/common/$Genome
. $CommonFolder/variables
for var in rRNA_MM hairpin_MM Genome_MM transposon_MM siRNA_bot siRNA_top piRNA_bot piRNA_top; do
    if [[ -z ${!var} ]]; then echo2 "value of $var is missing, please check your $CommonFolder/variables file" error; fi
done
declare -x GenomeFa=$CommonFolder/${Genome}.fa && assertFileExists ${GenomeFa}
declare ChromSize=$CommonFolder/${Genome}.ChromInfo.txt && assertFileExists ${ChromSize}
declare -x BOWTIE_INDEXES=$CommonFolder/BowtieIndex # for Bowtie

##############################
# beginning running pipeline #
##############################
echo2 "---------------------------------------------------------------------------------"
echo2 "Beginning running [${PACKAGE_NAME}] small RNA pipeline single library mode version $SmallRnaModuleVersion"

########################################
## Pre Processing before any Mapping ###
########################################
# convering fastq to insert; quality information will be lost
echo2 "Converting fastq format into insert format"
declare Insert=$ReadsDir/${Prefix}.insert
if [[ ! -f $SentinelDir/${RunUid}.${Step}.Done || SentinelDir/${RunUid}.${Step}.Done -ot ${RunUid}.$((Step-1)).Done ]]; then
    piPipes_fastq_to_insert "${InputFastq}" ${Insert} \
    && touch $SentinelDir/${RunUid}.${Step}.Done \
    || echo2 "fq2insert failed" error
fi
let Step+=1

####################
## pre filtering ###
####################
declare CurrentMM=0
declare CurrentInput=$Insert
declare CurrentPrefix=$(basename $CurrentInput) 
declare CurrentTargetName=""
declare CurrentTargetFa=""
declare CurrentOutdir=""
# parsing customer defined pre-genomic mapping variables
if [[ ! -z $FilterMappingFileList ]]; then
    echo2 "Mapping to customer defined filtering indexes"
    eval $(echo $FilterMappingFileList | awk 'BEGIN{FS=","}{printf "export FilterMappingFiles=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}')
    for target in "${FilterMappingFiles[@]}"; do
        CurrentTargetName=$(basename $target) && CurrentTargetName=${CurrentTargetName%.fa}
        CurrentTargetFa=$(readlink -f $target)
        if [[ ! -f $CurrentTargetFa ]]; then echo2 "File $target specified by -F do not exist" error; fi
        if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetName}_filtering_mapping ]]; then
            CurrentOutdir=$FilterDir/${CurrentTargetName} && mkdir -p $CurrentOutdir || echo2 "Cannot create directory $Outdir1, please check the permission. And try to only use letter and number to name the Fasta file" error
            bowtie-build $CurrentTargetFa $CurrentOutdir/$CurrentTargetName &>/dev/null || echo2 "Failed to build the bowtie index for $targetFa" error
            echo2 "Mapping to ${CurrentTargetName}"
            bowtie -r -v 0 -a --best --strata -p $Threads -S \
                --un ${CurrentInput%.insert}.x_${CurrentTargetName}.insert \
                $CurrentOutdir/$CurrentTargetName \
                $CurrentInput \
                1> /dev/null \
                2> ${LogDir}/${CurrentPrefix}.bowtie2${CurrentTargetName}.log \
            && rm -rf $CurrentOutdir/${CurrentTargetName}*ebwt \
            && touch $SentinelDir/${RunUid}.status.${Step}.${CurrentTargetName}_filtering_mapping
        fi
        CurrentInput=${CurrentInput%.insert}.x_${CurrentTargetName}.insert
        CurrentPrefix=$(basename $CurrentInput) 
    done
fi

#####################################
# Pre Processing before any Mapping #
#####################################
# getting rid of sequences mapping to rRNA, we use -k 1 option for speed purpose
echo2 "Mapping to rRNA, with $rRNA_MM mismatch(es) allowed"
declare x_rRNA_Insert=${CurrentInput%.insert}.x_rRNA.insert
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.rRNA_mapping ]]; then
    declare -i totalReads=$(awk '{a+=$2}END{printf "%d", a}' ${CurrentInput}) \
    && echo $totalReads > ${SentinelDir}/${RunUid}.totalReads \
    && bowtie -r -S -v $rRNA_MM -k 1 -p $Threads \
        --un $x_rRNA_Insert \
        rRNA \
        ${CurrentInput} \
        1> /dev/null \
        2> ${LogDir}/${Prefix}.bowtie2rRNA.log \
    && declare -i nonrRNAReads=$(awk '{a+=$2}END{printf "%d", a}' ${x_rRNA_Insert}) \
    && echo $nonrRNAReads > $SentinelDir/${RunUid}.nonrRNAReads \
    && declare -i rRNAReads=$((totalReads-nonrRNAReads)) \
    && echo $rRNAReads > $SentinelDir/${RunUid}.rRNAReads \
    && touch $SentinelDir/${RunUid}.status.${Step}.rRNA_mapping \
    || echo2 "mapping to rRNA failed" error
fi
let Step+=1
CurrentInput=${x_rRNA_Insert}

# reading values from file, this is for resuming the job, which won't run the previous step
totalReads=$(cat $SentinelDir/${RunUid}.totalReads)
rRNAReads=$(cat $SentinelDir/${RunUid}.rRNAReads)
nonrRNAReads=$(cat $SentinelDir/${RunUid}.nonrRNAReads)

#########################
# miRNA hairpin Mapping #
#########################
echo2 "Mapping to microRNA Hairpin, with $hairpin_MM mismatch(es) allowed; only keep unique mappers"
declare x_rRNA_HairpinInsert=${CurrentInput%insert}hairpin.insert # insert file storing reads that nonmappable to rRNA and mappable to hairpin
declare x_rRNA_x_hairpinInsert=${CurrentInput%insert}x_hairpin.insert # reads that nonmappable to rRNA or hairpin
declare x_rRNA_HairpinBed2=$miRnaDir/${Prefix}.x_rRNA.hairpin.v${hairpin_MM}m1.bed2 # bed2 format with hairpin mapper, with the hairpin as reference
declare x_rRNA_HairpinLendis=$miRnaDir/${Prefix}.x_rRNA.hairpin.v${hairpin_MM}m1.lendis # length distribution for hairpin mapper
declare x_rRNA_Hairpin_GenomeBed2=$GenomeMappingDir/${Prefix}.x_rRNA.hairpin.${Genome}v${Genome_MM}a.bed2 # bed2 format with hairpin mapper, with genome as reference
declare x_rRNA_Hairpin_GenomeLog=$LogDir/${Prefix}.x_rRNA.hairpin.${Genome}v${Genome_MM}a.log # log file for hairpin mapping
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.hairpin_mapping ]]; then
    bowtie -r --norc -v $hairpin_MM -m 1 --best --strata -p $Threads -S \
        --al $x_rRNA_HairpinInsert \
        --un $x_rRNA_x_hairpinInsert \
        hairpin \
        $x_rRNA_Insert \
        2> $LogDir/${Prefix}.hairpin_unique.log \
    | samtools view -bSF 0x4 - 2>/dev/null \
    | bedtools_piPipes bamtobed -i - > ${Prefix}.x_rRNA.hairpin.v${hairpin_MM}m1.bed \
    && piPipes_insertBed_to_bed2 $x_rRNA_Insert ${Prefix}.x_rRNA.hairpin.v${hairpin_MM}m1.bed > $x_rRNA_HairpinBed2 \
    && rm -rf ${Prefix}.x_rRNA.hairpin.v${hairpin_MM}m1.bed  \
    && bed2lendis $x_rRNA_HairpinBed2 > $x_rRNA_HairpinLendis \
    && bowtie -r -v $Genome_MM -a --best --strata -p $Threads \
        -S \
        genome \
        $x_rRNA_HairpinInsert \
        2> $x_rRNA_Hairpin_GenomeLog \
    | samtools view -uS -F0x4 - 2>/dev/null \
    | bedtools_piPipes bamtobed -i - > ${Prefix}.x_rRNA.hairpin.${Genome}v${Genome_MM}a.bed \
    && piPipes_insertBed_to_bed2 $x_rRNA_HairpinInsert ${Prefix}.x_rRNA.hairpin.${Genome}v${Genome_MM}a.bed > $x_rRNA_Hairpin_GenomeBed2 \
    && rm -rf ${Prefix}.x_rRNA.hairpin.${Genome}v${Genome_MM}a.bed \
    && declare -i hairpinReads=$(bedwc $x_rRNA_HairpinBed2) \
    && echo $hairpinReads > $SentinelDir/${RunUid}.hairpinReads \
    && touch $SentinelDir/${RunUid}.status.${Step}.hairpin_mapping
fi
let Step+=1 
declare -i hairpinReads=$(cat $SentinelDir/${RunUid}.hairpinReads)

exit # TODO: work frontier

# run miRNA heterogeneity analysis
echo2 "Calculate microRNA heterogeneity"
[ ! -f .${RunUid}.status.${Step}.miRNA_pipeline ] && \
    piPipes_calculate_miRNA_heterogeneity $CommonFolder/mature2hairpin.uniq.bed  ${x_rRNA_HairpinBed2} 1> ${x_rRNA_HairpinBed2%.bed*}.sum 2> ${x_rRNA_HairpinBed2%.bed*}.hetergeneity.log
    touch .${RunUid}.status.${Step}.miRNA_pipeline
Step=$((Step+1))

#############################
# custom pre-genome mapping #
#############################
INPUT=$x_rRNA_x_hairpinInsert
MM=0 # haven't implement method to take mismatch # from user
# parsing customer defined pre-genomic mapping variables
[[ ! -z $PreGenomeMappingFileList ]] && \
    echo2 "Mapping to customer defined pre-genome mapping indexes"
    eval `echo $PreGenomeMappingFileList | awk 'BEGIN{FS=","}{printf "export PRE_GENOME_MAPPING_FILES=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}'`
    for TARGET in "${PRE_GENOME_MAPPING_FILES[@]}"; do
        targetName1=`basename $TARGET`
        targetName=${targetName1%.fa}
        targetFa=`readlink -f $TARGET`
        [[ ! -f $targetFa ]] && echo2 "File $TARGET specified by -P do not exist" error
        if [[ ! -f .${RunUid}.status.${Step}.${targetName}_mapping ]]; then
            Outdir1=$PreGenomeDir/${targetName} && mkdir -p $Outdir1 || echo2 "Cannot create directory $Outdir1, please check the permission. And try to only use letter and number to name the Fasta file" error
            bowtie-build $targetFa $Outdir1/$targetName 1>/dev/null 2>/dev/null || echo2 "Failed to build the bowtie index for $targetFa" error
            faSize -tab -detailed $targetFa > $Outdir1/${targetName}.sizes
            Prefix1=`basename $INPUT` && Prefix1=${Outdir1}/${Prefix1%.insert} && \
            echo2 "Mapping to ${targetName}" && \
            bowtie -r -v 0 -a --best --strata -p $Threads -S \
                --un ${INPUT%.insert}.x_${targetName}.insert \
                $Outdir1/$targetName \
                $INPUT \
                2> ${Prefix1}.log | \
            samtools view -bSF 0x4 - 2>/dev/null | bedtools_piPipes bamtobed -i - > ${Prefix1}.${targetName}.v${MM}a.bed && \
            piPipes_insertBed_to_bed2 $INPUT ${Prefix1}.${targetName}.v${MM}a.bed > ${Prefix1}.${targetName}.v${MM}a.bed2 && \
            rm -rf ${Prefix1}.${targetName}.v${MM}a.bed && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.bed2 -c $Outdir1/${targetName}.sizes -o $Outdir1/${targetName}.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.summary $Outdir1/ $Threads 1 1>&2 && \
            bash $DEBUG piPipes_smallRNA_bed2_to_bw.sh \
                ${Prefix1}.${targetName}.v${MM}a.bed2 \
                $Outdir1/${targetName}.sizes \
                1 \
                $Threads \
                $Outdir1 && \
            para_file=$Outdir1/${RANDOM}${RANDOM}.para && \
            echo "awk '\$3-\$2>=$siRNA_bot && \$3-\$2<=$siRNA_top' ${Prefix1}.${targetName}.v${MM}a.bed2 > ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2" >  $para_file && \
            echo "awk '\$3-\$2>=$piRNA_bot && \$3-\$2<=$piRNA_top' ${Prefix1}.${targetName}.v${MM}a.bed2 > ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2" >> $para_file && \
            ParaFly -c $para_file -Threads $Threads -failed_cmds ${para_file}.failedCommands 1>&2 && \
            rm -rf ${para_file}* && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                if ($5==1) \
                { \
                    l=$3-$2; \
                    if (l>m) m=l; \
                    if ($6=="+") s[l]+=$4;\
                    else as[l]+=$4; \
                } \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }' ${Prefix1}.${targetName}.v${MM}a.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.unique.lendis && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                if ($5==1) \
                { \
                    l=$3-$2; \
                    if (l>m) m=l; \
                    if ($6=="+") s[l]+=$4;\
                    else as[l]+=$4; \
                } \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }'  ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.piRNA.unique.lendis && \
                awk 'BEGIN{FS=OFS="\t"}\
                { \
                    if ($5==1) \
                    { \
                        l=$3-$2; \
                        if (l>m) m=l; \
                        if ($6=="+") s[l]+=$4;\
                        else as[l]+=$4; \
                    } \
                }END\
                {\
                    for (d=1;d<=m;++d) \
                    {\
                        printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                    }\
                }'  ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.siRNA.unique.lendis && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                l=$3-$2; \
                if (l>m) m=l; \
                if ($6=="+") s[l]+=$4/$5;\
                else as[l]+=$4/$5; \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }' ${Prefix1}.${targetName}.v${MM}a.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.all.lendis && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                l=$3-$2; \
                if (l>m) m=l; \
                if ($6=="+") s[l]+=$4/$5;\
                else as[l]+=$4/$5; \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }'  ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.piRNA.all.lendis && \
                awk 'BEGIN{FS=OFS="\t"}\
                { \
                    l=$3-$2; \
                    if (l>m) m=l; \
                    if ($6=="+") s[l]+=$4/$5;\
                    else as[l]+=$4/$5; \
                }END\
                {\
                    for (d=1;d<=m;++d) \
                    {\
                        printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                    }\
                }'  ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.siRNA.all.lendis && \
            piPipes_local_ping_pong -a ${Prefix1}.${targetName}.v${MM}a.bed2 -b ${Prefix1}.${targetName}.v${MM}a.bed2 -p $Threads > ${Prefix1}.${targetName}.v${MM}a.pp && \
            piPipes_local_ping_pong -a ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 -b ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 -p $Threads > ${Prefix1}.${targetName}.v${MM}a.siRNA.pp && \
            piPipes_local_ping_pong -a ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -b ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -p $Threads > ${Prefix1}.${targetName}.v${MM}a.piRNA.pp && \
            ext_len=30 && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.3end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.3end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}'  ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.3end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.5end_60.reads.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.reads.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.reads.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.3end_60.reads.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.3end_60.reads.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.3end_60.reads.percentage && \
            Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features2.R \
                $Outdir1/${Prefix}".pre-genome."${targetName}.unique_species \
                ${Prefix1}.${targetName}.v${MM}a.unique.lendis \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.unique.lendis \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.unique.lendis \
                ${ext_len} \
                ${Prefix1}.${targetName}.v${MM}a.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.3end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.3end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.3end_60.percentage 1>&2 && \
            Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features.R \
                $Outdir1/${Prefix}".pre-genome."${targetName}.all_reads \
                ${Prefix1}.${targetName}.v${MM}a.all.lendis \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.all.lendis \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.all.lendis \
                ${ext_len} \
                ${Prefix1}.${targetName}.v${MM}a.5end_60.reads.percentage \
                ${Prefix1}.${targetName}.v${MM}a.pp \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.reads.percentage \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.pp \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.reads.percentage \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.pp 1>&2 && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 -c $Outdir1/${targetName}.sizes -o /dev/stdout | awk 'BEGIN{OFS="\t"}{$1=$1"-siRNA"; print $0}' > $Outdir1/${targetName}.siRNA.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.siRNA.summary $Outdir1/siRNA $Threads 1 1>&2 && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -c $Outdir1/${targetName}.sizes -o /dev/stdout | awk 'BEGIN{OFS="\t"}{$1=$1"-piRNA"; print $0}' > $Outdir1/${targetName}.piRNA.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.piRNA.summary $Outdir1/piRNA $Threads 1 1>&2 && \
            PDFs=$Outdir1/*pdf && \
            gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/`basename ${Prefix1}`.pre-genome.${targetName}.pdf ${PDFs} && \
            rm -rf $PDFs && \
            touch .${RunUid}.status.${Step}.${targetName}_mapping
        fi
        INPUT=${INPUT%.insert}.x_${targetName}.insert
        rm -f $Outdir1/${targetName}.1.ebwt $Outdir1/${targetName}.2.ebwt $Outdir1/${targetName}.3.ebwt $Outdir1/${targetName}.4.ebwt $Outdir1/${targetName}.rev.1.ebwt $Outdir1/${targetName}.rev.2.ebwt $Outdir1/${targetName}.sizes
    done

##################
# GENOME Mapping #
##################
# take the OUTPUT of last step as INPUT
Insert=`basename ${INPUT}`
# bed2 format storing all mappers for genomic mapping
GENOME_ALLMAP_BED2=$GenomeMappingDir/${Insert%.insert}.${Genome}v${Genome_MM}.all.bed2 # all mapper in bed2 format
GENOME_ALLMAP_LOG=$GenomeMappingDir/${Insert%.insert}.${Genome}v${Genome_MM}.all.log # log file
# bed2 format storing unique mappers for genomic mapping
GENOME_UNIQUEMAP_BED2=$GenomeMappingDir/${Insert%.insert}.${Genome}v${Genome_MM}.unique.bed2
# bed2 format storing unique mappers for genomic mapping and miRNA hairpin mapper
GENOME_UNIQUEMAP_HAIRPIN_BED2=$GenomeMappingDir/${Insert%.insert}.${Genome}v${Genome_MM}.unique.+hairpin.bed2
# mapping insert file to genome
echo2 "Mapping to genome, with ${Genome_MM} mismatch(es) allowed"
[ ! -f .${RunUid}.status.${Step}.genome_mapping ] && \
    bowtie -r -v $Genome_MM -a --best --strata -p $Threads \
        --al  ${INPUT%.insert}.${Genome}v${Genome_MM}a.al.insert \
        --un  ${INPUT%.insert}.${Genome}v${Genome_MM}a.un.insert \
        -S \
        genome \
        ${INPUT} \
        2> $Genome_ALLMAP_LOG | \
    samtools view -uS -F0x4 - 2>/dev/null | \
    bedtools_piPipes bamtobed -i - > ${Insert%.insert}.${Genome}v${Genome_MM}a.insert.bed && \
    piPipes_insertBed_to_bed2 $INPUT ${Insert%.insert}.${Genome}v${Genome_MM}a.insert.bed > ${GENOME_ALLMAP_BED2} && \
    rm -rf ${Insert%.insert}.${Genome}v${Genome_MM}a.insert.bed && \
    touch .${RunUid}.status.${Step}.genome_mapping
[ ! -f .${RunUid}.status.${Step}.genome_mapping ] && echo2 "Genome mapping failed" error
Step=$((Step+1))

# separating unique and multiple mappers
echo2 "Separating unique and multiple mappers"
[ ! -f .${RunUid}.status.${Step}.separate_unique_and_multiple ] && \
    awk 'BEGIN{OFS="\t"}{if ($5==1) print $0}' ${GENOME_ALLMAP_BED2} \
    1> ${GENOME_UNIQUEMAP_BED2}    && \
    totalMapCount=`bedwc ${GENOME_ALLMAP_BED2}` && echo $totalMapCount > .${RunUid}.totalMapCount && \
    uniqueMapCount=`bedwc ${GENOME_UNIQUEMAP_BED2}` && echo $uniqueMapCount > .${RunUid}.uniqueMapCount && \
    multipMapCount=$((totalMapCount-uniqueMapCount)) && echo $multipMapCount > .${RunUid}.multipMapCount && \
    cat $x_rRNA_Hairpin_GenomeBed2 ${GENOME_UNIQUEMAP_BED2} > $Genome_UNIQUEMAP_HAIRPIN_BED2 && \
    touch .${RunUid}.status.${Step}.separate_unique_and_multiple
Step=$((Step+1))
totalMapCount=`cat .${RunUid}.totalMapCount`
uniqueMapCount=`cat .${RunUid}.uniqueMapCount`
multipMapCount=`cat .${RunUid}.multipMapCount`

#############################
# custom post-genome mapping #
#############################
INPUT=${INPUT%.insert}.${Genome}v${Genome_MM}a.un.insert
# parsing customer defined post-genomic mapping variables
[[ ! -z $PostGenomeMappingFileList ]] && \
    echo2 "Mapping to customer defined post-genome mapping indexes"
    eval `echo $PostGenomeMappingFileList | awk 'BEGIN{FS=","}{printf "export POST_GENOME_MAPPING_FILES=(" ; ;for (i=1;i<=NF;++i) printf "\"%s\" ", $i; printf ")\n";}'`
    for TARGET in "${POST_GENOME_MAPPING_FILES[@]}"; do
        targetName1=`basename $TARGET`
        targetName=${targetName1%.fa}
        targetFa=`readlink -f $TARGET`
        if [[ ! -f .${RunUid}.status.${Step}.${targetName}_mapping ]]; then
            [[ ! -f $targetFa ]] && echo2 "File $TARGET specified by -P do not exist" error
            Outdir1=$PostGenomeDir/${targetName} && mkdir -p $Outdir1 || echo2 "Cannot create directory $Outdir1, please check the permission. And try to only use letter and number to name the Fasta file" error
            bowtie-build $targetFa $Outdir1/$targetName 1>/dev/null 2>/dev/null || echo2 "Failed to build the bowtie index for $targetFa" error
            faSize -tab -detailed $targetFa > $Outdir1/${targetName}.sizes
            Prefix1=`basename $INPUT` && Prefix1=${Outdir1}/${Prefix1%.insert} && \
            echo2 "Mapping to ${targetName}" && \
            bowtie -r -v 0 -a --best --strata -p $Threads -S \
                --un ${INPUT%.insert}.x_${targetName}.insert \
                $Outdir1/$targetName \
                $INPUT \
                2> ${Prefix1}.log | \
            samtools view -bSF 0x4 - 2>/dev/null | bedtools_piPipes bamtobed -i - > ${Prefix1}.${targetName}.v${MM}a.bed && \
            piPipes_insertBed_to_bed2 $INPUT ${Prefix1}.${targetName}.v${MM}a.bed > ${Prefix1}.${targetName}.v${MM}a.bed2 && \
            rm -rf ${Prefix1}.${targetName}.v${MM}a.bed && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.bed2 -c $Outdir1/${targetName}.sizes -o $Outdir1/${targetName}.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.summary $Outdir1/ $Threads 1 1>&2 && \
            bash $DEBUG piPipes_smallRNA_bed2_to_bw.sh \
                ${Prefix1}.${targetName}.v${MM}a.bed2 \
                $Outdir1/${targetName}.sizes \
                1 \
                $Threads \
                $Outdir1 && \
            para_file=$Outdir1/${RANDOM}${RANDOM}.para && \
            echo "awk '\$3-\$2>=$siRNA_bot && \$3-\$2<=$siRNA_top' ${Prefix1}.${targetName}.v${MM}a.bed2 > ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2" >  $para_file && \
            echo "awk '\$3-\$2>=$piRNA_bot && \$3-\$2<=$piRNA_top' ${Prefix1}.${targetName}.v${MM}a.bed2 > ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2" >> $para_file && \
            ParaFly -c $para_file -Threads $Threads -failed_cmds ${para_file}.failedCommands 1>&2 && \
            rm -rf ${para_file}* && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                if ($5==1) \
                { \
                    l=$3-$2; \
                    if (l>m) m=l; \
                    if ($6=="+") s[l]+=$4;\
                    else as[l]+=$4; \
                } \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }' ${Prefix1}.${targetName}.v${MM}a.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.lendis && \
            awk 'BEGIN{FS=OFS="\t"}\
            { \
                if ($5==1) \
                { \
                    l=$3-$2; \
                    if (l>m) m=l; \
                    if ($6=="+") s[l]+=$4;\
                    else as[l]+=$4; \
                } \
            }END\
            {\
                for (d=1;d<=m;++d) \
                {\
                    printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                }\
            }'  ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.piRNA.lendis && \
                awk 'BEGIN{FS=OFS="\t"}\
                { \
                    if ($5==1) \
                    { \
                        l=$3-$2; \
                        if (l>m) m=l; \
                        if ($6=="+") s[l]+=$4;\
                        else as[l]+=$4; \
                    } \
                }END\
                {\
                    for (d=1;d<=m;++d) \
                    {\
                        printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
                    }\
                }'  ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | sort -k1,1n > ${Prefix1}.${targetName}.v${MM}a.siRNA.lendis && \
            piPipes_local_ping_pong -a ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -b ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -p $Threads > ${Prefix1}.${targetName}.v${MM}a.piRNA.pp && \
            ext_len=30 && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.bed2       | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.3end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.piRNA.3end_60.percentage && \
            awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 | bedtools_piPipes getfasta -fi $targetFa -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${Prefix1}.${targetName}.v${MM}a.siRNA.3end_60.percentage && \
            Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features2.R \
                $Outdir1/${Prefix}".post-genome."${targetName} \
                ${Prefix1}.${targetName}.v${MM}a.lendis \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.lendis \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.lendis \
                ${ext_len} \
                ${Prefix1}.${targetName}.v${MM}a.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.3end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.siRNA.3end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.5end_60.percentage \
                ${Prefix1}.${targetName}.v${MM}a.piRNA.3end_60.percentage 1>&2 && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.siRNA.bed2 -c $Outdir1/${targetName}.sizes -o /dev/stdout | awk 'BEGIN{OFS="\t"}{$1=$1"-siRNA"; print $0}' > $Outdir1/${targetName}.siRNA.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.siRNA.summary $Outdir1/ $Threads 1 1>&2 && \
            piPipes_bed2Summary -5 -i ${Prefix1}.${targetName}.v${MM}a.piRNA.bed2 -c $Outdir1/${targetName}.sizes -o /dev/stdout | awk 'BEGIN{OFS="\t"}{$1=$1"-piRNA"; print $0}' > $Outdir1/${targetName}.piRNA.summary && \
            Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $Outdir1/${targetName}.piRNA.summary $Outdir1/ $Threads 1 1>&2 && \
            PDFs=$Outdir1/*pdf && \
            gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/`basename ${Prefix1}`.post-genome.${targetName}.pdf ${PDFs} && \
            rm -rf $PDFs && \
            touch .${RunUid}.status.${Step}.${targetName}_mapping
        fi
        INPUT=${INPUT%.insert}.x_${targetName}.insert
        rm -f $Outdir1/${targetName}.1.ebwt $Outdir1/${targetName}.2.ebwt $Outdir1/${targetName}.3.ebwt $Outdir1/${targetName}.4.ebwt $Outdir1/${targetName}.rev.1.ebwt $Outdir1/${targetName}.rev.2.ebwt $Outdir1/${targetName}.sizes
    done

#####################
# Length Separation #
#####################
echo2 "Separating siRNA, piRNA based on length"
[ -z "$siRNA_bot" -o -z "$siRNA_top" ]  && echo2 "length for siRNA is not defined! please check the \"variable\" file under common\$Genome" error
[ -z "$piRNA_bot" -o -z "$piRNA_top" ]  && echo2 "lengt for piRNA is not defined! please check the \"variable\" file under common\$Genome" error
[ ! -f .${RunUid}.status.${Step}.sep_length ] && \
    para_file=${RANDOM}${RANDOM}.para && \
    echo "awk '\$3-\$2>=$siRNA_bot && \$3-\$2<=$siRNA_top' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}siRNA.bed2" > $para_file && \
    echo "awk '\$3-\$2>=$piRNA_bot && \$3-\$2<=$piRNA_top' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}piRNA.bed2" >> $para_file && \
    ParaFly -c $para_file -Threads $Threads -failed_cmds ${para_file}.failedCommands 1>&2 && \
    rm -rf ${para_file}* && \
    touch  .${RunUid}.status.${Step}.sep_length
[ ! -f .${RunUid}.status.${Step}.sep_length ] && "separating siRNA, piRNA failed"
Step=$((Step+1))

# plotting length distribution
echo2 "Plotting length distribution"
[ ! -f .${RunUid}.status.${Step}.plotting_length_dis ] && \
    awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${GENOME_ALLMAP_BED2}  | sort -k1,1n > ${GENOME_ALLMAP_BED2}.lendis && \
    awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${GENOME_UNIQUEMAP_BED2}  | sort -k1,1n > ${GENOME_UNIQUEMAP_BED2}.lendis && \
    Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R ${GENOME_ALLMAP_BED2}.lendis $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.x_hairpin 1>&2 && \
    Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R ${GENOME_UNIQUEMAP_BED2}.lendis $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.x_hairpin 1>&2 && \
    awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${GENOME_ALLMAP_BED2}.lendis $x_rRNA_HairpinLendis | sort -k1,1n > ${GENOME_ALLMAP_BED2}.+hairpin.lendis && \
    awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${GENOME_UNIQUEMAP_BED2}.lendis $x_rRNA_HairpinLendis | sort -k1,1n > ${GENOME_UNIQUEMAP_BED2}.+hairpin.lendis && \
    Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R ${GENOME_ALLMAP_BED2}.+hairpin.lendis $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.+hairpin 1>&2 && \
    Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_lendis.R ${GENOME_UNIQUEMAP_BED2}.+hairpin.lendis $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.+hairpin 1>&2 && \
    touch .${RunUid}.status.${Step}.plotting_length_dis
Step=$((Step+1))

##################
# Print to table #
##################
# change dual library mode normalization method if change here
echo -e "total reads as input of the pipeline\t${totalReads}" > $Table && \
echo -e "rRNA reads with ${rRNA_MM} mismatches\t${rRNAReads}" >> $Table && \
echo -e "miRNA hairpin reads\t${hairpinReads}" >> $Table && \
echo -e "genome mapping reads (-rRNA; +miRNA_hairpin)\t$((totalMapCount+hairpinReads))" >> $Table && \
echo -e "genome mapping reads (-rRNA; -miRNA_hairpin)\t${totalMapCount}" >> $Table && \
echo -e "genome unique mapping reads (-rRNA; +miRNA_hairpin)\t$((uniqueMapCount+hairpinReads))" >> $Table && \
echo -e "genome unique mapping reads (-rRNA; -miRNA_hairpin)\t${uniqueMapCount}" >> $Table && \
echo -e "genome multiple mapping reads (-rRNA; -miRNA_hairpin)\t${multipMapCount}" >> $Table && \
export TOTAL_GENOME_MAPPING_READS=$((totalMapCount+hairpinReads))

# normalization method
# input | rRNA | unique | uniqueXmiRNA | all | allXmiRNA | miRNA
case "$Normmethod" in
input)
    NormScale=`head -1 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
rrna)
    NormScale=`head -2 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
mirna)
    NormScale=`head -3 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
all)
    NormScale=`head -4 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
allxmirna)
    NormScale=`head -5 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
unique)
    NormScale=`head -6 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
uniquexmirna)
    NormScale=`head -7 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
*)
    echo2 "unrecognized normalization option: $Normmethod; using the default method" "warning"
    NormScale=`head -6 $Table | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
;;
esac
echo $NormScale > .depth

####################################
# Intersecting with GENOME Feature #
####################################
echo2 "Intersecting with genomic features, make length distribution, nucleotide fraction for siRNA/piRNA assigned to each feature"
[ ! -f .${RunUid}.status.${Step}.intersect_with_genomic_features ] && \
bash $DEBUG piPipes_intersect_smallRNA_with_genomic_features.sh \
    ${GENOME_ALLMAP_BED2} \
    $SummaryDir/`basename ${GENOME_ALLMAP_BED2%.bed2}` \
    $Threads \
    $IntersectDir \
    1>&2 && \
    touch .${RunUid}.status.${Step}.intersect_with_genomic_features
Step=$((Step+1))

#######################
# Making BigWig Files #
#######################
# make BW files
echo2 "Making bigWig files for genome browser"
[ ! -f .${RunUid}.status.${Step}.make_bigWig_normalized_by_$Normmethod ] && \
    bash $DEBUG piPipes_smallRNA_bed2_to_bw.sh \
        ${GENOME_ALLMAP_BED2} \
        ${ChromSize} \
        ${NormScale} \
        $Threads \
        $BigwigDir && \
    bash $DEBUG piPipes_smallRNA_bed2_to_bw.sh \
        ${GENOME_ALLMAP_BED2%bed2}piRNA.bed2 \
        ${ChromSize} \
        ${NormScale} \
        $Threads \
        $BigwigDir && \
    touch .${RunUid}.status.${Step}.make_bigWig_normalized_by_$Normmethod
Step=$((Step+1))

##############################################
# Direct mapping to transposon/piRNA cluster #
##############################################
echo2 "Direct mapping to transposon and piRNA cluster and make distribution plot"
. $CommonFolder/genomic_features
Insert=`basename $x_rRNA_Insert`
[ ! -f .${RunUid}.status.${Step}.direct_mapping_normalized_by_$Normmethod ] && \
for t in "${DIRECT_MAPPING[@]}"; do \
    bowtie -r -v ${transposon_MM} -a --best --strata -p $Threads \
        -S \
        ${t} \
        ${x_rRNA_Insert} \
        2> ${TransposonDir}/${t}.log | \
    samtools view -uS -F0x4 - 2>/dev/null | \
    samtools sort -o -@ $Threads - foo | \
    bedtools_piPipes bamtobed -i - > $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.insert.bed && \
    piPipes_insertBed_to_bed2 $x_rRNA_Insert $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.insert.bed > $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.insert.bed2 && \
    piPipes_bed2Summary -5 -i $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.insert.bed2 -c $CommonFolder/BowtieIndex/${t}.sizes -o $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.summary && \
    Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.summary $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.normalized_by_$Normmethod $Threads $NormScale 1>&2 && \
    PDFs=$TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.normalized_by_${Normmethod}*pdf && \
    gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/${Insert%.insert}.${t}.pdf ${PDFs} && \
    rm -rf $PDFs && \
    rm -rf $TransposonDir/${Insert%.insert}.${t}.a${transposon_MM}.insert.bed
done && \
touch .${RunUid}.status.${Step}.direct_mapping_normalized_by_$Normmethod
Step=$((Step+1))

#####################################################
# Direct mapping to and quantification with eXpress #
#####################################################
# for accurate quantification, we map to the index of gene+cluster+repBase.
# echo2 "Quantification by direct mapping and eXpress"
# [ ! -f .${RunUid}.status.${Step}.direct_mapping_no_normalization ] && \
# awk '{for (j=0;j<$2;++j) print $1}' $x_rRNA_x_hairpinInsert | \
# bowtie \
#     -r \
#     -v ${transposon_MM} \
#     -a --best --strata \
#     -p $Threads \
#     -S \
#     gene+cluster+repBase \
#     - 2> $EXPRESS_DIR/${Prefix}.bowtie.gene+cluster+repBase.bowtie.log | \
#     samtools view -bS - > \
#     $EXPRESS_DIR/${Prefix}.bowtie.gene+cluster+repBase.bam && \
# touch .${RunUid}.status.${Step}.direct_mapping_no_normalization
# Step=$((Step+1))

# deprecated
# [ ! -f .${RunUid}.status.${Step}.quantification_by_eXpress ] && \
# express \
#     -B $eXpressBATCH \
#     -m $(( (siRNA + piRNA_top)/2 )) \
#     -s $(( (piRNA_top - 18)/2 )) \
#     --output-align-prob \
#     -o $EXPRESS_DIR \
#     --no-update-check \
#     $CommonFolder/${Genome}.gene+cluster+repBase.fa \
#     $EXPRESS_DIR/${Prefix}.bowtie.gene+cluster+repBase.bam \
#     1>&2 2> $EXPRESS_DIR/${Prefix}.eXpress.log && \
# touch .${RunUid}.status.${Step}.quantification_by_eXpress
# Step=$((Step+1))

################
# Joining Pdfs #
################
echo2 "Merging pdfs"
if [ -f $PdfDir/${Prefix}.features.pdf ]
then
    [ ! -f .${RunUid}.status.${Step}.merge_pdfs ] && \
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/${Prefix}.${PACKAGE_NAME}.small_RNA_pipeline.${PROG_VERSION}.pdf \
            $PdfDir/${Prefix}.pie.pdf \
            $PdfDir/${Prefix}.siRNA.pie.pdf \
            $PdfDir/${Prefix}.piRNA.pie.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.x_hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.x_hairpin.lendis.pdf  \
            $PdfDir/${Prefix}.features.pdf  && \
        rm -rf $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.x_hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.x_hairpin.lendis.pdf  \
            $PdfDir/${Prefix}.pie.pdf \
            $PdfDir/${Prefix}.siRNA.pie.pdf \
            $PdfDir/${Prefix}.piRNA.pie.pdf \
            $PdfDir/${Prefix}.features.pdf && \
        touch .${RunUid}.status.${Step}.merge_pdfs
else
    [ ! -f .${RunUid}.status.${Step}.merge_pdfs ] && \
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/${Prefix}.${PACKAGE_NAME}.small_RNA_pipeline.${PROG_VERSION}.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.x_hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.x_hairpin.lendis.pdf  && \
        rm -rf $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.+hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_UNIQUEMAP_BED2}`.x_hairpin.lendis.pdf \
            $PdfDir/`basename ${GENOME_ALLMAP_BED2}`.x_hairpin.lendis.pdf  && \
        touch .${RunUid}.status.${Step}.merge_pdfs
fi
Step=$((Step+1))

#############
# finishing #
#############
rm -f $GenomeMappingDir/*bed2
rm -f $TransposonDir/*bed2
echo2 "Finished running ${PACKAGE_NAME} small RNA pipeline version $SmallRnaModuleVersion"
echo2 "---------------------------------------------------------------------------------"
touch .${Genome}.${PROG_VERSION}
