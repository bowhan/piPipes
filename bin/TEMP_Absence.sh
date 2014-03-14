#!/bin/bash -x
# TEMP (Transposable Element Movement present in a Population)
# 2013-06-14
# Jiali Zhuang(jiali.zhuang@umassmed.edu)
# Zhiping Weng Lab
# Programs in Bioinformatics and Integrative Biology
# University of Massachusetts Medical School

#usage function
usage() {
echo -en "\e[1;36m"
cat <<EOF

usage: $0 -i input_file.sorted.bam -s scripts_directory -o output_directory -r transposon_rpmk.bed -t reference.2bit -f fragment_size -c CPUs -h 

TEMP is a software package for detecting transposable elements (TEs) 
insertions and excisions from pooled high-throughput sequencing data. 
Please send questions, suggestions and bug reports to:
jiali.zhuang@umassmed.edu

Options:
        -i     Input file in bam format with full path. Please sort and index the file before calling this program. 
               Sorting and indexing can be done by 'samtools sort' and 'samtools index'
        -s     Directory where all the scripts are
        -o     Path to output directory. Default is current directory
        -r     Annotated transposon positions in the genome (e.g., repeakMask) in bed6 format with full path
        -t     2bit file for the reference genome (can be downloaded from UCSC Genome Browser)
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -c     An integer specifying the number of CUPs used. Default is 4
        -h     Show help message

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hi:c:f:o:r:s:t:" OPTION
do
        case $OPTION in
                h)
                        usage && exit 1
		;;
                i)
                        BAM=$OPTARG
		;;
	        f)
		        INSERT=$OPTARG
		;;
                o)
                        OUTDIR=$OPTARG
                ;;
                c)
                        CPU=$OPTARG
                ;;
                s)
                        BINDIR=$OPTARG
                ;;
	        r)
		        TEBED=$OPTARG
		;;
	        t)
		        REF=$OPTARG
		;;
                ?)
                        usage && exit 1
                ;;
        esac
done

if [[ -z $BAM ]] || [[ -z $BINDIR ]] || [[ -z $TEBED ]] || [[ -z $REF ]]
then
        usage && exit 1
fi
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=4
[ ! -z "${INSERT##*[!0-9]*}" ] || INSERT=500
[ ! -z $OUTDIR ]  || OUTDIR=$PWD

mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
touch ${OUTDIR}/.writting_permission && rm -rf ${OUTDIR}/.writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

function checkExist {
        echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
        [[ $? != 0 ]] && echo -e "\e[1;36mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && exit 1
}
echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo"
checkExist "rm"
checkExist "mkdir"
checkExist "date"
checkExist "mv"
checkExist "sort"
checkExist "touch"
checkExist "awk"
checkExist "grep"
checkExist "bwa"
checkExist "samtools"
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

name=`basename $BAM`
i=${name/.sorted.bam/}
echo $name
echo $i
if [[ ! -s $name ]]
then
    cp $BAM ./
fi
if [[ ! -s $name.bai ]]
then cp $BAM.bai ./
fi

#Detect excision sites
samtools view -XF 0x2 $name > $i.unpair.sam
awk -F "\t" '{OFS="\t"; if ($9 != 0) print $0}' $i.unpair.sam > temp1.sam
perl $BINDIR/pickUniqIntervalPos.pl temp1.sam $INSERT > $i.unproper.uniq.interval.bed

rm temp1.sam $i.unpair.sam

# Map to transposons
bedtools intersect -a $TEBED -b $i.unproper.uniq.interval.bed -f 1.0 -wo > temp
perl $BINDIR/filterFalsePositive.ex.pl temp $INSERT $i.final.pairs.rpmk.bed
bedtools intersect -a $TEBED -b $i.final.pairs.rpmk.bed -f 1.0 -wo > temp2

perl $BINDIR/excision.clustering.pl temp2 $i.excision.cluster.rpmk
rm temp temp2 $i.unproper.uniq.interval.bed $i.final.pairs.rpmk.bed

# Identify breakpoints using soft-clipping information
perl $BINDIR/pickSoftClipping.over.pl $i.excision.cluster.rpmk $REF > $i.excision.cluster.rpmk.sfcp
perl $BINDIR/refine_breakpoint.ex.pl

# Estimate excision sites frequencies
perl $BINDIR/pickOverlapPair.ex.pl $i.excision.cluster.rpmk.refined.bp > $i.excision.cluster.rpmk.refined.bp.refsup
perl $BINDIR/summarize_excision.pl
