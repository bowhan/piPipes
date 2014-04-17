#!/bin/bash
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

usage: $0 -i input_file.sorted.bam -s scripts_directory -o output_directory -r transposon_database.fa -t annotated_TEs.bed -m MISMATCH -f fragment_size -c CPUs -h 

TEMP is a software package for detecting transposable elements (TEs) 
insertions and excisions from pooled high-throughput sequencing data. 
Please send questions, suggestions and bug reports to:
jiali.zhuang@umassmed.edu

Options:
        -i     Input file in bam format with full path. Please sort and index the file before calling this program. 
               Sorting and indexing can be done by 'samtools sort' and 'samtools index'
        -s     Directory where all the scripts are
        -o     Path to output directory. Default is current directory
        -r     Transposon sequence database in fasta format with full path
        -t     Annotated TEs in BED6 format with full path. Detectied snsertions that overlap with annoated TEs will be filtered. 
        -m     Number of mismatch allowed when mapping to TE concensus sequences
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -c     An integer specifying the number of CUPs used. Default is 8
        -h     Show help message

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hi:c:f:m:o:r:s:t:" OPTION
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
	        m)
		        MM=$OPTARG
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
		        TESEQ=$OPTARG
		;;
	        t)
                        ANNO=$OPTARG
                ;;
                ?)
                        usage && exit 1
                ;;
        esac
done

if [[ -z $BAM ]] || [[ -z $BINDIR ]] || [[ -z $TESEQ ]] || [[ -z $ANNO ]]
then
        usage && exit 1
fi
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "${INSERT##*[!0-9]*}" ] || INSERT=500
[ ! -z "${MM##*[!0-9]*}" ] || MM=3
[ ! -z $OUTDIR ]  || OUTDIR=$PWD

mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
touch .writting_permission && rm -rf .writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

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

cp $TESEQ ./
name=`basename $BAM`
te=`basename $TESEQ`
i=${name/.sorted.bam/}
echo $name
echo $i
if [[ ! -s $name ]]
then
    ln -s $BAM ./	# piPipes
fi
if [[ ! -s $name.bai ]]
then ln -s $BAM.bai ./	# piPipes
fi

# Get the mate seq of the uniq-unpaired reads
samtools view -XF 0x2  $name > $i.unpair.sam
perl $BINDIR/pickUniqPairFastq.pl $i.unpair.sam $i.unpair.uniq
perl $BINDIR/pickUniqPos.pl $i.unpair.sam > $i.unpair.uniq.bed

# Map to transposons
bwa index -a is $te
bwa aln -t $CPU -n $MM -l 100 -R 1000 $te $i.unpair.uniq.1.fastq > $i.unpair.uniq.1.sai
bwa aln -t $CPU -n $MM -l 100 -R 1000 $te $i.unpair.uniq.2.fastq > $i.unpair.uniq.2.sai
bwa sampe -P $te $i.unpair.uniq.1.sai $i.unpair.uniq.2.sai $i.unpair.uniq.1.fastq $i.unpair.uniq.2.fastq > $i.unpair.uniq.transposons.sam

#Summary
samtools view -hSXF 0x2 $i.unpair.uniq.transposons.sam > $i.unpair.uniq.transposons.unpair.sam
perl $BINDIR/pickUniqMate.pl $i.unpair.uniq.transposons.unpair.sam $i.unpair.uniq.bed > $i.unpair.uniq.transposons.bed
cp $i.unpair.uniq.transposons.bed $i.unpair.uniq.transposons.filtered.bed

# Throw out false positives (TEs that overlap with annotated TEs)
#perl $BINDIR/filterFalsePositive.in.pl $i.unpair.uniq.transposons.bed $ANNO > $i.unpair.uniq.transposons.fp.bed
#if [[ -s $i.unpair.uniq.transposons.fp.bed ]]
#then
#    bedtools subtract -a $i.unpair.uniq.transposons.bed -b $i.unpair.uniq.transposons.fp.bed -f 1.0 > $i.unpair.uniq.transposons.filtered.bed
#else 
#    cp $i.unpair.uniq.transposons.bed $i.unpair.uniq.transposons.filtered.bed
#fi

#Prepare for insertion breakpoints identification
awk -F "\t" -v sample=$i '{OFS="\t"; print $1,$2,$3,sample,$5,$6}' $i.unpair.uniq.transposons.filtered.bed >> tmp
perl $BINDIR/mergeTagsWithoutGap.pl tmp > $i.uniq.transposons.filtered.woGap.bed
perl $BINDIR/mergeTagsWithGap.pl $i.uniq.transposons.filtered.woGap.bed $INSERT > $i.uniq.transposons.filtered.wGap.bed
rm tmp
perl $BINDIR/get_class.pl $i.uniq.transposons.filtered.wGap.bed $i > $i.uniq.transposons.filtered.wGap.class.bed
perl $BINDIR/make.bp.bed.pl $i.uniq.transposons.filtered.wGap.class.bed

rm $i.unpair.sam $i.unpair.uniq.bed $i.unpair.uniq.?.fastq $i.unpair.uniq.?.sai 
rm $i.unpair.uniq.transposons.sam $i.unpair.uniq.transposons.unpair.sam $i.uniq.transposons.filtered.woGap.bed $i.uniq.transposons.filtered.wGap.bed

#Detect insertion breakpoints using soft-clipping information
samtools view -bf 0x2 $name > $i.pair.bam
samtools index $i.pair.bam
perl $BINDIR/pickClippedFastq.pl $i $te
perl $BINDIR/refine_breakpoint.in.pl

rm $i.pair.bam $i.pair.bam.bai

#Estimate insertion frequencies
perl $BINDIR/pickOverlapPair.in.pl $i.insertion.refined.bp $INSERT > $i.insertion.refined.bp.summary


#Remove called sites that overlap with annotated TEs
awk -F "\t" '{OFS="\t"; if ($5=="antisense") $5="-"; if ($5=="sense") $5="+"; if ($1 ~ /^chr/) print $1,$2,$3,$4,".",$5}' $i.insertion.refined.bp.summary > tmp
bedtools intersect -a tmp -b $ANNO -f 0.1 -wo -s > tmp1
awk -F "\t" '{OFS="\t"; if (($4==$10)&&($6==$12)) print $1,$2,$3,$4,$5,$6}' tmp1 > tmp2
if [[ -s "tmp2" ]]
then
    awk -F "\t" '{OFS="\t"; if ($1 ~ /^chr/) print}' $i.insertion.refined.bp.summary > tmp1
    bedtools subtract -a tmp1 -b tmp2 -f 1.0 > tmp3
    head -n 1 $i.insertion.refined.bp.summary > tmp4
    cat tmp4 tmp3 > $i.insertion.refined.bp.summary
fi
rm tmp*

################################
##End of processing insertions##
################################

