#! /bin/bash

##########
# Config #
##########
declare -xr PIPELINE_DIRECTORY=$(dirname `readlink -f $0`)
declare -xr MYBIN=${PIPELINE_DIRECTORY}/bin
declare -xr PATH=${PIPELINE_DIRECTORY}/bin:$PATH
. $PIPELINE_DIRECTORY/bin/colors.sh
. $PIPELINE_DIRECTORY/bin/piPipes_bash_functions.sh
. $PIPELINE_DIRECTORY/bin/usage.sh
declare -xr CONTACT_EMAILS=${UNDERLINE}"piPipesQ@gmail.com"$RESET

#######
# Run #
#######
if [ $# -lt 1 ]; then usage; exit 1; fi
PROGRAM=`echo ${1} | tr '[A-Z]' '[a-z]'`
case $PROGRAM in
i|install)
    shift && bash  piPipes_install_genomes.sh "$@" ;;
s|small|smallrna|smallrna-seq|sra)
    shift && bash  piPipes_smallRNA.sh "$@" ;;
s2|small2|smallrna2|smallrna-seq2|sra2)
    shift && bash  piPipes_smallRNA2.sh "$@" ;;
r|rnaseq|rna-seq|rsq|rna)
    shift && bash  piPipes_RNASeq.sh "$@" ;;
r2|rnaseq2|rna-seq2|rsq2|rna2)
    shift && bash  piPipes_RNASeq2.sh "$@" ;;
d|deg|degseq|deg-seq|degradome|race|cage|cage-seq|cageseq)
    shift && bash  piPipes_DegradomeSeq.sh "$@" ;;
c|chip|chipseq|chip-seq)
    shift && bash  piPipes_ChIPSeq.sh "$@" ;;
c2|chip2|chipseq2|chip-seq2)
    shift && bash  piPipes_ChIPSeq2.sh "$@" ;;
g|genomic|genome|genome-seq|genomeseq|dna|dnaseq|dna-seq)
    shift && bash  piPipes_GenomeSeq.sh "$@" ;;
*)
    echo2 "unrecognized option \"${1}\"! \nplease type \"${PACKAGE_NAME}\" without options to see all the options and usage." "error" ;;
esac
