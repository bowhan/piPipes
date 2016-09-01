#! /bin/bash
set -o pipefail
[[ -z ${PIPELINE_DIRECTORY} ]] \
    && declare PIPELINE_DIRECTORY=$(dirname $PWD) # if used outside of piPipes
declare MYSQLBIN=${PIPELINE_DIRECTORY}/bin/mysql

[[ ! -f ${1}.refSeq.Genes.bed12 ]] \
    && mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${1} < ${MYSQLBIN}/RefSeqGene.sql \
    | awk 'BEGIN{getline}{ \
        for(i=1;i<=10;++i) printf "%s\t", $i; \
        split($11,a,","); \
        split($12,b,","); \
        for(i=1;i<=$10;++i) printf "%d,", a[i]-b[i]; \
        printf "\t"; \
        for(i=1;i<=$10;++i) printf "%d,", b[i]-$2; \
        printf "\n"; \
        }' \
    > ${1}.refSeq.Genes.bed12

[[ ! -f ${1}.refSeq.Exons.bed6 ]] \
    && bedtools bed12tobed6 -i ${1}.refSeq.Genes.bed12 \
    | awk 'BEGIN{FS=OFS="\t"}{$4=$4".exon"c[$4]++; print $0}' \
    > ${1}.refSeq.Exons.bed6

[[ ! -f ${1}.refSeq.5UTR.bed6 ]] \
    && mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${1} < ${MYSQLBIN}/RefSeqFivePrimeUTR.sql \
    | awk 'BEGIN{FS=OFS="\t";getline}{ \
        $4=$4"."5UTR; \
        print $0; \
        }' \
    > ${1}.refSeq.5UTR.bed6

[[ ! -f ${1}.refSeq.3UTR.bed6 ]] \
    && mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D ${1} < ${MYSQLBIN}/RefSeqThreePrimeUTR.sql \
    | awk 'BEGIN{FS=OFS="\t";getline}{ \
        $4=$4"."3UTR; \
        print $0; \
        }' \
    > ${1}.refSeq.3UTR.bed6
