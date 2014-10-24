
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

while getopts "l:r:L:R:i:I:o:x:" OPTION; do
    case $OPTION in
        l)	dLEFT_IP_FQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
        r)	dRIGHT_IP_FQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
        L)	dLEFT_INPUT_FQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
        R)	dRIGHT_INPUT_FQ=`readlink -f $OPTARG`; PE_MODE=1 ;;
        i)	dIP_FQ=`readlink -f $OPTARG`; SE_MODE=1 ;;
        I)	dINPUT_FQ=`readlink -f $OPTARG`; SE_MODE=1 ;;
        o)	dDIRECTMAPPING_DIR=`readlink -f $OPTARG` ;;
        x)  dDIRECTMAPPING_INDEX=$OPTARG ;;
        *)	exit 1 ;;
    esac
done

TRANSCRIPTOME_SIZES=$BOWTIE2_INDEXES/${dDIRECTMAPPING_INDEX}.sizes
[[ ! -f $TRANSCRIPTOME_SIZES ]] && echo2 "Cannot find size file TRANSCRIPTOME_SIZES" "error"

if [[ -n "dLEFT_IP_FQ" ]]; then
# paired-end mode
    echo2 "Direct mapping IP to $dDIRECTMAPPING_INDEX"
    bowtie2 \
        -x $dDIRECTMAPPING_INDEX \
        -1 ${dLEFT_IP_FQ} \
        -2 ${dRIGHT_IP_FQ} \
        -q \
        $bowtie2PhredOption \
        -a \
        -X 800 \
        --no-mixed \
        --quiet \
        -p $CPU \
        2> ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.log | \
    samtools view -bS - > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.bam && \
    samtools sort -o -@ $CPU ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.bam ${dDIRECTMAPPING_DIR}/foo | \
    bedtools_piPipes bamtobed -i - | \
    awk -v MAPQ=10 '$5 > MAPQ' > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bed

    echo2 "Direct mapping INPUT to $dDIRECTMAPPING_INDEX"
    bowtie2 \
        -x $dDIRECTMAPPING_INDEX \
        -1 ${dLEFT_INPUT_FQ} \
        -2 ${dRIGHT_INPUT_FQ} \
        -q \
        $bowtie2PhredOption \
        -a \
        -X 800 \
        --no-mixed \
        --quiet \
        -p $CPU \
        2> ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.log | \
    samtools view -bS - > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.bam && \
    samtools sort -o -@ $CPU ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.bam ${dDIRECTMAPPING_DIR}/foo | \
    bedtools_piPipes bamtobed -i - | \
    awk -v MAPQ=10 '$5 > MAPQ' > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bed
else
# single-end mode
    echo2 "Direct mapping IP to $dDIRECTMAPPING_INDEX"
    bowtie2 \
        -x $dDIRECTMAPPING_INDEX \
        -U ${dIP_FQ} \
        -q \
        $bowtie2PhredOption \
        -a \
        --quiet \
        -p $CPU \
        2> ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.log | \
    samtools view -bS - > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.bam && \
    samtools sort -o -@ $CPU ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.bam ${dDIRECTMAPPING_DIR}/foo | \
    bedtools_piPipes bamtobed -i - | \
    awk -v MAPQ=10 '$5 > MAPQ' > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bed

    echo2 "Direct mapping INPUT to $dDIRECTMAPPING_INDEX"
    bowtie2 \
        -x $dDIRECTMAPPING_INDEX \
        -U ${dINPUT_FQ} \
        -q \
        $bowtie2PhredOption \
        -a \
        --quiet \
        -p $CPU \
        2> ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.log | \
    samtools view -bS - > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.bam && \
    samtools sort -o -@ $CPU ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.bam ${dDIRECTMAPPING_DIR}/foo | \
    bedtools_piPipes bamtobed -i - | \
    awk -v MAPQ=10 '$5 > MAPQ' > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bed
fi

echo2 "Making summary graph"
grep -v 'NM_' $TRANSCRIPTOME_SIZES | grep -v 'NR_' > ${dDIRECTMAPPING_DIR}/transposon.sizes && \
bedtools_piPipes genomecov -i ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bed    -g $TRANSCRIPTOME_SIZES -bg -scale $NormScaleIP    > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bedGraph && \
bedtools_piPipes genomecov -i ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bed -g $TRANSCRIPTOME_SIZES -bg -scale $NormScaleINPUT > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bedGraph && \
bedGraphToBigWig ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bedGraph    $TRANSCRIPTOME_SIZES ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig && \
bedGraphToBigWig ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bedGraph $TRANSCRIPTOME_SIZES ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig && \
rm -f ${dDIRECTMAPPING_DIR}/*.bedGraph && \
paraFile=${dDIRECTMAPPING_DIR}/${RANDOM}${RANDOM}.para && \
bgP=${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig && \
bgM=${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig && \
awk -v bgP=${bgP} -v bgM=${bgM} -v binSize=${BINSIZE} '{print "bigWigSummary", bgP, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x027 >", bgP"."$1; print "bigWigSummary", bgM, $1, 0, $2, binSize, "| sed -e \x27s/n\\/a/0/g\x027 >", bgM"."$1;}' ${dDIRECTMAPPING_DIR}/transposon.sizes > $paraFile && \
ParaFly -c $paraFile -CPU $CPU && \
paraFile=${OUTDIR}/drawFigures && \
rm -f ${dDIRECTMAPPING_DIR}/${PREFIX}.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.summary
for i in `cut -f1 ${dDIRECTMAPPING_DIR}/transposon.sizes`; do \
    [ ! -s ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i  ] && echo | awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i
    [ ! -s ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i ] && echo | awk -v binSize=${BINSIZE} 'BEGIN{for (i=0;i<binSize-1;++i){printf "%d\t", 0} print 0;}' > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i
    awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i    > ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t
    awk -v name=$i '{for (i=1;i<=NF;++i){printf "%s\t%d\t%d\n", name, i, $i}}' ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i > ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t
    paste ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t | cut -f1,2,3,6 >> ${dDIRECTMAPPING_DIR}/${PREFIX}.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.summary
    rm -rf ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i \
           ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i \
           ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t \
           ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.$i.t
done

Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R ${dDIRECTMAPPING_DIR}/${PREFIX}.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.summary ${dDIRECTMAPPING_DIR}/${PREFIX}.${dDIRECTMAPPING_INDEX}.sorted.unique.bigWig.summary $CPU $NormScale 1>&2 && \
PDFs=${dDIRECTMAPPING_DIR}/${PREFIX}.${dDIRECTMAPPING_INDEX}.sorted*pdf && \
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.${dDIRECTMAPPING_INDEX}.unique.pdf ${PDFs} && \
rm -rf ${PDFs}


echo2 "quantification IP using eXpress"
express \
    -B $eXpressBATCH \
    -o $dDIRECTMAPPING_DIR \
    --library-size $EFFECTIVE_DEPTH_IP \
    --no-update-check \
    $COMMON_FOLDER/${GENOME}.${dDIRECTMAPPING_INDEX}.fa \
    ${dDIRECTMAPPING_DIR}/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.bam 1>&2 && \
awk -v depth=$NormScale 'BEGIN{OFS="\t"; getline; print}{$8*=depth; print}' $dDIRECTMAPPING_DIR/results.xprs > $dDIRECTMAPPING_DIR/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.results.xprs.normalized
    
echo2 "quantification input using eXpress"
express \
    -B $eXpressBATCH \
    -o $dDIRECTMAPPING_DIR \
    --library-size $EFFECTIVE_DEPTH_INPUT \
    --no-update-check \
    $COMMON_FOLDER/${GENOME}.${dDIRECTMAPPING_INDEX}.fa \
    ${dDIRECTMAPPING_DIR}/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.bam 1>&2 && \
awk -v depth=$NormScale 'BEGIN{OFS="\t"; getline; print}{$8*=depth; print}' $dDIRECTMAPPING_DIR/results.xprs > $dDIRECTMAPPING_DIR/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.results.xprs.normalized

echo -e "target_id\teff_counts" > $dDIRECTMAPPING_DIR/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts && \
echo -e "target_id\teff_counts" > $dDIRECTMAPPING_DIR/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts && \
cat $dDIRECTMAPPING_DIR/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.results.xprs.normalized    | cut -f2,8 | grep -v id | sort -k1,1 | bedtools_piPipes groupby -i stdin -g 1 -c 2 -o mean >> $dDIRECTMAPPING_DIR/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts && \
cat $dDIRECTMAPPING_DIR/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.results.xprs.normalized | cut -f2,8 | grep -v id | sort -k1,1 | bedtools_piPipes groupby -i stdin -g 1 -c 2 -o mean >> $dDIRECTMAPPING_DIR/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts && \
Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_scatter_plot_eXpress_counts.R \
    $dDIRECTMAPPING_DIR/${PREFIX}.IP.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts \
    $dDIRECTMAPPING_DIR/${PREFIX}.input.${dDIRECTMAPPING_INDEX}.results.xprs.normalized_counts \
    ${PREFIX}.IP \
    ${PREFIX}.input \
    $PDF_DIR/${PREFIX}.IP.vs.input.effective_counts
