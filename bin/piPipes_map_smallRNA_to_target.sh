
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

# this script uses exported variables so will NOT work alone 

echo2 "Building index for ${CurrentTargetNamePrefix}"
bowtie-build \
    $CurrentTarget \
    $BowtieIndexName \
    &> ${LogDir}/bowtie_build.${CurrentTargetNamePrefix}.log \
|| echo2 "Failed to build the bowtie index for $CurrentTargetNamePrefix" error

faSize -tab -detailed $CurrentTarget > $CurrentOutdir/${CurrentTargetNamePrefix}.sizes \
|| echo2 "Failed to run faSize to obtain the size of ${CurrentTargetName}" error

bowtie -r -v $CurrentMM \
    -a --best --strata -p $Threads -S \
    --un ${NextInput} \
    $BowtieIndexName \
    $CurrentInput \
    2> ${LogDir}/${CurrentInputNamePrefix}.bowtie2${CurrentTargetNamePrefix}.log \
| samtools view -bSF 0x4 - \
    2> /dev/null \
| bedtools_piPipes bamtobed -i - \
    > ${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}.bed \
&& piPipes_insertBed_to_bed2 \
    $CurrentInput \
    ${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}.bed \
    > $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.bed2 \
&& awk -v l=$siRNA_bot -v h=$siRNA_top '$3-$2>=l && $3-$2<=h' \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.bed2 \
    > $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.bed2 \
&& awk -v l=$piRNA_bot -v h=$piRNA_top '$3-$2>=l && $3-$2<=h' \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.bed2 \
    > $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.bed2 \
&& rm -f ${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}.bed \
&& rm -f $CurrentOutdir/${CurrentTargetName}*ebwt 

mkdir -p $PdfDir/${CurrentTargetNamePrefix}
declare bed=$CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.bed2
piPipes_bed2Summary \
    -5 \
    -i $bed \
    -c $CurrentOutdir/${CurrentTargetNamePrefix}.sizes \
    -o ${bed%.bed*}.summary \
&& Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_summary.R \
    ${bed%.bed*}.summary \
    $PdfDir/$CurrentTargetNamePrefix/ \
    $Threads \
    1 \
    &> /dev/null \
&& piPipes_smallRNA_bed2_to_bw.sh \
    $bed \
    $CurrentOutdir/${CurrentTargetNamePrefix}.sizes \
    1 \
    $Threads \
    $BigwigDir \
|| echo2 "failed to generate summary plot for ${CurrentTargetNamePrefix}" error

for bed in \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.bed2 \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.bed2 \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.bed2 \
; do
    echo2 "Calculating length distribution for ${CurrentTargetNamePrefix}"

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
    }' $bed \
    | sort -k1,1n \
    > ${bed%.bed2}.unique.lendis

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
    }' $bed \
    | sort -k1,1n \
    > ${bed%.bed2}.all.lendis

    echo2 "Calculating cis-Ping-Pong score for ${CurrentTargetNamePrefix}"

    piPipes_local_ping_pong \
        -a ${bed} \
        -b ${bed} \
        -p $Threads \
        > ${bed%bed2}pp
    
    # 5' end
    ## species
    echo2 "Calculating nucleotide composition for ${CurrentTargetNamePrefix}"
    awk -v ext_len=$NucCompExtLen \
    'BEGIN{OFS="\t"} { \
        if (($5==1)&&(!printed[$7])) { \
            printed[$7]=1; \
            if ($2>=ext_len) { \
                for (i=1;i<=1;++i) { \
                    if ($6=="+") { \
                        print $1, $2-ext_len,   $2+ext_len+1,$4,$5,$6; \
                    } else { \
                        print $1, $3-ext_len-1, $3+ext_len,  $4,$5,$6; \
                    } \
                } \
            } \
        } \
    }'  ${bed} \
    | bedtools_piPipes getfasta \
        -fi $CurrentTarget \
        -bed stdin \
        -fo stdout \
        -s \
        -name \
        -tab \
    | piPipes_nuc_percentage.py \
        $NucCompExtLen \
    > ${bed%.bed*}.5end_$((NucCompExtLen*2)).species.percentage
    ## reads
    awk -v ext_len=$NucCompExtLen 'BEGIN{OFS="\t"} { \
        if (($5==1)&&(!printed[$7])) {\
            printed[$7]=1; if ($2>=ext_len) { \
                for (i=1;i<=$4;++i) { \
                    if ($6=="+") { \
                        print $1,$2-ext_len,$2+ext_len+1, $4, $5, $6; \
                    } else { \
                        print $1,$3-ext_len-1,$3+ext_len, $4, $5, $6; \
                    }\
                }\
            }\
        }\
    }' ${bed} \
    | bedtools_piPipes getfasta \
        -fi $CurrentTarget \
        -bed stdin \
        -fo stdout \
        -s \
        -name \
        -tab \
    | piPipes_nuc_percentage.py \
        $NucCompExtLen \
    > ${bed%.bed*}.5end_$((NucCompExtLen*2)).reads.percentage

    # 3' end
    ## species 
    awk -v ext_len=$NucCompExtLen \
    'BEGIN{OFS="\t"} { \
        if (($5==1)&&(!printed[$7])) {\
            printed[$7]=1; \
            if ($2>=ext_len) { \
                for (i=1;i<=1;++i) { \
                    if ($6=="-") { \
                        print $1, $2-ext_len,   $2+ext_len+1, $4, $5, $6; \
                        } else { \
                        print $1, $3-ext_len-1, $3+ext_len,   $4, $5, $6; \
                    } \
                } \
            } \
        } \
    }'  ${bed} \
    | bedtools_piPipes getfasta \
        -fi $CurrentTarget \
        -bed stdin \
        -fo stdout \
        -s \
        -name \
        -tab \
    | piPipes_nuc_percentage.py \
        $NucCompExtLen \
    > ${bed%.bed*}.3end_$((NucCompExtLen*2)).species.percentage
    ## reads 
        awk -v ext_len=$NucCompExtLen \
    'BEGIN{OFS="\t"} { \
        if (($5==1)&&(!printed[$7])) {\
            printed[$7]=1; \
            if ($2>=ext_len) { \
                for (i=1;i<=$4;++i) { \
                    if ($6=="-") { \
                        print $1, $2-ext_len,   $2+ext_len+1, $4, $5, $6; \
                        } else { \
                        print $1, $3-ext_len-1, $3+ext_len,   $4, $5, $6; \
                    } \
                } \
            } \
        } \
    }'  ${bed} \
    | bedtools_piPipes getfasta \
        -fi $CurrentTarget \
        -bed stdin \
        -fo stdout \
        -s \
        -name \
        -tab \
    | piPipes_nuc_percentage.py \
        $NucCompExtLen \
    > ${bed%.bed*}.3end_$((NucCompExtLen*2)).reads.percentage
done

Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features2.R \
    $PdfDir/${CurrentTargetNamePrefix}/0.unique_species \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.unique.lendis \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.unique.lendis \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.unique.lendis \
    ${NucCompExtLen} \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.5end_$((NucCompExtLen*2)).species.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.3end_$((NucCompExtLen*2)).species.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.5end_$((NucCompExtLen*2)).species.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.3end_$((NucCompExtLen*2)).species.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.5end_$((NucCompExtLen*2)).species.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.3end_$((NucCompExtLen*2)).species.percentage \
|| echo2 "failed to draw unique/species pdf for ${CurrentTargetNamePrefix}" error

Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features2.R \
    $PdfDir/${CurrentTargetNamePrefix}/1.all_length \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.all.lendis \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.all.lendis \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.all.lendis \
    ${NucCompExtLen} \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.5end_$((NucCompExtLen*2)).reads.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.3end_$((NucCompExtLen*2)).reads.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.5end_$((NucCompExtLen*2)).reads.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.siRNA.3end_$((NucCompExtLen*2)).reads.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.5end_$((NucCompExtLen*2)).reads.percentage \
    $CurrentOutdir/${CurrentInputNamePrefix}.${CurrentTargetNamePrefix}_v${CurrentMM}a.piRNA.3end_$((NucCompExtLen*2)).reads.percentage \
|| echo2 "failed to draw all/reads pdf for ${CurrentTargetNamePrefix}" error

PDFs=$PdfDir/${CurrentTargetNamePrefix}/*pdf \
&& gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite \
    -sOutputFile=$PdfDir/PreGenome.${CurrentTargetNamePrefix}.pdf \
    ${PDFs} \
|| echo2 "failed to combine all the pdf for ${CurrentTargetNamePrefix}" error