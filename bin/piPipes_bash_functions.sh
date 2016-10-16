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
declare -x ISO_8601='%Y-%m-%d %H:%M:%S %Z'
declare -xi DEFAULT_THREADS=8

function echo2 {
if [[ $# -eq 1 ]]; then
    echo -e $COLOR_GREEN"[`date "+$ISO_8601"`] $1${COLOR_END}";
else
    case $2 in
        error)      echo -e $COLOR_RED_BOLD"[`date "+$ISO_8601"`] Error: $1${COLOR_END}" && exit 1 ;;
        warning)    echo -e $COLOR_MAGENTA"[`date "+$ISO_8601"`] Warning: $1${COLOR_END}" ;;
        *)          echo -e $COLOR_GREEN"[`date "+$ISO_8601"`] $1${COLOR_END}";;
    esac
fi
}
export -f echo2

function assertBinExists {
    if ! which $1 &>/dev/null; then
        echo2 "Required program \"$1\" is not available" error 
    fi
}
export -f assertBinExists

function assertFileExists {
    if [[ ! -f "${1}" ]]; then echo2 "Required file \"${1}\" doesn't exist" error; fi
}
export -f assertFileExists

function assertDirExists {
    if [[ ! -d "${1}" ]]; then echo2 "Required directory \"${1}\" doesn't exist" error; fi
}
export -f assertDirExists

function assertDirWritable {
    if [[ ! -w "${1}" ]]; then echo2 "Directory \"${1}\" is not writable" error; fi
}
export -f assertDirWritable

function fileExists {
    if [[ ! -f "${1}" ]]; then return 1; fi
    return 0
}
export -f fileExists

function dirExists {
    if [[ ! -d "${1}" ]]; then return 1; fi
    return 0
}
export -f dirExists

function bowtieCheck {
    for sfx in "1.ebwt" "2.ebwt" "3.ebwt" "4.ebwt" "rev.1.ebwt" "rev.2.ebwt"; do
        if [[ ! -f ${1}.$sfx ]]; then return 1; fi
    done
    return 0
}
export -f bowtieCheck

function bowtie2Check {
    for sfx in "1.bt2" "2.bt2" "3.bt2" "4.bt2" "rev.1.bt2" "rev.2.bt2"; do
        if [[ ! -f ${1}.$sfx ]]; then return 1; fi
    done
    return 0
}
export -f bowtie2Check

function bwaCheck {
    for sfx in "amb" "ann" "bwt" "pac" "sa"; do
        if [[ ! -f ${1}.$sfx ]]; then return 1; fi
    done
    return 0
}
export -f bwaCheck

function StarCheck {
    for f in "chrLength.txt" "chrNameLength.txt" "chrName.txt" "chrStart.txt" "Genome" "genomeParameters.txt" "SA" "SAindex" "sjdbInfo.txt" "sjdbList.out.tab"; do
        if [[ ! -f ${1}/${f} ]]; then return 1; fi
    done
    return 0
}
export -f StarCheck

# check whether the genome has been installed
function check_genome {
	for genome_supported in $(cat $PIPELINE_DIRECTORY/common/genome_supported.txt); do
		if [[ $genome_supported == $1 ]]; then return; fi
	done
    echo2 "The genome \"${1}\" is not currently supported. Currently supported genomes: $(cat $PIPELINE_DIRECTORY/common/genome_supported.txt); Please use \"piPipes install\" to install new genome" error
}
export -f check_genome

# count reads for bed2 format
function bedwc {
    awk '{a[$7]=$4}END{COUNTER=0; for (b in a) {COUNTER+=a[b]} printf "%d" , COUNTER}' $1
}
export -f bedwc

# produce length distribution for bed2 format
function bed2lendis {
    awk '{l=$3-$2; if (l>m) m=l; c[l]+=$4/$5;}END{for (d=1;d<=m;++d) {printf "%d\t%.0f\n", d, (c[d]?c[d]:0)} }' $1
}
export -f bed2lendis

# merge bed2
# the reason of +=$4 is because the same ($1,$2,$3) combination could appear more than once if mismatch is allowed
function mergeMirBed2 {
    awk -v d1=$3 -v d2=$4 'BEGIN{OFS="\t";} \
        { \
            if (ARGIND==1) { \
                mir[$1]=1; \
                ct1[$1"_"$2"_"$3"_"$5]+=$4; \
            } else { \
                mir[$1]=1;\
                ct2[$1"_"$2"_"$3"_"$5]+=$4;\
            }\
        }\
        END {\
                for (m in mir) {\
                    for (x=-6;x<7;++x) {\
                        for (y=-6;y<7;++y) {\
                            print m, x, y, ct1[m"_"x"_"y"_"1]?(ct1[m"_"x"_"y"_"1]*d1):0, ct2[m"_"x"_"y"_"1]?(ct2[m"_"x"_"y"_"1]*d2):0, ct1[m"_"x"_"y"_"3]?(ct1[m"_"x"_"y"_"3]*d1):0, ct2[m"_"x"_"y"_"3]?(ct2[m"_"x"_"y"_"3]*d2):0;\
                        }\
                    }\
                }\
        }' $1  $2
}
export -f mergeMirBed2
