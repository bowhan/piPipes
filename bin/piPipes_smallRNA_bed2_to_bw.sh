
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


# this script takes 
input_bed2=$1 # input bed2 format
input_bed_name=$(basename $input_bed2)
input_bed_name_prefix=${input_bed_name%bed2}
chrom_size_file=$2 # file with chrom size information
norm_scale=$3 # normalization scale
output_dir=$4 && mkdir -p $output_dir # output dir
para_file=$output_dir/${RANDOM}${RANDOM}${RANDOM}.make_bigWig.para
touch $para_file
# sorting is required for bedtools genomecov; using sort supporing MT
( sort \
    --parallel=$Threads --temporary-directory=$output_dir \
    -k1,1 $input_bed2 | tee $output_dir/${input_bed_name_prefix}sorted.bed2 | awk '$5==1' > $output_dir/${input_bed_name_prefix}sorted.uniq.bed2 ) \
|| \
( sort                                               \
    -k1,1 $input_bed2 | tee $output_dir/${input_bed_name_prefix}sorted.bed2 | awk '$5==1' > $output_dir/${input_bed_name_prefix}sorted.uniq.bed2 ) \
|| echo2 "Failed to sort $input_bed2 file" error

# making bedgraph
piPipes_bed2_to_bedGraph -i $output_dir/${input_bed_name_prefix}sorted.bed2      -c $chrom_size_file -p $Threads  -o $output_dir/${input_bed_name_prefix}sorted.
piPipes_bed2_to_bedGraph -i $output_dir/${input_bed_name_prefix}sorted.uniq.bed2 -c $chrom_size_file -p $Threads  -o $output_dir/${input_bed_name_prefix}sorted.uniq.
# making bigWig
if [[ -s $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph ]]; then
    echo -e "awk -v depth=$norm_scale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph > $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph.normalized && bedGraphToBigWig $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph.normalized $chrom_size_file $output_dir/${input_bed_name_prefix}Watson.bigWig" >> $para_file
    echo -e "awk -v depth=$norm_scale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $output_dir/${input_bed_name_prefix}sorted.uniq.Watson.bedGraph > $output_dir/${input_bed_name_prefix}sorted.uniq.Watson.bedGraph.normalized && bedGraphToBigWig $output_dir/${input_bed_name_prefix}sorted.uniq.Watson.bedGraph.normalized $chrom_size_file $output_dir/${input_bed_name_prefix}uniq.Watson.bigWig" >> $para_file
fi
if [[ -s $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph ]]; then
    echo -e "awk -v depth=$norm_scale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph  > $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph.normalized && bedGraphToBigWig $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph.normalized $chrom_size_file $output_dir/${input_bed_name_prefix}Crick.bigWig" >> $para_file
    echo -e "awk -v depth=$norm_scale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $output_dir/${input_bed_name_prefix}sorted.uniq.Crick.bedGraph  > $output_dir/${input_bed_name_prefix}sorted.uniq.Crick.bedGraph.normalized && bedGraphToBigWig $output_dir/${input_bed_name_prefix}sorted.uniq.Crick.bedGraph.normalized $chrom_size_file $output_dir/${input_bed_name_prefix}uniq.Crick.bigWig" >> $para_file
fi

if [[ -s $para_file ]]; then
    ParaFly -c $para_file -CPU $Threads -failed_cmds ${para_file}.failedCommands &>/dev/null \
    || echo2 "Failed to generate bigWig files using ParaFly" error
fi

rm -rf $para_file ${para_file}".completed" \
    $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph \
    $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph \
    $output_dir/${input_bed_name_prefix}sorted.uniq.Watson.bedGraph  \
    $output_dir/${input_bed_name_prefix}sorted.uniq.Crick.bedGraph \
    $output_dir/${input_bed_name_prefix}sorted.Watson.bedGraph.normalized \
    $output_dir/${input_bed_name_prefix}sorted.Crick.bedGraph.normalized \
    $output_dir/${input_bed_name_prefix}sorted.uniq.Watson.bedGraph.normalized \
    $output_dir/${input_bed_name_prefix}sorted.uniq.Crick.bedGraph.normalized \
    $output_dir/${input_bed_name_prefix}sorted.bed2  \
    $output_dir/${input_bed_name_prefix}sorted.uniq.bed2