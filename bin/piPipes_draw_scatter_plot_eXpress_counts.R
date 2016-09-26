
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


source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/piPipes.R",sep=""))

library(pacman)
p_load(gdata)
p_load(ggplot2)
p_load(grid)
p_load(ggthemes)
p_load(scales)

argv = commandArgs (TRUE)
sample1 = read.table (argv[1], T)
sample2 = read.table (argv[2], T)
name1 = argv[3]
name2 = argv[4]

main = basename (argv[5])
main = gsub ("\\."," ",main)
main = paste (strwrap(main, width = 80), collapse = "\n")

sample = merge (sample1, sample2, by="target_id")
lim = roundUp (10*(max (sample$eff_counts.x, sample$eff_counts.y)))/10
# for fly dm6 genome, NM and NR annotation is current unavailable; so flybase annotation is used here. FBtr includes both mRNA and ncRNA...
Group = ifelse ( startsWith(sample[,1],"NM_") | startsWith(sample[,1],"FBtr") , "NM", ifelse (startsWith(sample[,1], "NR_"), "NR", "Transposon"))
sample_trn = sample[Group == "Transposon",]

write.csv (sample, file=paste (argv[5],".csv", sep=""))

pdf (paste (argv[5],"1.pdf", sep=""), onefile=TRUE, width=10, height=10, title=main )
ggplot( sample, aes(x = eff_counts.x, y = eff_counts.y, colour = Group, size = Group) ) + 
    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
    theme_few() + 
    scale_colour_few ("dark") + 
    geom_point(alpha=0.85, na.rm=T) +
    scale_x_log10 ( limits = c(0.1,lim), 
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x) ) ) +
    scale_y_log10 ( limits = c(0.1,lim),
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x) ) ) +
    annotation_logticks () +
	xlab ( substitute( paste( italic(name1), "normalized number of reads (log10)"), list(name1=name1) )) +
	ylab ( substitute( paste( italic(name2), "normalized number of reads (log10)"), list(name2=name2) )) +
    coord_fixed()
g = dev.off ()

pdf (paste (argv[5],"2.pdf", sep=""), onefile=TRUE, width=10, height=10, title=main )
ggplot() + 
    theme_few () + 
    scale_colour_few ("dark") + 
    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
    geom_point (data = sample_trn, aes ( x = eff_counts.x, y = eff_counts.y , size = 3, alpha=0.85, na.rm=T ))  + 
    geom_density2d (data = sample, aes ( x = eff_counts.x, y = eff_counts.y , size = 0.1, colour = Group , geom='polygon', alpha=0.85 )) +
    scale_x_log10 ( limits = c(0.1,lim), 
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x) ) ) +
    scale_y_log10 ( limits = c(0.1,lim),
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x) ) ) +
    annotation_logticks () +
    xlab ( substitute( paste( italic(name1), "normalized number of reads (log10)"), list(name1=name1) )) +
    ylab ( substitute( paste( italic(name2), "normalized number of reads (log10)"), list(name2=name2) )) +
    coord_fixed()
g = dev.off ()


