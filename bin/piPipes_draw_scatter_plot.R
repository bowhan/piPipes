
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
pkgTest ("gdata")
pkgTest ("ggplot2")
pkgTest ("grid")
pkgTest ("ggthemes")
pkgTest ("scales")
pkgTest ("gridExtra")

argv = commandArgs (TRUE)
sample1 = read.table (argv[1])
sample2 = read.table (argv[2])
name1 = argv[3]
name2 = argv[4]

main = basename (argv[5])
main = gsub ("\\."," ",main)
main = paste (strwrap(main, width = 80), collapse = "\n")

pdf (paste (argv[5],".pdf", sep=""), title=main )

if (ncol(sample1)==4) { # with grouping information 
        colnames (sample1) = c("name", "group", "S1", "AS1")
        colnames (sample2) = c("name", "group", "S2", "AS2")
        sample = merge (sample1, sample2, by="name")

		lim = roundUp (10*(max (sample$S1, sample$S2, sample$AS1, sample$AS2)))/10

		gg=ggplot( sample, aes(x = S1, y = AS1, color= factor (group.x)) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
				) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.75, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    guides(colour = guide_legend (title=expression (paste (italic("Li., Cell, 2009"), " Transposon group")), title.position = "top")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name1), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)
		
		gg=ggplot( sample, aes(x = S2, y = AS2, color= factor (group.x)) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.75, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    guides (colour = guide_legend (title=expression (paste (italic("Li., Cell, 2009"), " Transposon group")), title.position = "top")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name2) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)
		
		gg=ggplot( sample, aes(x = S1, y = S2, color= factor (group.x)) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.75, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    guides(colour = guide_legend (title=expression (paste (italic("Li., Cell, 2009"), " Transposon group")), title.position = "top")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)
		
		gg=ggplot( sample, aes(x = AS1, y = AS2, color= factor (group.x)) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.75, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    guides(colour = guide_legend (title=expression (paste (italic("Li., Cell, 2009"), " Transposon group")), title.position = "top")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
			print (gg)

} else {
        colnames (sample1) = c("name", "S1", "AS1")
        colnames (sample2) = c("name", "S2", "AS2")
        sample = merge (sample1, sample2, by="name")
		lim = roundUp (10*(max (sample$S1, sample$S2, sample$AS1, sample$AS2)))/10
		
		gg=ggplot( sample, aes(x = S1, y = AS1) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
				) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.5, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name1), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)

		gg=ggplot( sample, aes(x = S2, y = AS2) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.5, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name2) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)

		gg=ggplot( sample, aes(x = S1, y = S2) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.5, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  sense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)

		gg=ggplot( sample, aes(x = AS1, y = AS2) ) + 
			theme (
				plot.margin=unit(c(1,1,0,0),"lines"),
				axis.text=element_text (size=4), 
				axis.title=element_text(size=6), 
				legend.margin=unit(0,"lines"), 
				panel.margin=unit(0, "lines"), 
				axis.ticks.margin=unit(0,"lines"),
				legend.key.size=unit(0.5,"lines")
		 		) +
		    scale_size_manual( values=c(0,8) ) + 
		    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') + 
		    theme_few () + 
		    scale_fill_continuous(guide = "legend") + 
		    geom_point(size=5, alpha=0.5, na.rm=T) +
		    scale_colour_manual(values=c("lightblue","black","darkgreen","red")) + 
		    scale_x_log10 ( limits = c(1,lim), 
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    scale_y_log10 ( limits = c(1,lim),
		                    breaks = trans_breaks("log10", function(x) 10^x),
		                    labels = trans_format("log10", math_format(10^.x) ) ) +
		    annotation_logticks () +
		    xlab ( substitute ( paste(italic(name1) ,"  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    ylab ( substitute ( paste(italic(name2), "  antisense, normalized number of reads (log10)"), list(name1=name1, name2=name2))) +
		    coord_fixed()
		print (gg)
}

gc = dev.off()
write.csv (sample, file=paste (argv[5],".csv", sep=""))



