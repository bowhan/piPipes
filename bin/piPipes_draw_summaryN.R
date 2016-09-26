
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

p_load("parallel")
p_load("RColorBrewer")
p_load("ggplot2")
p_load("grid")

# https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

my_colors = addalpha(brewer.pal(7,"Set1"), 0.75)

draw_summaryN = function (p, pdfPrefix, normScales, name1, name2) {
	pdf (paste (pdfPrefix, p[1,1], ".pdf", sep=""))
	par (bty="n")
	sample_size = ncol(p) - 2
	sample_size_ = length(normScales)
	stopifnot(sample_size==sample_size_)
	for (i in 3:ncol(p)) {
		p[,i] = p[,i] * normScales[i-2]
	}
	plot (p$V2,p$V3, xlim=c(0,nrow(p)), ylim=c(0, 1.2*max(p[,c(3,4,5,6)])) , type='n', xlab=paste("Gene body", nrow(p), sep=" "), ylab="Signal", tck=0.01, main=p[1,1])
	for (i in seq(3, ncol(p)-1, 2)) {
		points (p$V2, p[,i], col=my_colors[i-2], type="s", lwd=3)
	}
	for (i in seq(4, ncol(p), 2)) {
		points (p$V2, p[,i], col=my_colors[i-2], type="s", lwd=1)
	}
	
	# TODO: only support N == 2
	legend("topright", inset=.05, title="sample", c(paste(name1, "IP"),paste(name1, "input"),paste(name2, "IP"), paste(name2, "input")), fill=my_colors[1:4], border=FALSE, bty="n") 
	gc = dev.off()
	
	
	# g = ggplot(p[,-1], aes(x=p[,2])) +
	# theme_few () +
	# scale_colour_few() +
	# theme( panel.border = element_blank () ,
	# 	panel.grid.major=element_blank(),
	# 	panel.grid.minor=element_blank(),
	# 	plot.title=element_text(family="Helvetica", lineheight=.8) ) +
	# ggtitle("ChIP") +
	# geom_line( size=1, alpha=0.75 ) +
	# xlab ("Position (bp)") +
	# ylab("ChIP-seq signal")
	# 
	# for (i in seq(3, ncol(p)-1, 2)) {
	# 	g + geom_polygon(aes(y=p[,3], fill=my_colors[i]), size=1) 
	# }
	# for (i in seq(4, ncol(p), 2)) {
	# 	g + geom_points (aes(y = p[,i]), fill=my_colors[i])
	# }
}

argv  = commandArgs (TRUE)
summaryTableFile = argv[1]
outPdfPrefix = argv[2]
numOfCore = argv[3]
normScales = as.double(unlist(strsplit(argv[4], ",")))
name1 = argv[5]
name2 = argv[6]

# normScale = as.double (argv[4])
summaryTable = read.table (summaryTableFile, F)
summaryTableSplited = split (summaryTable, summaryTable$V1)
mclapply (summaryTableSplited, draw_summaryN, mc.cores=numOfCore, pdfPrefix=outPdfPrefix, normScales=normScales, name1 = name1, name2 = name2)
