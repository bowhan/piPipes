
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

.libPaths(c(paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep='/'), .libPaths()))

# this function is modified from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE, lib=paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep="/"), repos='http://cran.us.r-project.org')
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
}

# for bioconductor package test and installation
bioConductorTest = function (x) {
	if (!require(x,character.only = TRUE)) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(x, lib=paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep="/") )
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
}

# copied from http://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
	if(length(x) != 1) stop("'x' must be of length 1")
	10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# function to draw small RNA ggplot lendis
draw_smRNA_lendis = function (file, main) {
	lendis = read.table (file,FALSE)
	minRow = min ( min(lendis[lendis[,2]!=0,1]), min(lendis[lendis[,3]!=0,1]) )
	maxRow = max ( max(lendis[lendis[,2]!=0,1]), max(lendis[lendis[,3]!=0,1]) )
	lendis = lendis[seq(minRow, maxRow),]
	lendis[,3]=lendis[,3]*-1
	ru=roundUp( (max(lendis$V2)-min(lendis$V3) )/20 ) 
	gg = ggplot (lendis, aes (V1,V2)) + 
			theme_minimal() + 
			geom_bar (stat="identity", colour="blue", fill="blue") +
			geom_bar (aes(V1,V3), stat="identity", colour="red", fill="red") +
			scale_x_discrete (breaks=c(minRow:maxRow)) + 
			coord_cartesian(xlim = c(minRow, maxRow)) + 
			scale_y_continuous(labels = comma, breaks=seq(ru*-20, ru*20, ru)) +
			labs(title=paste("Length Distribution",main,sep="\t")) + 
			xlab("Length, nt") + 
			ylab("Reads") + 
			theme(title=element_text(size=6, colour='black'), plot.margin=unit(c(1,1,0,0),"lines"), legend.margin=unit(0,"lines"), panel.margin=unit(0, "lines"), axis.ticks.margin=unit(0,"lines"), axis.text=element_text (size=4), axis.title=element_text(size=6), axis.ticks = element_line(size = 0.5))
	return (gg)
}

# function to draw ping pong
draw_ping_pong = function (ppbedfile, main) {
	ppbed = read.table (ppbedfile, F)
	minRow = min(ppbed$V1)
	maxRow = max(ppbed$V1)
	zScore=(ppbed[10,2]-mean(ppbed[-10,2]))/sd(ppbed[-10,2])
	gg = ggplot (ppbed, aes (V1,V2)) + 
			theme_minimal() + 
			geom_bar (stat="identity") + 
			scale_x_discrete (breaks=c(minRow:maxRow)) + 
			coord_cartesian(xlim = c(minRow, maxRow)) + 
			scale_y_continuous(labels = comma, breaks=seq(0,max(ppbed$V2),roundUp(max(ppbed$V2)/10))) +
			labs(title=paste("5' to 5' Overlap",main,paste("Z=",signif (zScore,3),sep=""), sep="\t")) + 
			xlab("Length, nt") + 
			ylab("Reads x Reads (scores)") + 
			theme(title=element_text(size=7, colour='black',family="Helvetica"), plot.margin=unit(c(1,1,0,0),"lines"), legend.margin=unit(0,"lines"), panel.margin=unit(0, "lines"), axis.ticks.margin=unit(0,"lines"), legend.key.size=unit(0.5,"lines")) +
			theme(legend.title=element_blank(), legend.position = "bottom", axis.text=element_text (size=5,family="Helvetica"), axis.text.x=element_text (angle=25, size=5,family="Helvetica"), legend.text=element_text(size=5,family="Helvetica"), axis.title=element_text(size=6), axis.ticks = element_line(size = 0.5))
	return (gg)
}

# function to draw small RNA ggplot percentage
draw_smRNA_percentage = function (file, ext, main) {
	t = read.table (file, FALSE )
	colnames (t) = c ("A","C","G","T")
	rownames (t) = seq (-1*ext,ext,1)
	tm = melt(cbind(t,pos=rownames(t)),is.var = c('bind'))
	levs = levels (tm$variable)
	levs = sort (levs, T)
	tm$variable = factor(tm$variable, levels=levs)
	tm2 = tm[with (tm, order (tm$variable)),]
	tm2$pos = factor (tm2$pos, levels=seq (-1*ext,ext,1))
	gg = ggplot(tm2,aes(x = pos, y = value, fill=variable)) + 
			theme_minimal() + 
			geom_bar(position = "fill", stat="identity") + 
			scale_y_continuous(labels = percent_format()) + 
			labs(title=paste("Nucleotide Percentage", main, sep="\t")) + 
			xlab("Relative Position") + 
			ylab("Percentage of Nucleotides") +
			scale_fill_manual(values=c("red","darkgreen", "yellow","blue")) +
			theme(title=element_text(size=6, colour='black',family="Helvetica"), plot.margin=unit(c(1,1,0,0),"lines"), legend.margin=unit(0,"lines"), panel.margin=unit(0, "lines"), axis.ticks.margin=unit(0,"lines"), legend.key.size=unit(0.5,"lines")) +
			theme(legend.title=element_blank(), legend.position = "bottom", axis.text=element_text (size=4,family="Helvetica"), axis.text.x=element_text (angle=25, size=5,family="Helvetica"), legend.text=element_text(size=5,family="Helvetica"), axis.title=element_text(size=6), axis.ticks = element_line(size = 0.5))
	return (gg)
}

# function to draw gene mode from a single table
draw_summary = function (p, pdfPrefix, normScale) {
	pdf (paste (pdfPrefix, p[1,1], ".pdf", sep=""))
	par (bty="n")
	p$V3 = p$V3 * normScale
	p$V4 = p$V4 * normScale
	plot (p$V2,p$V3, xlim=c(0,nrow(p)), ylim=c(1.4*min(p$V4), 1.4*max(p$V3)) , type='n', xlab=paste("scaled gene body, scaled to", nrow(p), sep=" "), ylab="normalized signal", tck=0.01, main=p[1,1])
	points (p$V2, p$V3, col="blue", type="s")
	points (p$V2, p$V4, col="red", type="s")
	abline (h=0, lty=2)
	gc = dev.off()
}

# function to draw balloon plot
draw_microRNA_balloon = function (t1, hetName, mutName, outDir) {
	sum4th = sum(t1[,4])
	if (sum4th > 0) { t1[,4] = 100*t1[,4]/sum(t1[,4]); }
	sum5th = sum(t1[,5])
	if (sum5th > 0) { t1[,5] = 100*t1[,5]/sum(t1[,5]); }
	sum6th = sum(t1[,6])
	if (sum6th > 0) { t1[,6] = 100*t1[,6]/sum(t1[,6]); }
	sum7th = sum(t1[,7])
	if (sum7th > 0) { t1[,7] = 100*t1[,7]/sum(t1[,7]); }
	
	fivePrimeArm_het =  xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V3)
	fivePrimeArm_mut =  xtabs (as.integer (t1$V5) ~ t1$V2 + t1$V3)
	threePrimeArm_het = xtabs (as.integer (t1$V6) ~ t1$V2 + t1$V3)
	threePrimeArm_mut = xtabs (as.integer (t1$V7) ~ t1$V2 + t1$V3)
	
	hetName=gsub ("\\."," ",hetName)
	mutName=gsub ("\\."," ",mutName)
	pdf (paste (outDir, '/', hetName, mutName, t1$V1[1], ".miRNAballoonPlot.pdf", sep=''), family="Helvetica")
	par (mfrow=c(2,2),mar=c(5,2,2,1))
	
	main = paste (t1$V1[1], hetName," 5' arm:", sum4th, sep=' ')
	main = paste (strwrap(main, width = 35), collapse = "\n") 
	balloonplot (fivePrimeArm_het,  main=main,xlab="5'",ylab="3'",sorted=T,label.size=.6,text.size=.5,rowmar=1,show.zeros=F)
	
	main = paste (t1$V1[1], mutName," 5' arm:", sum5th, sep=' ')
	main = paste (strwrap(main, width = 35), collapse = "\n") 
	balloonplot (fivePrimeArm_mut,  main=main,xlab="5'",ylab="3'",sorted=T,label.size=.6,text.size=.5,rowmar=1,show.zeros=F, dotcolor="lightgreen")
	
	main = paste (t1$V1[1], hetName," 3' arm:", sum6th, sep=' ')
	main = paste (strwrap(main, width = 35), collapse = "\n") 
	balloonplot (threePrimeArm_het, main=main,xlab="5'",ylab="3'",sorted=T,label.size=.6,text.size=.5,rowmar=1,show.zeros=F)
	
	main = paste (t1$V1[1], mutName," 3' arm:", sum7th, sep=' ')
	main = paste (strwrap(main, width = 35), collapse = "\n") 
	balloonplot (threePrimeArm_mut, main=main,xlab="5'",ylab="3'",sorted=T,label.size=.6,text.size=.5,rowmar=1,show.zeros=F, dotcolor="lightgreen")
	
	invisible(dev.off())
}

# function to write aggregate plot
draw_agg = function (t1, name) {
	colors = c("dodgerblue3","darkorchid4","tomato1","cyan4","darkolivegreen4","deeppink4","thistle4","yellow4")
	plots = read.table(t1, F, sep="\t")
	colnames(plots) = c('Feature','ChIP','Position','Signal')
	ggplot(plots, aes(x=Position,y=Signal,color=ChIP)) + 
		ggtitle(name) + 
		theme( plot.title=element_text(lineheight=.8, face="bold" ) ) + 
		geom_line( size=1 ) + 
		theme_minimal() + 
		scale_color_manual(values= colors[1:length(levels (plots$ChIP))]) + 
		ylab("ChIP-seq Signal Enriched Over Input")
}



