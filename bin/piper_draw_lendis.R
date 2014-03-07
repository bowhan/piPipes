# small RNA pipeline from pipipe: https://github.com/bowhan/pipipe.git
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bowhan@me.com)
# the Zamore lab and the Weng lab
# University of Massachusetts Medical School

source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/pipipe.R",sep=""))

pkgTest ("ggplot2")
pkgTest ("scales")

argv = commandArgs (TRUE)
lendis = read.table (argv[1],FALSE)
main = argv[2]
pdf (paste (main, ".lendis.pdf", sep=''))
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 50), collapse = "\n") 
minRow = min(lendis[lendis[,2]!=0,1])
maxRow = max(lendis[lendis[,2]!=0,1])
lendis = lendis[seq(minRow, maxRow),]
ggplot (lendis, aes (V1,V2)) + 
	geom_bar (stat="identity") + 
	scale_x_discrete (breaks=c(minRow:maxRow)) + 
	coord_cartesian(xlim = c(minRow, maxRow)) + 
	scale_y_continuous(labels = comma, breaks=seq(0,max(lendis$V2),roundUp(max(lendis$V2)/10))) +
	labs(title=paste("Length Distribution",main,sep="\n")) + 
	xlab("Length, nt") + 
	ylab("Reads")
noprint = dev.off ()