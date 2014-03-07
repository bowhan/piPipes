# small RNA pipeline from piper: https://github.com/bowhan/piper.git
# piper: https://github.com/bowhan/piper.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Zlab
# University of Massachusetts Medical School


source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/piper.R",sep=""))

pkgTest ("ggplot2")
pkgTest ("scales")
pkgTest ("reshape")
pkgTest ("gridExtra")

argv  = commandArgs (TRUE)
# reading options
main = argv[1]
ext_len = as.numeric (argv[2])
deg_5end_per = argv[3]

main=basename (main)
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 80), collapse = "\n") 
pdf (paste (argv[1], ".pdf", sep=''), onefile=TRUE, width=8.5, height=11, title=main)

deg_5end_per_gg = draw_smRNA_percentage (deg_5end_per, ext_len, "degradome 5' end ext")

grid.arrange(
		deg_5end_per_gg,
		ncol=1,
		as.table=TRUE, 
		main = textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1))
)
gc = dev.off()
