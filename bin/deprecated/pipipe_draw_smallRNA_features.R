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
total_lendis = argv[2]
siRNA_lendis = argv[3]
piRNA_lendis = argv[4]
ext_len = as.numeric (argv[5])
total_5end_per = argv[6]
total_3end_per = argv[7]
siRNA_5end_per = argv[8]
siRNA_3end_per = argv[9]
piRNA_5end_per = argv[10]
piRNA_3end_per = argv[11]

main=basename (main)
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 80), collapse = "\n") 
pdf (paste (argv[1], ".pdf", sep=''), onefile=TRUE, width=8.5, height=11, title=main)

total_lendis_gg = draw_smRNA_lendis (total_lendis, "smRNA lendis")
siRNA_lendis_gg = draw_smRNA_lendis (siRNA_lendis, "siRNA lendis")
piRNA_lendis_gg = draw_smRNA_lendis (piRNA_lendis, "piRNA lendis")

total_5end_per_gg = draw_smRNA_percentage (total_5end_per, ext_len, "smRNA 5' end ext")
total_3end_per_gg = draw_smRNA_percentage (total_3end_per, ext_len, "smRNA 3' end ext")

siRNA_5end_per_gg = draw_smRNA_percentage (siRNA_5end_per, ext_len, "siRNA 5' end ext")
siRNA_3end_per_gg = draw_smRNA_percentage (siRNA_3end_per, ext_len, "siRNA 3' end ext")

piRNA_5end_per_gg = draw_smRNA_percentage (piRNA_5end_per, ext_len, "piRNA 5' end ext")
piRNA_3end_per_gg = draw_smRNA_percentage (piRNA_3end_per, ext_len, "piRNA 3' end ext")

grid.arrange(
		total_3end_per_gg,
		siRNA_3end_per_gg,
		piRNA_3end_per_gg,
		total_5end_per_gg, 
		siRNA_5end_per_gg, 
		piRNA_5end_per_gg, 
		total_lendis_gg, 
		siRNA_lendis_gg, 
		piRNA_lendis_gg, 
		ncol=3,
		as.table=TRUE, 
		main = textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1))
)
gc = dev.off()
