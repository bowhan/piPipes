# pipipe: https://github.com/bowhan/pipipe.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Zlab
# University of Massachusetts Medical School

source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/pipipe.R",sep=""))

pkgTest ("ggplot2")
pkgTest ("extrafont")
pkgTest ("RColorBrewer")
argv  = commandArgs (TRUE)
pdf (paste(argv[2], ".pdf", sep=''))
draw_agg (argv[1], argv[3])
invisible (dev.off )

