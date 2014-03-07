# pipipe: https://github.com/bowhan/pipipe.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Zlab
# University of Massachusetts Medical School

source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/pipipe.R",sep=""))

pkgTest ("multicore")

argv  = commandArgs (TRUE)
summaryTableFile = argv[1]
outPdfPrefix = argv[2]
numOfCore = argv[3]
normScale = as.real (argv[4])
summaryTable = read.table (summaryTableFile, F)
summaryTableSplited = split (summaryTable, summaryTable$V1)
mclapply (summaryTableSplited, draw_summary, mc.cores=numOfCore, pdfPrefix=outPdfPrefix, normScale=normScale)
