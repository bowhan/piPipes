# small RNA pipeline from pipipe: https://github.com/bowhan/pipipe.git
# pipipe: https://github.com/bowhan/pipipe.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Zlab
# University of Massachusetts Medical School


source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/pipipe.R",sep=""))

pkgTest("gplots")
pkgTest("multicore")

argv = commandArgs (TRUE)
mirRelativePos = read.table (argv[1], F, sep="\t", stringsAsFactors=F)
hetName = argv[2]
mutName = argv[3]
numOfCore = argv[4]
outDir = argv[5]
mirRelativePosSplitted = split (mirRelativePos, mirRelativePos$V1)
mclapply (mirRelativePosSplitted, draw_microRNA_balloon, mc.cores=numOfCore, hetName = hetName, mutName = mutName, outDir = outDir)