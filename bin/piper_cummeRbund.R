# piper: https://github.com/bowhan/piper.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Zlab
# University of Massachusetts Medical School

source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/piper.R",sep=""))

bioConductorTest ("cummeRbund")

topN = 50
argv = commandargv (TRUE)

cuff = readCufflinks (argv[1], package = "cummeRbund")
pdf (paste (argv[2],'.genes.csDensity.pdf',sep=''))
csDensity (genes(cuff))
invisible(dev.off())

pdf (paste (argv[2],'.genes.csBoxplot.pdf',sep=''))
csBoxplot (genes(cuff))
invisible(dev.off())

names = samples (cuff)$sample_name

pdf(paste (argv[2],"_",names[1],"_vs_",names[2],".genes.csScatter.pdf",sep=""))
csScatter(genes(cuff),names[1],names[2],smooth=T)
invisible(dev.off())

pdf(paste (argv[2],"_",names[1],"_vs_",names[2],".genes.csVolcano.pdf",sep=""))
csVolcano(genes(cuff),names[1],names[2])
invisible(dev.off())

gene.diff = diffData(genes(cuff))
gene.diff.top = gene.diff[order(gene.diff$q_value),][1:topN,]
myGeneIds = gene.diff.top$gene_id
myGenes = getGenes(cuff, myGeneIds)

pdf (paste (argv[2],'.genes.csHeatMap_top', topN,'.pdf',sep=''))
csHeatmap(myGenes, cluster="both")
invisible(dev.off())

pdf (paste (argv[2],'.genes.csExpressionBarplot_top', topN, '.pdf',sep=''))
expressionBarplot(myGenes)
invisible(dev.off())

pdf (argv[2],'.isoforms.csDensity.pdf',sep=''))
csDensity (isoforms(cuff))
invisible(dev.off())

pdf (paste (argv[2],'.isoforms.csBoxplot.pdf',sep=''))
csBoxplot (isoforms(cuff))
invisible(dev.off())
