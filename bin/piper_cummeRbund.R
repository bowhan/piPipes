
# piper, a pipeline collection for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) <2014>  <Bo Han, Wei Wang, Phillip Zamore, Zhiping Weng>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


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
