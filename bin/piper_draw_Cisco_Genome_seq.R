
# piper, a pipeline collection for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) 2014  Bo Han, Wei Wang, Phillip Zamore, Zhiping Weng
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


source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/piper.R",sep=""))

pkgTest ("RCircos")
argv  = commandArgs (TRUE)
cyto.info = read.table (argv[1],F)
piRNA_cluster = read.table(gzfile(argv[2]), header=F)
TEMP.out = read.table (argv[3],T)
retroSeq.out = read.table (argv[4],T)
BrD.out = read.table (argv[5],T)
pdf (argv[6], heigh=10, width=10)
tracks.inside  = 6
tracks.outside = 0
RCircos.Set.Core.Components(cyto.info, NULL, tracks.inside, tracks.outside);
RCircos.Set.Plot.Area()
par(mai=c(.1,.1,.1,.1))
plot.window(c(-2,2), c(-2, 2))
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot (piRNA_cluster, 1 ,"in")
RCircos.Gene.Name.Plot (piRNA_cluster, 4, 2 ,"in")
RCircos.Tile.Plot (TEMP.out, 4 ,"in")
RCircos.Tile.Plot (retroSeq.out, 5 ,"in")
RCircos.Link.Plot (BrD.out, 6, "in")
gc = dev.off ()