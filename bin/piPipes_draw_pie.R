
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

source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/piPipes.R",sep=""))

pkgTest ("RColorBrewer")

argv = commandArgs (TRUE)
pdf (paste (argv[1], ".pdf", sep=''), width=10, height=10)
par (mar=c(5, 5, 5, 5))
counts = read.table (argv[2], FALSE, sep="\t")
main = basename (argv[1])
main = gsub ("\\."," ",main)
main = paste (strwrap(main, width = 50), collapse = "\n")
percent = round (counts[,2]/sum(counts[,2])*100, digits=1)
lbls = paste (counts[,1], percent)
lbls = paste (lbls,"%",sep="")
lbls = gsub ('Five', '5', lbls)
lbls = gsub ('Three', '3', lbls)
lbls = gsub ('Prime', '\'', lbls)
lbls = gsub ('Prime', '\'', lbls)
pie(counts [,2], labels = lbls, main=main, col= brewer.pal(nrow(counts),"Set1"), clockwise=TRUE)
noprint = dev.off ()
