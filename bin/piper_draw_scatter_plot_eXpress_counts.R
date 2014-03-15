
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
pkgTest ("gdata")

argv = commandArgs (TRUE)
sample1 = read.table (argv[1], T)
sample2 = read.table (argv[2], T)
name1 = argv[3]
name2 = argv[4]

main = basename (argv[5])
main = gsub ("\\."," ",main)
main = paste (strwrap(main, width = 80), collapse = "\n")

pdf (paste (argv[5],".pdf", sep=""), onefile=TRUE, width=10, height=10, title=main )

sample = merge (sample1, sample2, by="target_id")
color = ifelse (startsWith(sample[,1],"NM_"), "black", ifelse (startsWith(sample[,1], "NR_"), "darkgrey", "red"))
size = ifelse (startsWith(sample[,1],"NM_"), 1, ifelse (startsWith(sample[,1], "NR_"), 1, 1.5))

lim = roundUp (10*log2(max (sample$eff_counts.x, sample$eff_counts.y)))/10

g = plot ( log2(sample$eff_counts.x), log2(sample$eff_counts.y), 
		xlim=c(-lim, lim), 
		ylim=c(-lim, lim), 
		xlab=paste(name1,"normalized count,log2"), 
		ylab=paste(name2,"normalized count,log2"), 
		main = main, 
		pch=21, col="white", bg=color, cex=size, xaxt='n', yaxt='n', frame=F ) + 
abline (0,1, lty=2) + 
axis (1, tck=0.01, lwd=1, cex.axis=1) + 
axis (2, tck=0.01, lwd=1, cex.axis=1)
legend ("topleft", c("NM_", "NR_", "transposon & cluster"), col=c("black","darkgrey","red"), pch=16, cex=1.25, box.lwd=0, box.col="transparent", bg="transparent")
g = dev.off ()
