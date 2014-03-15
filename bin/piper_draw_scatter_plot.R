
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
argv = commandArgs (TRUE)
sample1 = read.table (argv[1])
sample2 = read.table (argv[2])
name1 = argv[3]
name2 = argv[4]

main = basename (argv[5])
main = gsub ("\\."," ",main)
main = paste (strwrap(main, width = 80), collapse = "\n")

pdf (paste (argv[5],".pdf", sep=""), onefile=TRUE, width=10, height=10, title=main )

par (mfrow=c(2,2), mai=c(0.8, 0.8, 0.4, 0.1))

if (ncol(sample1)==4) { # with grouping information 
        colnames (sample1) = c("name", "group", "S1", "AS1")
        colnames (sample2) = c("name", "group", "S2", "AS2")
        sample = merge (sample1, sample2, by="name")
        color = ifelse (
                sample1$group == 1, "black", ifelse (
                        sample1$group== 2, "green", ifelse (
                                sample1$group == 3, "red", "goldenrod3"
        )))
} else {
        color = "black"
        colnames (sample1) = c("name", "S1", "AS1")
        colnames (sample2) = c("name", "S2", "AS2")
        sample = merge (sample1, sample2, by="name")
}
#lim1 = roundUp (10*log2(min (sample$S1, sample$S2, sample$AS1, sample$AS2))-1)/10
lim1 = 0
lim2 = roundUp (10*log2(max (sample$S1, sample$S2, sample$AS1, sample$AS2)))/10
g = plot (log2(sample$S1), log2(sample$AS1),   xlim=c(lim1, lim2), ylim=c(lim1, lim2), xlab=paste(name1,"sense count,log2"), ylab=paste(name1,"antisense count,log2"), pch=21, col="white", bg=color, cex=1.5, xaxt='n', yaxt='n', frame=F) + abline (0,1, lty=2) + axis (1, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1) + axis (2, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1)
g = plot (log2(sample$S2), log2(sample$AS2),   xlim=c(lim1, lim2), ylim=c(lim1, lim2), xlab=paste(name2,"sense count,log2"), ylab=paste(name2,"antisense count,log2"), pch=21, col="white", bg=color, cex=1.5, xaxt='n', yaxt='n', frame=F) + abline (0,1, lty=2) + axis (1, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1) + axis (2, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1)
g = plot (log2(sample$S1), log2(sample$S2),    xlim=c(lim1, lim2), ylim=c(lim1, lim2), xlab=paste(name1,"sense count,log2"), ylab=paste(name2,"sense count,log2"), pch=21, col="white", bg=color, cex=1.5, xaxt='n', yaxt='n', frame=F) + abline (0,1, lty=2) + axis (1, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1) + axis (2, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1)
g = plot (log2(sample$AS1),log2(sample$AS2),   xlim=c(lim1, lim2), ylim=c(lim1, lim2), xlab=paste(name1,"antisense count,log2"), ylab=paste(name2,"antisense count,log2"), pch=21, col="white", bg=color, cex=1.5, xaxt='n', yaxt='n', frame=F) + abline (0,1, lty=2) + axis (1, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1) + axis (2, tck=0.01, lwd=1, at=as.integer(seq(lim1,lim2,length.out=5)), cex.axis=1)
g = dev.off ()
