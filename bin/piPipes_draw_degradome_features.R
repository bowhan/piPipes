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

pkgTest ("ggplot2")
pkgTest ("scales")
pkgTest ("reshape")
pkgTest ("gridExtra")

argv  = commandArgs (TRUE)
# reading options
main = argv[1]
ext_len = as.numeric (argv[2])
deg_5end_per = argv[3]

main=basename (main)
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 80), collapse = "\n") 
pdf (paste (argv[1], ".pdf", sep=''), onefile=TRUE, width=8.5, height=11, title=main)

deg_5end_per_gg = draw_smRNA_percentage (deg_5end_per, ext_len, "degradome 5' end up/down stream")

if (length(argv)==4) {
	deg_smRNA_pp = total_ping_pong_gg = draw_ping_pong (argv[4], "smRNA-degradome")
	grid.arrange(
		deg_5end_per_gg,
		deg_smRNA_pp,
		ncol=1,
		as.table=TRUE, 
		main = textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1))
	)
} else {
	grid.arrange(
		deg_5end_per_gg,
		ncol=1,
		as.table=TRUE, 
		main = textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1))
	)
}
gc = dev.off()
