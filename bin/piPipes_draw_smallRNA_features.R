
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
pkgTest ("ggthemes")
pkgTest ("scales")
pkgTest ("reshape")
pkgTest ("gridExtra")

argv  = commandArgs (TRUE)
# reading options
main = argv[1]
total_lendis = argv[2]
siRNA_lendis = argv[3]
piRNA_lendis = argv[4]
ext_len = as.numeric (argv[5])
total_5end_per = argv[6]
total_ping_pong = argv[7]
siRNA_5end_per = argv[8]
siRNA_ping_pong = argv[9]
piRNA_5end_per = argv[10]
piRNA_ping_pong = argv[11]

main=basename (main)
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 80), collapse = "\n") 
pdf (paste (argv[1], ".pdf", sep=''), onefile=TRUE, width=8.5, height=11, title=main)

total_lendis_gg = draw_smRNA_lendis (total_lendis, "smRNA lendis")
siRNA_lendis_gg = draw_smRNA_lendis (siRNA_lendis, "siRNA lendis")
piRNA_lendis_gg = draw_smRNA_lendis (piRNA_lendis, "piRNA lendis")

total_5end_per_gg = draw_smRNA_percentage (total_5end_per, ext_len, "smRNA 5' end ext")
total_ping_pong_gg = draw_ping_pong (total_ping_pong, "smRNA")

siRNA_5end_per_gg = draw_smRNA_percentage (siRNA_5end_per, ext_len, "siRNA 5' end ext")
siRNA_ping_pong_gg = draw_ping_pong (siRNA_ping_pong, "siRNA")

piRNA_5end_per_gg = draw_smRNA_percentage (piRNA_5end_per, ext_len, "piRNA 5' end ext")
piRNA_ping_pong_gg = draw_ping_pong (piRNA_ping_pong, "piRNA")

grid.arrange(
		total_ping_pong_gg,
		siRNA_ping_pong_gg,
		piRNA_ping_pong_gg,
		total_5end_per_gg, 
		siRNA_5end_per_gg, 
		piRNA_5end_per_gg, 
		total_lendis_gg, 
		siRNA_lendis_gg, 
		piRNA_lendis_gg, 
		ncol=3,
		as.table=TRUE, 
		main = textGrob(main, vjust = 1, gp = gpar(fontface = "bold", cex = 1))
)
gc = dev.off()
