
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
library(pacman)
p_load(gplots)
p_load(parallel)

argv = commandArgs (TRUE)
mirRelativePos = read.table (argv[1], F, sep="\t", stringsAsFactors=F)

hetName = argv[2]
mutName = argv[3]
#hetName = gsub ("\\.", "_", argv[2]) 
#mutName = gsub ("\\.", "_", argv[3]) 
numOfCore = argv[4]
outDir = argv[5]
mirRelativePosSplitted = split (mirRelativePos, mirRelativePos$V1)
mclapply (mirRelativePosSplitted, draw_microRNA_balloon, mc.cores=numOfCore, hetName = hetName, mutName = mutName, outDir = outDir)
