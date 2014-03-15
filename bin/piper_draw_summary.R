
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

pkgTest ("multicore")

argv  = commandArgs (TRUE)
summaryTableFile = argv[1]
outPdfPrefix = argv[2]
numOfCore = argv[3]
normScale = as.real (argv[4])
summaryTable = read.table (summaryTableFile, F)
summaryTableSplited = split (summaryTable, summaryTable$V1)
mclapply (summaryTableSplited, draw_summary, mc.cores=numOfCore, pdfPrefix=outPdfPrefix, normScale=normScale)
