
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

pkgTest ("ggplot2")
#pkgTest ("extrafont")
pkgTest ("gplots")
pkgTest ("multicore")
pkgTest ("scales")
pkgTest ("reshape")
pkgTest ("gridExtra")
pkgTest ("RColorBrewer")
pkgTest ("gdata")


bioConductorTest ("cummeRbund")
#bioConductorTest ("edgeR")
