
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

if(!require("pacman")) { 
    install.packages("pacman", repos='http://cran.us.r-project.org') 
}
library(pacman)
p_load(readr)
p_load(dplyr)
p_load(tidyr)
p_load(ggplot2)
p_load(ggthemes)
p_load(gplots)
p_load(grid)
p_load(gridExtra)
p_load(gdata)
p_load(parallel)
p_load(rmarkdown)
p_load(plotly)
p_load(RColorBrewer)
