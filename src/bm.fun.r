#{{{ load required libraries, define common variables
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
source("circlePlot.R")
options(stringsAsFactors = F)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
source('~/projects/maize.expression/src/me.fun.r')
source('~/projects/genomes/src/ge.fun.r')


dirp = '~/projects/biomap'
dird = file.path(dirp, 'data')
#}}}




