library(ggplot2)

df <- read.table('/data/RDCO/kevbrick/GL_Sorting_Paper/nxfOut1/peakCalls/Ab8580_H3K4me3_ChIPSeq.allPeakCounts.forR.tab',header=TRUE)

ggplot(df,aes(x=Nuclei,y=peaks)) + geom_point(shape=23,size=.7) + geom_smooth() +
xlab('Nuclei (#)') + ylab('Peaks')

ggsave('test.png')
