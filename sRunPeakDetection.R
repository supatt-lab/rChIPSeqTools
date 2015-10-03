source("sPeakDetection.R")

getChIPSeqPeaks("example.data/IRF2BP2_ChIP.sorted.chrX.bam","example.data/IRF2BP2_IgG_control.sorted.chrX.bam",fragment.size=200,"chrX","example.data/IRF2BP2.chrX.peaks.txt")


