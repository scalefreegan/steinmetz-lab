# make GRanges Object for Yeast Genome

library( GenomicRanges )

# yeast genome sequence
library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome = BSgenome.Scerevisiae.UCSC.sacCer3