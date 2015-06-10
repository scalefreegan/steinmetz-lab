# make GRanges Object for Yeast Genome

library( GenomicRanges )

# yeast genome sequence
library("BSgenome.Scerevisiae.UCSC.sacCer3")
genome = BSgenome.Scerevisiae.UCSC.sacCer3

genome_seqs = seqlengths(genome)

yeast_gr = unlist( GRangesList( lapply(1:length(genome_seqs),function(i){
	GRanges(
	seqnames = names( genome_seqs )[i],
	ranges = IRanges( start=1, end = genome_seqs[i]),
	strand = strand( "*" ),
	GC = signif( sum(alphabetFrequency(getSeq(genome,names( genome_seqs )[i]))[c("G","C")]) / sum(alphabetFrequency(getSeq(genome,names( genome_seqs )[i]))), 2)
	) } ) ) )
 
 save( yeast_gr, file="")