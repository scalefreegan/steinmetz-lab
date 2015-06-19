## + is minus because in 3 Tagseq the strand is flipped

find ../alignments_unique -name "*bam" > allfiles.txt

xargs -a allfiles.txt -n 1 -P 15 -I{} sh -c 'bedtools genomecov -strand - -5 -bg -ibam $1 > "$1.plus.bedGraph"' -- {}
xargs -a allfiles.txt -n 1 -P 15 -I{} sh -c 'bedtools genomecov -strand + -5 -bg -ibam $1 > "$1.minus.bedGraph"' -- {}

# Tool:    bedtools genomecov (aka genomeCoverageBed)
# Version: v2.17.0
# Summary: Compute the coverage of a feature file among a genome.
# 	-strand		Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6). - (STRING): can be + or -
# -5		Calculate coverage of 5" positions (instead of entire interval)
# -bg		Report depth in BedGraph format. For details, see: genome.ucsc.edu/goldenPath/help/bedgraph.html
# -ibam		The input file is in BAM format. Note: BAM _must_ be sorted by position


#bedtools unionbedg  -header -names S96_bio1 S96_bio2 S96_bio3 S96_bio4 -i S96_bio1_plus.bedGraph S96_bio2_plus.bedGraph S96_bio3_plus.bedGraph S96_bio4_plus.bedGraph > plus.bedGraph
#bedtools unionbedg  -header -names S96_bio1 S96_bio2 S96_bio3 S96_bio4 -i S96_bio1_minus.bedGraph S96_bio2_minus.bedGraph S96_bio3_minus.bedGraph S96_bio4_minus.bedGraph > minus.bedGraph

bedtools unionbedg  -header -names $(ls *plus.bedGraph) -i $(ls *plus.bedGraph) > plus.bedGraph 
bedtools unionbedg  -header -names $(ls *minus.bedGraph) -i $(ls *minus.bedGraph) > minus.bedGraph

# Tool:    bedtools unionbedg (aka unionBedGraphs)
# Version: v2.17.0
# Summary: Combines multiple BedGraph files into a single file,
# 	 allowing coverage comparisons between them.
# -header		Print a header line.
# 		(chrom/start/end + names of each file).
# -names		A list of names (one/file) to describe each file in -i.
# 	These names will be printed in the header line.