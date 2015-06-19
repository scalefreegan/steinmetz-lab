#infile="../rawdata/C6LRUACXX_SxY_Tagseq_P1A_15s009414-1-1_Tekkedil_lane510Abiorep1_sequence.txt.gz"

#bowtie2 -p 1 --sensitive --met-file "stats/`basename $infile`.stat" -x /g/steinmetz/project/GenPhen/data/genome/s288c-R64_ERCC -U $infile | samtools view -uS - | samtools sort -@ 1 -f - "`basename $infile`.bam"

# bowtie2 options
# Bowtie 2 outputs alignments in SAM format
# -p/--threads <int> number of alignment threads to launch (1)
# --sensitive  -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
# --met-file <path>  send metrics to file at <path> (off)
# -X <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
#               NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
# -U  <r>  Files with unpaired reads.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#
# samtools options
# view        SAM<->BAM<->CRAM conversion
# -u       uncompressed BAM output (implies -b)
# -S       ignored (input format is auto-detected)
# sort        sort alignment file
# -@ INT     Set number of sorting and compression threads [1]
# -f         Use <out.prefix> as full final filename rather than prefix
# -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)

find ../rawdata/ -name "*SxY*txt.gz" > allfiles.txt
## bowtie way
xargs -a allfiles.txt -n 1 -P 5 -I{} sh -c 'bowtie2 -p 6 --very-sensitive -x /g/steinmetz/project/GenPhen/data/genome/s288c-R64_ERCC --trim3 1 -U $1 | samtools view -uS - | samtools sort -@ 1 - "`basename $1`"' -- {}

# --arg-file=file, -a file Read items from file instead of standard input.  If you use this option, stdin remains unchanged when commands are run.	Otherwise, stdin is redirected from /dev/null.
# --max-args=max-args, -n max-args. Use  at  most  max-args  arguments per command line.  Fewer than max-args arguments will be used if the size (see the -s  option)is  exceeded, unless the -x option is given, in which case xargs will exit.
# --max-procs=max-procs, -P max-procs. Run up to max-procs processes at a time; the default is  1.   If max-procs	 is 0, xargs will run as many processes as possible at a time.  Use the -n option with -P; otherwise chances  are  that only one exec will be done.
# -I replace-str. Replace occurrences of replace-str in the initial-arguments with names read from standard input.  Also, unquoted  blanks  do  not terminate	 input	items;	instead	 the  separator is the newline character.  Implies -x and -L 1.
## bwa way
xargs -a allfiles.txt -n 1 -P 1 -I{} sh -c '../bwa.sh $1 `basename $1`' -- {}

## filter for unique alignment
xargs -a allfiles.txt -n 1 -P 20 -I{} sh -c 'samtools view -q 10 -b $1 >"`basename $1`"' -- 
# -q INT   only include reads with mapping quality >= INT [0]
# -b       output BAM

