# ARM-seq fetch

## (1) Get Data

Fetch ARM-seq data from NCBI using SRA-toolkit. Process using command line tools

`cd /g/steinmetz/brooks/data/iesy/armseq/untreated/`

Download untreated fastq files

```
fastq-dump SRR1874045
fastq-dump SRR1874048
fastq-dump SRR1875056

```

`cd /g/steinmetz/brooks/data/iesy/armseq/alkb_treated/`

Download AlkB-treated fastq files

```
fastq-dump SRR1874029
fastq-dump SRR1874032
fastq-dump SRR1874034
```

## (2) Get Genome

Download BY4741 genome from SGD and make Bowtie 2 index

`cd /g/steinmetz/brooks/data/iesy/armseq/genome/`

`wget http://downloads.yeastgenome.org/sequence/strains/BY4741/BY4741_Stanford_2014_JRIS00000000/BY4741_Stanford_2014_JRIS00000000.fsa.gz`

`cd ./bowtie2_index`

`bowtie2-build ../BY4741_Stanford_2014_JRIS00000000.fsa BY4741`

## (3) Align to Genome

`cd /g/steinmetz/brooks/data/iesy/armseq/aligments`

`bowtie2 -x ../genome/bowtie2_index/BY4741 -U ../untreated/SRR1874045.fastq -S untreated-SRR1874045.sam`

bowtie2 -x /g/steinmetz/project/GenPhen/data/genome/s288c-R64_ERCC -U ../untreated/SRR1874045.fastq -S untreated-SRR1874045.sam --sensitive
