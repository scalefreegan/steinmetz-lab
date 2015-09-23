# various functions used to desgin guide RNA
# Author: czhu
###############################################################################

design_basic_gRNA = function(seqs) {
    ## given a sequence, find the PAM site (with PAM.MOTIF) and retrieve the sequence as gRNA 
    ## if the given sequence is not long enough (at the boundary), it ignores these PAM sites 
    require(Biostrings)
    GRNA.LENGTH=20
    ## NGG
    PAM.MOTIF ="GG"
    mysubseq = function(x,ir) {
        ## given ir, IRanges 
        ## x, DNAStringSet
        ## retrieve sequences 
        if(length(ir)>0 & (sum((start(ir) -GRNA.LENGTH)>0)>0)) {
            sel = which((start(ir) - 1 - GRNA.LENGTH)>0)
            
            rv = data.frame(
                pStart =  start(ir)[sel] - 1 - GRNA.LENGTH,
                pEnd   =  end(ir)[sel] - nchar(PAM.MOTIF) - 1, ## minus the N site
                stringsAsFactors=FALSE)            
            rv$seq = as.character(subseq(rep(x, length(sel)), 
                    start=rv$pStart, end=rv$pEnd))
            rv
                
        } else {
            rv = NULL
        }
        rv
    } 
#    if(is.null(names(seqs))) {
#        names(seqs) = paste0("target",1:length(seqs))
#    }
    #seqs = fsa
    pam_pos = vmatchPattern(PAM.MOTIF,seqs,fixed=FALSE)
    ## merge the polymers
    pam_pos = lapply(pam_pos, function(x){
            rv = reduce(x)
            sel = width(rv)>nchar(PAM.MOTIF)
            start(rv[sel]) = start(rv[sel]) + 1 
            end(rv[sel]) = start(rv[sel]) + 1
            rv
        })
    gRNA.seq = lapply(1:length(pam_pos),  function(i){
            mysubseq(seqs[i], pam_pos[[i]])
        }
    )
    names(gRNA.seq) = names(seqs)
    rv = data.frame(
        index=rep(1:length(seqs), sapply(gRNA.seq,nrow)),
        gSeq=unlist(sapply(gRNA.seq, "[","seq")),
        pStart=unlist(sapply(gRNA.seq, "[","pStart")),
        pEnd=unlist(sapply(gRNA.seq, "[","pEnd")),
        stringsAsFactors=FALSE)    
    rv
}


assess_bowtie2 = function(candi) {    
    ################################################################################
    ## method 1 use bowtie
    ## limitation, only consider cases with 1 MM
    ## and bowtie2 might not be sensitive enough
    N.MISMATCH = 1
    BW.INDEX= "/g/steinmetz/czhu/genome/S288c/indexes/bowtie2/s288c_20110203"
    require(GenomicAlignments)
    
    fafile = tempfile()
    alnfile = tempfile()
    
    writeXStringSet(candi, filepath=fafile)
    cat("Writing fasta to", fafile, "\n")
    cmd = paste("bowtie2 -x", BW.INDEX, 
        "-f", fafile ,"--very-sensitive -N 1 --all -L 23 | samtools view -b - >", alnfile) 
    cat("Running",cmd,"\n")
    system(cmd)
    myparam = ScanBamParam(tag=c("XM", "MD"),what=c("qname"))
    ga = readGAlignments(alnfile, param=myparam)
    
    unlink(c(fafile,alnfile))
    ## only reads with 1 mismatch, no gaps
    ga= ga[cigar(ga)==paste0(qwidth(ga),"M")]
    ga= ga[values(ga)$XM<2]
    
    if(length(ga)==0){
        return(NULL)
    } else {
        spl = strsplit(values(ga)$MD, "\\D")
        
        mmpos = rep(0,length(ga))
        hasMM = sapply(spl, length)>1
        mmpos[hasMM] = as.integer(sapply(spl[hasMM],"[",1)) + 1
        
        da = data.frame(
            subject = values(ga)$qname,
            chr= seqnames(ga),
            start = start(ga),
            end=end(ga),
            strand=strand(ga),
            width = width(ga),
            nMM = values(ga)$XM,
            MMpos = mmpos,
            stringsAsFactors=FALSE)
        return(da)
    }
    
}

## ##############################################################################
## method 2 use R's pattern recogniction. This is much slower but should be more accurate 
## this is not implemented
#lapply(names(genome), function(chr){
#        plus_hits = matchPDict(candi,genome[[chr]],max.mismatch=N.MISMATCH)
#    })
offtarget_analysis = function(candi){
    
    if(class(candi) != "DNAStringSet") {
        candi = DNAStringSet(candi)
    }
    ## if the method is implemented, this is the place to choose which method
    rp = assess_bowtie2(candi)
    if(is.null(rp)) {
        stop("nothing mapped\n")
    } else {
        fa = factor(rp$subject,levels=unique(rp$subject))
        finalrp =by(rp, fa, function(x){
                rv = data.frame(
                    nUni=sum(x$nMM==0), nMM=sum(x$nMM==1),
                    unisite="",offtargets="",mutPos ="", stringsAsFactors=FALSE)
                if(rv$nUni>1) {
                    sel = which(x$nMM ==0)
                    rv$unisite = paste(paste(x$chr,x$start,x$end,x$strand,sep="_")[sel], collapse=",")
                }
                if(rv$nMM>0){
                    sel = which(x$nMM ==1)
                    rv$offtargets = paste(paste(x$chr,x$start,x$end,x$strand,sep="_")[sel], collapse=",")
                    rv$mutPos = paste(x$MMpos[sel], collapse=",")
                } 
                rv
            })
        finalrp=do.call(rbind,finalrp)
        return(finalrp)        
    }
}
## for backward compatibility
offsite_analysis = offtarget_analysis 

## extract sequence from a given genome
getseq = function(dat,genome){
    ## dat, defines the region with chr, strand, start and end
    ## genome is the genome sequence as DNAStringSet object
    require(Biostrings)
    rv = DNAStringSet(lapply(1:nrow(dat), function(i){
                subseq(genome[[dat$chr[i]]], start=dat$start[i], end=dat$end[i]) 
            }))
    isMinus = which(dat$strand == "-")
    if(length(isMinus)>0) {
        rv[isMinus] = reverseComplement(rv[isMinus])       
    }
    rv
}


gcContent <-function(x) {
    if(is.character(x)) 
        x = DNAStringSet(x)
    alf <- alphabetFrequency(x, as.prob=TRUE)
    rowSums(alf[,c("G", "C"),drop=FALSE])
}



