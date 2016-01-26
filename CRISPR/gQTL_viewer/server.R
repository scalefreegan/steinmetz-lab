# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#
#-------------------------------------------------------------------#
# Shiny interface for ploting/exploring mQTLs
#
#-------------------------------------------------------------------#

.author = "Aaron Brooks"
.copyright = "Copyright 2015"
.credits = "Aaron Brooks"
.license = "WTFPL"
.version = "0.0.1"
.maintainer = "Aaron Brooks"
.email = "aaron.brooks@embl.de"
.status = "Development"
.plot = FALSE


shinyServer(function(input, output, session) {

  options(DT.options = list(pageLength = 10, searching = TRUE))

  gbrowse_link = function(chr,start,end,flanking = c(-2000,2000)) {
    start = start - abs(flanking[1])
    if (start < 0) {
      start = 1
    }
    end = end + abs(flanking[2])
    val = paste("chr",as.roman(chr),":",start,"..",end, sep = "")
    #sprintf('<a href="http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
    #        val)
    sprintf('<a href="http://steinmetzlab.embl.de/mQTL/?loc=%s" target="_blank" class="btn btn-primary">Go to gene</a>',
            val)
  }

  df_query = reactive({
    nenv = apply(gQTL[,4:dim(gQTL)[2]],1,function(i){sum(i<=0.05)})
    wenv = apply(gQTL[,4:dim(gQTL)[2]],1,function(i){paste(colnames(gQTL)[4:dim(gQTL)[2]][which(i<0.05)],collapse=", ")})
    gQTL_df = cbind(gQTL[,1:3],N_env = nenv, env = wenv)
    colnames(gQTL_df) = c("Chr","Start","Stop","N.Envs","Envs")
    nn = sapply(gQTL_df[,"Chr"],function(i){
      paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
    })
    gQTL_df[,"Chr"] = nn
    gQTL_df = gQTL_df[order(gQTL_df[,"N.Envs"],decreasing=T),]
    gQTL_df = cbind(QTL = seq(1,dim(gQTL_df)[1]),gQTL_df)
    return(data.frame(gQTL_df,stringsAsFactors = F))
  })
  
  df_full = reactive({
    s = as.numeric(input$dt1_rows_selected[length(input$dt1_rows_selected)])
    print(s)
    print(values$old_selection)
    #print(s)
    if (length(s)) {
      gQTL_df = df_query()
      # print(gQTL_df[s,])
      # map gQTL locs to genomic ranges with genes, etc
      # check for update to input$cr
      if (length(values$old_selection)==0) {
        # update slider with default QTL coords
        updateSliderInput(session, "cr",
                          value = abs(diff(c(gQTL_df[s,"Start"],gQTL_df[s,"Stop"]))))
      } else if (s != values$old_selection) {
        updateSliderInput(session, "cr",
                          value = abs(diff(c(gQTL_df[s,"Start"],gQTL_df[s,"Stop"]))))
      }
      mid = round(mean(gQTL_df[s,"Start"],gQTL_df[s,"Stop"]))
      nStart =  round(mid - input$cr/2)
      if (nStart < 1) {
        nStart = 1
      }
      nStop = round(mid + input$cr/2)
      QTL_range = IRanges::IRanges(start = nStart, end = nStop)
      QTLgr = GenomicRanges::GRanges(seqnames = gQTL_df[s, "Chr"] , ranges = QTL_range)
      g2QTL = as.data.frame(cdsByOverlaps(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, QTLgr, type = "any", columns = "gene_id"))
      if (dim(g2QTL)[1]>0) {
        # GENE INFOS
        #qtl_df$gene_id = unlist(qtl_df$gene_id)
        gname_t = unlist(gname[unlist(g2QTL$gene_id)])
        gname_t = data.frame(gene_id = names(gname_t), name = gname_t)
        dname_t = unlist(dname[unlist(g2QTL$gene_id)])
        if (is.null(dname_t)) {
          dname_t = data.frame(gene_id = names(dname)[1], alias = unlist(dname)[1])
        } else {
          dname_t = data.frame(gene_id = names(dname_t), alias = dname_t)
        }
        dname_t_long = unlist(dname_long[unlist(g2QTL$gene_id)])
        if (is.null(dname_t)) {
          dname_t = data.frame(gene_id = names(dname_long)[1], alias = unlist(dname_long)[1])
        } else {
          dname_t_long = data.frame(gene_id = names(dname_t_long), desc = dname_t_long)
        }
        g2QTL = merge(g2QTL,gname_t,by="gene_id",sort=F,all.x=T)
        g2QTL = merge(g2QTL,dname_t,by="gene_id",sort=F,all.x=T)
        g2QTL = merge(g2QTL,dname_t_long,by="gene_id",sort=F,all.x=T)
        g2QTL = g2QTL[,c("gene_id","name","seqnames","start","end","strand","alias","desc")]
        colnames(g2QTL) = c("Sys.Name","Name","Chr","Start","End","Strand","Alias","Desc")
        return(g2QTL)
      } else {
        return(data.frame())
      }
    } else {
      return(data.frame())
    }
  })
  
  df_var = reactive({
    # VAR INFOS
    g2QTL = df_full()
    if (dim(g2QTL)[1]>0) {
      var_df = filter(var_info, SNPEFF_TRANSCRIPT_ID%in%unlist(g2QTL$Sys.Name)) %>%
        group_by(.,SNPEFF_TRANSCRIPT_ID) %>%
        do({
          data.frame(
            SNPS = sum(.$id=="snp"),
            INDELS = sum(.$id=="indel"),
            UPSTREAM = sum(.$SNPEFF_EFFECT=="UPSTREAM"),
            DOWNSTREAM = sum(.$SNPEFF_EFFECT=="DOWNSTREAM"),
            INTRONS = sum(.$SNPEFF_EFFECT=="INTRONS"),
            CODING = sum(.$SNPEFF_EFFECT%in%c("UPSTREAM","DOWNSTREAM","INTRONS")==FALSE),
            HIGH = sum(.$SNPEFF_IMPACT=="HIGH"),
            START_LOST = sum(.$SNPEFF_EFFECT=="START_LOST"),
            STOP_GAINED = sum(.$SNPEFF_EFFECT=="STOP_GAINED"),
            STOP_LOST = sum(.$SNPEFF_EFFECT=="STOP_LOST"),
            FRAME_SHIFT = sum(.$SNPEFF_EFFECT=="FRAME_SHIFT"),
            MODERATE = sum(.$SNPEFF_IMPACT=="MODERATE"),
            NON_SYNONYMOUS_CODING = sum(.$SNPEFF_EFFECT=="NON_SYNONYMOUS_CODING"),
            CODON_DELETION = sum(.$SNPEFF_EFFECT=="CODON_DELETION"),
            CODON_INSERTION = sum(.$SNPEFF_EFFECT=="CODON_INSERTION"),
            CODON_CHANGE_PLUS_CODON_DELETION = sum(.$SNPEFF_EFFECT=="CODON_CHANGE_PLUS_CODON_DELETION"),
            CODON_CHANGE_PLUS_CODON_INSERTION = sum(.$SNPEFF_EFFECT=="CODON_CHANGE_PLUS_CODON_INSERTION"),
            LOW = sum(.$SNPEFF_IMPACT=="LOW"),
            SYNONYMOUS_CODING = sum(.$SNPEFF_EFFECT=="SYNONYMOUS_CODING"),
            SYNONYMOUS_STOP = sum(.$SNPEFF_EFFECT=="SYNONYMOUS_STOP"),
            NON_SYNONYMOUS_START = sum(.$SNPEFF_EFFECT=="NON_SYNONYMOUS_START")
          )
        }) %>% ungroup(.)
      return(var_df)
    } else (
      return(data.frame())
    )
  })
  
  df_combine = reactive({
    df_full = df_full()
    if (dim(df_full)[1]>0) {
      tor_df = merge(df_full(),df_var(),by.x="Sys.Name",by.y="SNPEFF_TRANSCRIPT_ID",sort=F,all.x=T)
      tor_df = tor_df[,c("Sys.Name","Name","SNPS", "INDELS", "UPSTREAM",
                         "DOWNSTREAM", "INTRONS", "CODING","HIGH","MODERATE","LOW",
                         "Chr","Start","End","Strand","Alias","Desc")]
      return(tor_df)
    } else {
      return(data.frame())
    }
  })

  output$snptype = renderPlot({
    # reorder
    var_df = df_var()
    if(dim(var_df)[1]>0) {
      high = c("START_LOST","STOP_GAINED","STOP_LOST","FRAME_SHIFT")
      moderate = c("NON_SYNONYMOUS_CODING", "CODON_DELETION","CODON_INSERTION",
                   "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION")
      low = c("SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_START")
      var_df_melt = melt(var_df[,c("SNPEFF_TRANSCRIPT_ID",high,moderate,low)],id.vars="SNPEFF_TRANSCRIPT_ID") %>% filter(.,value>0,variable%in%c(high,moderate,low))
      var_df_melt$impact = sapply(var_df_melt$variable,function(i){
        if (i %in% high) {
          return("High Impact")
        } else if (i %in% moderate) {
          return("Moderate Impact")
        } else {
          return ("Low Impact")
        }
      })
      var_df_melt2 = do.call(rbind,lapply(seq(1,dim(var_df_melt)[1]),function(i){
        do.call(rbind,lapply(seq(1,as.numeric(var_df_melt[i,"value"])),function(j){
          var_df_melt[i,]
        }))
      }))
      var_df_melt2$impact = factor(var_df_melt2$impact, levels = c("High Impact","Moderate Impact","Low Impact"))
      var_df_melt2$SNPEFF_TRANSCRIPT_ID = factor(var_df_melt2$SNPEFF_TRANSCRIPT_ID, levels = sort(unique(var_df_melt2$SNPEFF_TRANSCRIPT_ID)))
      ggplot(var_df_melt2, aes(x = factor(SNPEFF_TRANSCRIPT_ID), fill=variable)) + geom_bar(width=.8) +
        coord_flip() + facet_wrap(~ impact) + scale_x_discrete(limits=rev(levels(var_df_melt2$SNPEFF_TRANSCRIPT_ID))) +
        xlab("Gene") + ylab("# SNPs/Indels") + theme(legend.position="bottom")
    } else {
      return(NULL)
    }
  })

  # DATA TABLES
  output$dt1 = DT::renderDataTable(
    df_query(), server = TRUE, selection = "multiple",
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)
  
  output$dt2 = DT::renderDataTable(
    df_combine(), server = TRUE, selection = "none",
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)

  # REACTIVE VALUES
  values = reactiveValues(
    old_selection = 0,
    old_interval = NULL,
    link = NULL
  )

  # MONITOR OLD SELECTION
  session$onFlush(once=FALSE, function(){
    isolate({ values$old_selection <- as.numeric(input$dt1_rows_selected[length(input$dt1_rows_selected)]) })
  })

  output$image <- renderText({
   paste0('<td align="middle"><img src="http://scalefreegan.github.io/steinmetz-lab/CRISPR/data/ge_vals.png" ',
           'valign="middle" style="width: 40%;max-height: 100%"></td>')
    })
  
  output$link = renderText({
    s = as.numeric(input$dt1_rows_selected[length(input$dt1_rows_selected)])
    d = df_query()
    print(s)
    if (length(s)) {
      chr = d[s, "Chr"]
      start = d[s, "Start"]
      end = d[s, "Stop"]
      mid = round(mean(c(start,end)))
      #a = curlEscape(paste0('[{ "seq_id":"chrII", "start": ', start, ', "end": ',end,', "name": "QTL"}]'))
      a = paste0(chr,"%3A",start,"..",end)
      val = paste(chr,"%3A",mid-input$cr,"..",mid+input$cr, sep = "")
      val = paste('http://steinmetzlab.embl.de/mQTL/?loc=',val,"&highlight=", a, sep="")
      # qtl segment to add
      print(val)
    } else {
      val = "http://steinmetzlab.embl.de/mQTL/"
    }
    paste("<div style='width: 100%; height: 600px'><iframe style='border: 1px solid black' src='", val ,"'width='100%' height='100%'></iframe></div>",sep="")
  })

})
