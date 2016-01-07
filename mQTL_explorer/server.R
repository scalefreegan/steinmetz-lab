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

  df_full = reactive({

    # select chromosomes
    chrs = unique(data[[input$m]]$qtl[data[[input$m]]$qtl[,type]>=summary(data[[input$m]]$permout[,type],input$co/100)[1],"chr"])
    chrs = levels(chrs)[chrs]

    lodcolumn = if(type=="mlod"){ 2 } else { 1 }
    qtl_intervals = list()
    if (length(chrs)>0) {
      for (i in chrs) {
        qtl_intervals[[i]] = try(mrk[rownames(bayesint(data[[input$m]]$qtl, chr = str_pad(i, 2, pad = "0"), prob=input$bci/100, lodcolumn=lodcolumn))],silent = T)
        if (class(qtl_intervals[[i]])=="try-error") {
          qtl_intervals[[i]] = NULL
        } else {
          nn = sapply(as.character(seqnames(qtl_intervals[[i]])),function(i){
            paste(substr(i,1,3),as.roman(substr(i,4,5)),sep="")
          })
          qtl_intervals[[i]] = renameSeqlevels(qtl_intervals[[i]],nn)
          qtl_intervals[[i]] = keepSeqlevels(qtl_intervals[[i]],unique(nn))
          qtl_intervals[[i]] = range(qtl_intervals[[i]])
          qtl_intervals[[i]] = as.data.frame(cdsByOverlaps(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,qtl_intervals[[i]], type = "any", columns = "gene_id"))
        }
      }
    }
    qtl_df = do.call(rbind,qtl_intervals)
    if (length(qtl_df) != 0) {
      qtl_df$gene_id = unlist(qtl_df$gene_id)
      gname_t = unlist(gname[unlist(qtl_df$gene_id)])
      gname_t = data.frame(gene_id = names(gname_t), name = gname_t)
      dname_t = unlist(dname[unlist(qtl_df$gene_id)])
      dname_t = data.frame(gene_id = names(dname_t), alias = dname_t)
      dname_t_long = unlist(dname_long[unlist(qtl_df$gene_id)])
      dname_t_long = data.frame(gene_id = names(dname_t_long), desc = dname_t_long)
      qtl_df = merge(qtl_df,gname_t,by="gene_id",sort=F,all.x=T)
      qtl_df = merge(qtl_df,dname_t,by="gene_id",sort=F,all.x=T)
      qtl_df = merge(qtl_df,dname_t_long,by="gene_id",sort=F,all.x=T)
      qtl_df = qtl_df[,c("gene_id","name","seqnames","start","end","strand","alias","desc")]
      colnames(qtl_df) = c("Sys.Name","Name","Chr","Start","End","Strand","Alias","Desc")
      #rownames(qtl_df) = qtl_df[,"Sys.Name"]
      qtl_df = qtl_df[!duplicated(qtl_df),]
      # add stitch predictions
      stitch = as.data.frame(filter(genphen_stitch,
                protein%in%unlist(qtl_df$Sys.Name),alias==input$m)%>%ungroup())[,c("protein","score")]
      colnames(stitch)[2] = "STITCH"
      qtl_df = merge(qtl_df,stitch,by.x="Sys.Name",by.y="protein",sort=F,all.x=T)
      # add snps
      #
      # TODO: break out sub categrories as pie charts or something similar
      return(qtl_df)
    }
  })

  df_var = reactive({
    qtl_df = df_full()
    var_df = filter(var_info, SNPEFF_TRANSCRIPT_ID%in%unlist(qtl_df$Sys.Name)) %>%
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
  })

  df_dt = reactive({
    # reorder
    qtl_df = merge(df_full(),df_var(),by.x="Sys.Name",by.y="SNPEFF_TRANSCRIPT_ID",sort=F,all.x=T)
    qtl_df = qtl_df[,c("Sys.Name","Name","STITCH","SNPS", "INDELS", "UPSTREAM",
                       "DOWNSTREAM", "INTRONS", "CODING","HIGH","MODERATE","LOW",
                       "Chr","Start","End","Strand","Alias","Desc")]
  })

  output$snptype = renderPlot({
    # reorder
    qtl_df = df_var()
    high = c("START_LOST","STOP_GAINED","STOP_LOST","FRAME_SHIFT")
    moderate = c("NON_SYNONYMOUS_CODING", "CODON_DELETION","CODON_INSERTION",
                 "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION")
    low = c("SYNONYMOUS_CODING","SYNONYMOUS_STOP","NON_SYNONYMOUS_START")
    var_df_melt = melt(qtl_df,id.vars="SNPEFF_TRANSCRIPT_ID") %>% filter(.,value>0,variable%in%c(high,moderate,low))
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
  })

  # DATA TABLE
  output$dt = DT::renderDataTable(
    df_dt(), server = TRUE, selection = "single",
    rownames = FALSE, extensions = 'Responsive', escape = FALSE)

  # REACTIVE VALUES
  values = reactiveValues(
    old_selection = NULL,
    link = NULL
  )

  # MONITOR OLD SELECTION
  session$onFlush(once=FALSE, function(){
    isolate({ values$old_selection <- input$dt_rows_selected })
  })

  #output$link = renderPrint("hi")
  output$link = renderText({
    s = input$dt_rows_selected[length(input$dt_rows_selected)]
    d = df_full()
    s = which(d[,"Sys.Name"]==s)[1]
    print(s)
    if (length(s)) {
      chr = levels(d[s, "Chr"])[d[s, "Chr"]]
      start = d[s, "Start"]-10000
      end = d[s, "End"]+10000
      val = paste(chr,"%3A",start,"..",end, sep = "")
      val = paste('http://steinmetzlab.embl.de/mQTL/?loc=',val,"&tracks=",input$m,"%2C", d[s, "Sys.Name"], sep="")
      print(val)
    } else {
      val = "http://steinmetzlab.embl.de/mQTL/"
    }
    paste("<div style='width: 100%; height: 600px'><iframe style='border: 1px solid black' src='", val ,"'width='100%' height='100%'></iframe></div>",sep="")
  })



  alpha_10 = reactive({
    #summary(data[[input$m]]$permout[,type],input$co/100)
    as.integer(input$co)[1]
  })

   output$manhattan = renderText(
   {
   paste0('<td align="middle"><img src="http://steinmetzlab.embl.de/GenPhen/mQTL_plots/',
    input$m, "/mQTL_", input$m, "_FDR_", alpha_10(), '.png" ',
    'valign="middle" style="width: 100%;max-height: 100%"></td>')
   })


})
