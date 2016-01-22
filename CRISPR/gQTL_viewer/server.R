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
    return(data.frame())
  })

  df_var = reactive({
    qtl_df = df_full()
    if (dim(qtl_df)[1]>0) {
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
    } else {
      return(data.frame())
    }
  })

  df_dt = reactive({
    # reorder
    df_full = df_full()
    if (dim(df_full)[1]>0) {
      qtl_df = merge(df_full,df_var(),by.x="Sys.Name",by.y="SNPEFF_TRANSCRIPT_ID",sort=F,all.x=T)
      qtl_df = qtl_df[,c("Sys.Name","Name","STITCH","SNPS", "INDELS", "UPSTREAM",
                         "DOWNSTREAM", "INTRONS", "CODING","HIGH","MODERATE","LOW",
                         "Chr","Start","End","Strand","Alias","Desc")]
    } else {
      return(data.frame())
    }
  })

  output$snptype = renderPlot({
    # reorder
    qtl_df = df_var()
    if(dim(qtl_df)[1]>0) {
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
    } else {
      return(NULL)
    }
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

  output$image <- renderText({
   paste0('<td align="middle"><img src="http://scalefreegan.github.io/steinmetz-lab/CRISPR/data/ge_vals.png" ',
           'valign="middle" style="width: 40%;max-height: 100%"></td>')
    })
  
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

})
