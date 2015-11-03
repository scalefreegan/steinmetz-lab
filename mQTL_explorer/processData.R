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

# Import packages ---------------------------------------------------
library(funqtl)

# tmp resource location / will be changed
DDIR = "/Users/brooks/Documents/steinmetz_local/genphen/metabolome"

load(file.path(DDIR,"endometabolite_full_12102015.rda"))
load(file.path(DDIR,"genotypes_S288c_R64.rda"))
load(file.path(DDIR,"mQTLs_comball_funqtl_2014.rda"))

data = lapply(seq(1,length(mQTLs_funqtl_2014)), function(i){
  if (class(mQTLs_funqtl_2014[[i]])=="try-error") {
    return(NULL)
  } else {
    #print(i)
    o = list()
    o[["qtl"]] = mQTLs_funqtl_2014[[i]][["qtls_alt"]]
    o[["permout"]] = mQTLs_funqtl_2014[[i]]$permout
    return(o)
  }
})
names(data) = names(mQTLs_funqtl_2014)

# remove NULLs
data = data[!sapply(data,is.null)]
