# The purpose of this script is to:
# - join records from the sample into a single column
# - set the COG ID as an indicator name for components 
# - export the count table and the ortholog structure as a text format
#
# This script requires a .biom output file that includes all samples/records of an experiment.
# To create a .biom output file in MEGAN: 
#   1) loading all datasets into one Megan explorer, (using the compare datasets function).
#   2) selecting/highlight the full ortholog structure.
#   3) click File->Export-> as Biom File. 
# Generating this output might take a moment. 
#
# The second part of this script also include the possibility to read and join the count 
# table from multiple .txt/.csv court sheets. However, this procedure is much more complicated and 
# requires a secondary matching of the components to a ortholog structure export to recover the ortholog 
# information. 
# 


####################################################################################
# -- General
library(phyloseq) # required for load via .biom file

source("tools.R")

## set paths
# input directory
path.in <- './Data'




####################################################################################
# -- Export (and join) count table from .biom file.  

# load file
otu <- Biome2Table(paste(path.in, "/readcounts/Comparison-EGGNOG.biom",sep=""), 
                    sample.reg="ID-[0-9]+-", rename="[0-9]+", sort.rows = T, sort.cols = T )

# remove empty rows
otu <- otu[rowSums(otu[,1:26])>0,]

# set ortholog ID as row name
c.names <- stringr::str_extract(otu[,ncol(otu)], "^(\\w+)")
c.names[is.na(c.names)] <- "Not_assigned"
rownames(otu) <- c.names
# sort alphabetically
otu <- otu[order(rownames(otu)),]

# split count table from ortholog structure
ds <- otu[,1:26]
ortholog <- otu[,27:30]
# rename ortholog columns
colnames(ortholog) <- c("L0", "L1", "L2", "L3")

# write count table to file 
write.table(ds, file=paste(path.in,"MT-OTU-all.txt",sep = "/"), 
            quote = F, sep='\t', eol='\n', na='NA', row.names = T, col.names = T)
# write ortholog table to file
write.table(ortholog, file=paste(path.in,"MT-ortholog-all.txt",sep = "/"), 
            quote = F, sep='\t', eol='\n', na='NA', row.names = T, col.names = T)




####################################################################################
# -- Alternative approach: compute count table from (multiple) count sheets.


# # set regular expressions to search for in folder and record name. optional
# re.levels <- c(".*level1\\.txt$", ".*level2\\.txt$", ".*level3\\.txt$")
# re.samples <- ".*?([0-9]+).*"
# # specify any column to drop
# drop.col <- NA
# 
# # get input file from level 3
# files <- list.files(path = paste(path.in,'/readcounts',sep='') , 
#                    full.names = TRUE, recursive = TRUE, pattern= re.levels[3] )
# 
# # create count table 
# ds <- ReadCountTable(files = files, re.samples = re.samples, re.ID.separator = " ")
# # sort count table by sample number 
# ds <- ds[,order(as.numeric(names(ds)))]
# 
# # write count table to file 
# write.table(ds, file=paste(path.in,"MT-OTU-L3-all.txt",sep = "/"), 
#             quote = F, sep='\t', eol='\n', na='NA', row.names = T, col.names = T)









