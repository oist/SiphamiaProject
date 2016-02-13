library(ggplot2)
library(GOstats)
library(GSEABase)

#read in go gene reference
go <-read.csv("~/Desktop/geneGOref.csv", header=FALSE)
# add column - GOstats expects this column to be there
go$evidence <- "ISS"
#add column names
names(go) <- c("isoform", "GO", "evidence")
#reorder columns
go <- go[,c("GO","evidence","isoform")]
#GOframe is a function in GOstats - did not call only unique isoforms or from a particular 
#organism because there are not repeats here and everything should be P.mandapamensis
goframe=GOFrame(go)
# commands from Sasha 
goAllFrame=GOAllFrame(goframe)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
