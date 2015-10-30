`M2.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/M2.genes.results")
`M3.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/M3.genes.results")
`M4.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/M4.genes.results")
`M5.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/M5.genes.results")
`LO2.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/LO2.genes.results")
`LO3.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/LO3.genes.results")
`LO4.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/LO4.genes.results")
`LO5.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/fish/LO5.genes.results")
Genes<- data.frame(M2.genes$gene_id,M2.genes$expected_count,M3.genes$expected_count,M4.genes$expected_count,M5.genes$expected_count, LO2.genes$expected_count,LO3.genes$expected_count, LO4.genes$expected_count,LO5.genes$expected_count,LO5.genes$expected_count)

View (Genes)

names(Genes)<-c("gene_id","M2","M3","M4", "M5","LO2","LO3","LO4", "LO5", "LO5a")

## write data files to csv files
#write.csv(Genes, file="/Users/brisbin/desktop/Siphamia/siph.genes.results/Genes.csv")
#write.csv(FPKM, file="/Users/brisbin/desktop/Siphamia/siph.genes.results/FPKM.csv")

FPKM<- data.frame(M2.genes$gene_id,M2.genes$FPKM,M3.genes$FPKM,M4.genes$FPKM, M5.genes$FPKM, LO2.genes$FPKM, LO3.genes$FPKM,LO4.genes$FPKM,  LO5.genes$FPKM )
ercc_conc<-cms_095046 <- read.delim("~/Desktop/Siphamia/siph.genes.results/cms_095046.txt", stringsAsFactors=FALSE)
ercc_conc<-ercc_conc[,2:7]
ercc_conc<-ercc_conc[order(ercc_conc$ERCC.ID),]
ERCC<-FPKM[1:92,]

names(Genes)<-c("gene_id","M2","M3","M4", "M5","LO2","LO3","LO4", "LO5", "LO5b")
names(ERCC)<-c("gene_id","M2","M3","M4", "M5","LO2","LO3","LO4", "LO5" )

#Make table for Mix 1 and Mix2
Mix1<- data.frame(ERCC$gene_id, ERCC$LO2, ERCC$LO5, ERCC$M3, ERCC$M4)
Mix2<- data.frame(ERCC$gene_id, ERCC$LO3, ERCC$LO4, ERCC$M2, ERCC$M5)
# make a column for average of sample FPKM for each row in each mix's dataframs
Mix1$Mix1avg<-rowMeans(Mix1[2:5])
Mix2$Mix2avg<-rowMeans(Mix2[2:5])
# add the mix two average column to the mix 1 data frame
Mix1$Mix2avg<-Mix2$Mix2avg
#add column with expected ratios for ech gene from ERCC website
Mix1$control<-ercc_conc$log2.Mix.1.Mix.2.
# filter any gene with less than 1 FPKM
Mix1<- Mix1[Mix1$Mix1avg>=1,]
Mix1<-Mix1[Mix1$Mix2avg>=1,]

# calculate the observed ratio by dividing Mix1avg by Mix2avg  
Mix1$ratio<- Mix1$Mix1avg / Mix1$Mix2avg
#take log of the ratio
Mix1$logratio<-log2(Mix1$ratio)
# filter data frame to remove undefined values
Mix1<-Mix1[is.finite(Mix1$logratio), ]


library("ggplot2", lib.loc="~/Desktop/R-3.2.2/library")

expectobserve<-ggplot(Mix1, aes(x=control,y=logratio)) +geom_point()+ggtitle(expression(paste("Log"[2]," mix1/mix2", sep="")))+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," expected ratio", sep=""))) +
  labs(y=expression(paste("Log"[2]," observed ratio", sep=""))) +
  scale_x_continuous(limits = c(-2, 2)) 


## Make data frames with MIX 1 and MIX 2 FPKM from **ONE sample TYPE** at a time 
Mix1lo<- data.frame(ERCC$gene_id, ERCC$LO2, ERCC$LO5)
Mix2lo<- data.frame(ERCC$gene_id, ERCC$LO3, ERCC$LO4)
Mix1m<- data.frame(ERCC$gene_id, ERCC$M3, ERCC$M4)
Mix2m<- data.frame(ERCC$gene_id, ERCC$M2, ERCC$M5)
# make a column for average of sample FPKM for each row in each mix's dataframs
Mix1lo$Mix1avg<-rowMeans(Mix1lo[2:3])
Mix2lo$Mix2avg<-rowMeans(Mix2lo[2:3])
Mix1m$Mix1avg<-rowMeans(Mix1m[2:3])
Mix2m$Mix2avg<-rowMeans(Mix2m[2:3])
# add the mix two average column to the mix 1 data frame
Mix1lo$Mix2avg<-Mix2lo$Mix2avg
Mix1m$Mix2avg<-Mix2m$Mix2avg
#add column with expected ratios for ech gene from ERCC website
Mix1lo$control<-ercc_conc$log2.Mix.1.Mix.2.
Mix1m$control<-ercc_conc$log2.Mix.1.Mix.2.
# filter any gene with less than 1 FPKM
Mix1m<- Mix1m[Mix1m$Mix1avg>=10,]
Mix1m<-Mix1m[Mix1m$Mix2avg>=10,]
Mix1lo<- Mix1lo[Mix1lo$Mix1avg>=10,]
Mix1lo<-Mix1lo[Mix1lo$Mix2avg>=10,]
# calculate the observed ratio by dividing Mix1avg by Mix2avg  
Mix1m$ratio<- Mix1m$Mix1avg / Mix1m$Mix2avg
Mix1lo$ratio<- Mix1lo$Mix1avg / Mix1lo$Mix2avg
#take log of the ratio
Mix1m$logratio<-log2(Mix1m$ratio)
Mix1lo$logratio<-log2(Mix1lo$ratio)
# filter data frame to remove undefined values
Mix1m1<-Mix1m[is.finite(Mix1m$logratio), ]
Mix1lo1<-Mix1lo[is.finite(Mix1lo$logratio), ]

loANDm<- merge(Mix1lo1, Mix1m1, by ="ERCC.gene_id")
keeps<- c("control.x","logratio.x","logratio.y")
loANDm<- loANDm[keeps]
names(loANDm)<-c("expected", "LO", "M")

  
library("reshape2", lib.loc="~/Desktop/R-3.2.2/library")

mloANDm <- melt(loANDm, id=c("expected"))

library("ggplot2", lib.loc="~/Desktop/R-3.2.2/library")

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

ggplot(mloANDm, aes(x=expected, y=value)) + 
  geom_point(aes(colour=variable)) +ggtitle(expression(paste("Log"[2]," mix1/mix2 (FPKM filter <10)", sep="")))+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  labs(x=expression(paste("Log"[2]," expected ratio", sep=""))) +
  labs(y=expression(paste("Log"[2]," observed ratio", sep=""))) +
  geom_smooth(aes(expected,value, colour=variable), method=lm, se=FALSE)+
  scale_colour_manual(values=cbPalette)+
  scale_x_continuous(limits = c(-1, 2))+
  scale_y_continuous(limits = c(-2, 4.5))
  





#Make data frame with just sample FPKM values for ERCC genes and concentration of ERCC genes in spike in mix
# make sure items in data frame are integers

#LO2 mix 1
LO2spkindf<-data.frame(ERCC$LO2,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(LO2spkindf)<-c("FPKM", "conc")
LO2spkindf$FPKM<-as.integer(LO2spkindf$FPKM)
LO2spkindf$conc<-as.integer(LO2spkindf$conc)
row_sub = apply(LO2spkindf, 1, function(row) all(row !=0 ))
LO2spkindf=LO2spkindf[row_sub,]
LO2spkindf$FPKM<-log2(LO2spkindf$FPKM)
LO2spkindf$conc<-log2(LO2spkindf$conc)

#LO3 mix 2
LO3spkindf<-data.frame(ERCC$LO3,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(LO3spkindf)<-c("FPKM", "conc")
LO3spkindf$FPKM<-as.integer(LO3spkindf$FPKM)
LO3spkindf$conc<-as.integer(LO3spkindf$conc)
row_sub = apply(LO3spkindf, 1, function(row) all(row !=0 ))
LO3spkindf=LO3spkindf[row_sub,]
LO3spkindf$FPKM<-log2(LO3spkindf$FPKM)
LO3spkindf$conc<-log2(LO3spkindf$conc)

#LO4 mix 2
LO4spkindf<-data.frame(ERCC$LO4,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(LO4spkindf)<-c("FPKM", "conc")
LO4spkindf$FPKM<-as.integer(LO4spkindf$FPKM)
LO4spkindf$conc<-as.integer(LO4spkindf$conc)
row_sub = apply(LO4spkindf, 1, function(row) all(row !=0 ))
LO4spkindf=LO4spkindf[row_sub,]
LO4spkindf$FPKM<-log2(LO4spkindf$FPKM)
LO4spkindf$conc<-log2(LO4spkindf$conc)

#LO5 mix 1
LO5spkindf<-data.frame(ERCC$LO5,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(LO5spkindf)<-c("FPKM", "conc")
LO5spkindf$FPKM<-as.integer(LO5spkindf$FPKM)
LO5spkindf$conc<-as.integer(LO5spkindf$conc)
row_sub = apply(LO5spkindf, 1, function(row) all(row !=0 ))
LO5spkindf=LO5spkindf[row_sub,]
LO5spkindf$FPKM<-log2(LO5spkindf$FPKM)
LO5spkindf$conc<-log2(LO5spkindf$conc)

#M2 mix 2
M2spkindf<-data.frame(ERCC$M2,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(M2spkindf)<-c("FPKM", "conc")
M2spkindf$FPKM<-as.integer(M2spkindf$FPKM)
M2spkindf$conc<-as.integer(M2spkindf$conc)
row_sub = apply(M2spkindf, 1, function(row) all(row !=0 ))
M2spkindf=M2spkindf[row_sub,]
M2spkindf$FPKM<-log2(M2spkindf$FPKM)
M2spkindf$conc<-log2(M2spkindf$conc)

#M3 mix 1
M3spkindf<-data.frame(ERCC$M3,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(M3spkindf)<-c("FPKM", "conc")
M3spkindf$FPKM<-as.integer(M3spkindf$FPKM)
M3spkindf$conc<-as.integer(M3spkindf$conc)
row_sub = apply(M3spkindf, 1, function(row) all(row !=0 ))
M3spkindf=M3spkindf[row_sub,]
M3spkindf$FPKM<-log2(M3spkindf$FPKM)
M3spkindf$conc<-log2(M3spkindf$conc)

#M4 mix 1
M4spkindf<-data.frame(ERCC$M4,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(M4spkindf)<-c("FPKM", "conc")
M4spkindf$FPKM<-as.integer(M4spkindf$FPKM)
M4spkindf$conc<-as.integer(M4spkindf$conc)
row_sub = apply(M4spkindf, 1, function(row) all(row !=0 ))
M4spkindf=M4spkindf[row_sub,]
M4spkindf$FPKM<-log2(M4spkindf$FPKM)
M4spkindf$conc<-log2(M4spkindf$conc)

#M5 mix 2
M5spkindf<-data.frame(ERCC$M5,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(M5spkindf)<-c("FPKM", "conc")
M5spkindf$FPKM<-as.integer(M5spkindf$FPKM)
M5spkindf$conc<-as.integer(M5spkindf$conc)
row_sub = apply(M5spkindf, 1, function(row) all(row !=0 ))
M5spkindf=M5spkindf[row_sub,]
M5spkindf$FPKM<-log2(M5spkindf$FPKM)
M5spkindf$conc<-log2(M5spkindf$conc)

#function creating text for line equation and r squared value
#enter name of data frame
lm_eqn <- function(df){
  y<-df$FPKM
  x<-df$conc
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

# plot ERCC genes FPKM values v. the concentration of genes from ERCC spikein
#use Log2 x and y axis
#add line equation and rsquared to the plot
## tried to add linear regression line and could not get past the error messages

library("ggplot2", lib.loc="~/Desktop/R-3.2.2/library")
#LO2
spikeinPlotLO2<-ggplot(LO2spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO2/spike-in 1")+ geom_smooth(method=lm, se=FALSE, color="black") +
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black'))+ geom_text(x = 6, y = 14, label = lm_eqn(LO2spkindf), parse = TRUE) +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO3
spikeinPlotLO3<-ggplot(LO3spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO3/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(LO3spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

# LO4
spikeinPlotLO4<-ggplot(LO4spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO4/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(LO4spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO5
spikeinPlotLO5<-ggplot(LO5spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO5/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(LO5spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#M2
spikeinPlotM2<-ggplot(M2spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("M2/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(B29spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#M3
spikeinPlotM3<-ggplot(M3spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("M3/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(B30spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#M4
spikeinPlotM4<-ggplot(M4spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("M4/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(B31spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#M5
spikeinPlotM5<-ggplot(M5spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("M5/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 6, y = 14, label = lm_eqn(B32spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#put plots together
library("gridExtra", lib.loc="~/Desktop/R-3.2.2/library")
grid.arrange(spikeinPlotLO2, spikeinPlotLO3, spikeinPlotLO4, spikeinPlotLO5, ncol=2, nrow=2)
grid.arrange(spikeinPlotM2, spikeinPlotM3, spikeinPlotM4, spikeinPlotM5, ncol=2, nrow=2)


#get list of genes that average FPKM across all samples  <1
# will make life easier if remove from count df at the same time!
# add column to df FPKM that is the average of all other rows
# remove all rows that are <1 in the average column
#export column with gene IDs as a text file?



## getting low expressed gene list! for python script to further filter IsoformFiltered_Trinity.fasta to FPKMfiltered_Trinity.fasta
# data frame: geneid, gene counts for all samples, FPKM for all samples
Counts_n_FPKM<- data.frame(M2.genes$gene_id,M2.genes$expected_count,M3.genes$expected_count,M4.genes$expected_count,M5.genes$expected_count, LO2.genes$expected_count,LO3.genes$expected_count, LO4.genes$expected_count,LO5.genes$expected_count,M2.genes$FPKM,M3.genes$FPKM,M4.genes$FPKM, M5.genes$FPKM, LO2.genes$FPKM, LO3.genes$FPKM,LO4.genes$FPKM,  LO5.genes$FPKM )
#remove ERCC genes from dataframe 
Counts_n_FPKM<-Counts_n_FPKM[93:239373,]
# add column that is the average of FPKM for all samples
Counts_n_FPKM$FPKMavg<-rowMeans(Counts_n_FPKM[10:17])
#remove all rows that have a value less than 1 in column FPKMavg
Counts_n_FPKM_new <- Counts_n_FPKM[Counts_n_FPKM$FPKMavg>1,]
#now only 120,026 genes!! (down from 239,281)
#make list of just gene ids
genelist<- data.frame(Counts_n_FPKM_new$M2.genes.gene_id)
names(genelist)<-c("gene_id")
# write genelist to a csv
write.csv(genelist, file="/Users/brisbin/desktop/Siphamia/siph.genes.results/filtered_1FPKM_genelist.csv")




#Load edgeR
library("edgeR", lib.loc="~/Desktop/R-3.2.2/library")
#set working directory
setwd("/Users/brisbin/desktop/Siphamia/edgeR_wd")


#remove ERCC genes from dataframe Genes
Genes<-Genes[93:239373,]

#read in data
counts <- Genes[ , -c(1,ncol(Genes)) ]
rownames(counts)<-Genes[ , 1] #gene names
colnames(counts)<-paste(c("M2","M3","M4","M5","LO2","LO3","LO4", "LO5"))
View(counts)

##Basic info about the data
#size of data frame
dim(Genes)
#column sums - library size
colSums(counts)
#library size in millions
colSums(counts)/1e06
# Number of genes with low counts
table( rowSums( counts ) )[ 1:30 ] 

#convert count matrix to edgeR Object
#create a group variable that tells edgeR which samples belong to which group
group<- c(rep("M",4), rep("LO",4))
cds<-DGEList(counts, group=group)
names(cds)
#original count data
head(cds$counts)
#contains sample information
head(cds$samples)
# How many genes have 0 counts across all samples
sum( cds$all.zeros )


##need to filter out low count reads bc impossible to detect differential expression
#keep only genes with atleast 1 read per million reads in at least 3 samples
#once this is done, can calculate normalisation factors which correct for different compositions of samples
# effective library size = product of actual library size and these factors

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
dim( cds )
cds <- calcNormFactors( cds )
cds$samples

# effective library sizes
cds$samples$lib.size * cds$samples$norm.factors

#Multi-Dimensional Scaling Plot
#measures similarity of samples and projects measure into 2 dimensions
plotMDS( cds , main = "MDS Plot for Fish Count Data", labels = colnames( cds$counts ) )

plotMD(cpm(cds, log=TRUE), column=8)
abline(h=0, col="red", lty=2, lwd=2)

#Estimating Dispersions
#1st: calculate common dispersion
#each gene get assigned the same dispersion estimate
#output of the estimation includes estimate andother elements added to object cds
cds <- estimateCommonDisp( cds )
names( cds )
#The estimate
cds$common.dispersion

#with common dispersion, can estimate tagwise dispersions
#each gene will get its own dispersion estimate
#tagwise dispersions are squeezed towards a common value
#amt of squeezing is governed by parameter prior.n
#higher prior.n --> closer the estimates will be to common dispersion
#default value is nearest interger to 50/(#samples-#groups)
# 50/(6-2)= 12.5 -->13

cds <- estimateTagwiseDisp( cds)
names( cds )
summary( cds$tagwise.dispersion )

#Mean-Variance Plot
#see how well the dispersion factors fir the data
# see raw variance of counts (grey dots), 
#the variance using tagwise dispersions(light blue dots)
#variances using common dispersion (solid blue line)
#variance=mean aka poisson variance (solid black line)

meanVarPlot <- plotMeanVar( cds , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

#Testing
# exactTest() performs pairwise tests for diff exp between 2 groups
# if common.disp parameter is true, uses common dispersion estimate for testing
# if false, uses tagwise dispersion estimates
# pair indicates which groups should be compared
#output is a list of elements, one of which is a table of results
# de.poi sets dispersion parameter to zero - poisson model - can compare negative binomial results to this

de.cmn <- exactTest( cds , dispersion="common", pair = c( "M" , "LO" ) )
de.tgw <- exactTest( cds , dispersion="tagwise", pair = c( "M" , "LO" ) )
de.poi <- exactTest( cds , dispersion = 1e-06 , pair = c( "M" , "LO" ) )

names( de.tgw )
de.tgw$comparison # which groups have been compared
head( de.tgw$table ) # results table in order of your count matrix.
head( cds$counts )

# this does not contain p-values adjusted for multiple testin
# topTags() takes output from exactTest(), adjusts raw p-values using False Discover Rate (FDR) correction
# returns top differentially expressed genes
# sort.by argument allws sorting by p-value, concentration, or fold change
# topTags() returns original counts of the top differentially expressed genes
# set n parameter to number of genes, save the entire topTags() results table

# Top tags for tagwise analysis
options( digits = 3 ) # print only 3 digits
topTags( de.tgw , n = 20 , sort.by = "p.value" ) # top 20 DE genes
# Back to count matrix for tagwise analysis
cds$counts[ rownames( topTags( de.tgw , n = 15 )$table ) , ]
# Sort tagwise results by Fold-Change instead of p-value
resultsByFC.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) , sort.by = "logFC" )$table
head( resultsByFC.tgw )
# Store full topTags results table
resultsTbl.cmn <- topTags( de.cmn , n = nrow( de.cmn$table ) )$table
resultsTbl.tgw <- topTags( de.tgw , n = nrow( de.tgw$table ) )$table
resultsTbl.poi <- topTags( de.poi , n = nrow( de.poi$table ) )$table
head( resultsTbl.tgw )

#compare p-values to significance level and determine # of diff exp genes
# sig level = .05
# decideTestsDGE() to determine how many diff exp genes are up or down regulated compared to control

# Names/IDs of DE genes
de.genes.cmn <- rownames( resultsTbl.cmn )[ resultsTbl.cmn$PValue <= 0.05 ]
de.genes.tgw <- rownames( resultsTbl.tgw )[ resultsTbl.tgw$PValue <= 0.05 ]
de.genes.poi <- rownames( resultsTbl.poi )[ resultsTbl.poi$PValue <= 0.05 ]
# Amount significant
length( de.genes.cmn )
length( de.genes.tgw )
length( de.genes.poi )
# Percentage of total genes
length( de.genes.cmn ) / nrow( resultsTbl.cmn ) * 100
length( de.genes.tgw ) / nrow( resultsTbl.tgw ) * 100
length( de.genes.poi ) / nrow( resultsTbl.poi ) * 100
# Up/Down regulated summary for tagwise results
summary( decideTestsDGE( de.tgw , p.value = 0.05 ) ) # the adjusted p-values are used here

# compare analyses to see how many diff exp genes are in common


#MA plot shows relationship between concentration and fold-change across genes
# DE genes are red
# non DE genes are black
# organge dots = counts were zero in all samples in one group
# blue line is at log-FC of 2 to represent a level for biological significance

par( mfrow=c(2,1) )
plotSmear( cds , de.tags=de.genes.poi , main="Fish Poisson MA plot" ,
           pair = c("M","LO") ,
           cex = .35 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h = c(-2, 2) , col = "dodgerblue" )
plotSmear( cds , de.tags=de.genes.tgw , main="Fish Tagwise MA plot" ,
           pair = c("M","LO") ,
           cex = .35 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
par( mfrow=c(1,1) )

# log difference between DE genes under negative binomial model and poisson model in MA plot
# plot top 500 DE genes are red, rest are black

par( mfrow = c(2,1) )
plotSmear( cds , de.tags=de.genes.poi[1:500] , main="Fish Poisson MA plot" ,
           pair=c("M","LO") ,
           cex=.5 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
plotSmear( cds , de.tags=de.genes.tgw[1:500] , main="Fish Tagwise MA plot, top 500 DE" ,
           pair=c("M","LO") ,
           cex = .5 ,
           xlab="Log Concentration" , ylab="Log Fold-Change" )
abline( h=c(-2,2) , col="dodgerblue" )
par( mfrow=c(1,1) )

#Output results
#make a table or csv file containing results with concentrations, fold-changes, p-values
#up/down regulated variable, dispersions, and the count matrix

# Change column names to be specific to the analysis, logConc and logFC are the same in both.
colnames( resultsTbl.cmn ) <- c( "logConc" , "logFC" , "pVal.Cmn" , "adj.pVal.Cmn" )
colnames( resultsTbl.tgw ) <- c( "logConc" , "logFC" , "pVal.Tgw" , "adj.pVal.Tgw" )
# Below provides the info to re-order the count matrix to be in line with the order of the results.
wh.rows.tgw <- match( rownames( resultsTbl.tgw ) , rownames( cds$counts ) )
wh.rows.cmn <- match( rownames( resultsTbl.cmn ) , rownames( cds$counts ) )
head( wh.rows.tgw )
# Tagwise Results
combResults.tgw <- cbind( resultsTbl.tgw ,
                          "Tgw.Disp" = cds$tagwise.dispersion[ wh.rows.tgw ] ,
                          "UpDown.Tgw" = decideTestsDGE( de.tgw , p.value = 0.05 )[ wh.rows.tgw ] ,
                          cds$counts[ wh.rows.tgw , ] )
head( combResults.tgw )
# Common Results
combResults.cmn <- cbind( resultsTbl.cmn ,
                          "Cmn.Disp" = cds$common.dispersion ,
                          "UpDown.Cmn" = decideTestsDGE( de.cmn , p.value = 0.05 )[ wh.rows.cmn ] ,
                          cds$counts[ wh.rows.cmn , ] )
head( combResults.cmn )

#combine both common and tagwise results together

wh.rows <- match( rownames( combResults.cmn ) , rownames( combResults.tgw ) )
combResults.all <- cbind( combResults.cmn[,1:4] ,
                          combResults.tgw[wh.rows,3:4] ,
                          "Cmn.Disp" = combResults.cmn[,5],
                          "Tgw.Disp" = combResults.tgw[wh.rows,5],
                          "UpDown.Cmn" = combResults.cmn[,6],
                          "UpDown.Tgw" = combResults.tgw[wh.rows,6],
                          combResults.cmn[,7:ncol(combResults.cmn)] )
head( combResults.all )
# Ouput csv tables of results
write.table( combResults.tgw , file = "combResults_tgw_ex1.csv" , sep = "," , row.names = TRUE )
write.table( combResults.cmn , file = "combResults_cmn_ex1.csv" , sep = "," , row.names = TRUE )

#Visualize results
#spread of expression levels for DE genes
#plot concentrations for top 100 DE genes for each analysis

par( mfrow=c(3 ,1) )
hist( resultsTbl.poi[de.genes.poi[1:100],"logConc"] , breaks=10 , xlab="Log Concentration" ,
      col="red" , xlim=c(-18,-6) , ylim=c(0,0.4) , freq=FALSE , main="Poisson: Top 100" )

hist( resultsTbl.cmn[de.genes.cmn[1:100],"logConc"] , breaks=25 , xlab="Log Concentration" ,
      col="green" , xlim=c(-18,-6) , ylim=c(0,0.4) , freq=FALSE , main="Common: Top 100" )

hist( resultsTbl.tgw[de.genes.tgw[1:100],"logConc"] , breaks=25 , xlab="Log Concentration" ,
      col="blue" , xlim=c(-18,-6) , ylim=c(0,0.4) , freq=FALSE , main="Tagwise: Top 100" )
par( mfrow=c(1,1) )
write.table( combResults.all , file = "combResults_all_ex1.csv" , sep = "," , row.names = TRUE )