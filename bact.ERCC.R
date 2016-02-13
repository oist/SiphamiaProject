# load data for bacteria DE
`B29.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/B29.genes.results")
`B30.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/B30.genes.results")
`B31.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/B31.genes.results")
`B32.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/B32.genes.results")
`LO2.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/LO2.genes.results")
`LO3.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/LO3.genes.results")
`LO4.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/LO4.genes.results")
`LO5.genes` <- read.delim("~/Desktop/Siphamia/siph.genes.results/bacteria/LO5.genes.results")

#make data frames
Genes<- data.frame(B29.genes$gene_id,B29.genes$expected_count,B30.genes$expected_count, B31.genes$expected_count,B32.genes$expected_count, LO2.genes$expected_count,LO3.genes$expected_count, LO4.genes$expected_count,LO5.genes$expected_count,LO5.genes$expected_count)
FPKM<- data.frame(B29.genes$gene_id,B29.genes$FPKM,B30.genes$FPKM,B31.genes$FPKM, B32.genes$FPKM, LO2.genes$FPKM, LO3.genes$FPKM,LO4.genes$FPKM,  LO5.genes$FPKM )
ercc_conc<-cms_095046 <- read.delim("~/Desktop/Siphamia/siph.genes.results/cms_095046.txt", stringsAsFactors=FALSE)
ercc_conc<-ercc_conc[,2:7]
ercc_conc<-ercc_conc[order(ercc_conc$ERCC.ID),]
ERCC<-FPKM[1:92,]

## ERCC Expected v. Observed RATIOs
FPKM<- data.frame(B29.genes$gene_id,B29.genes$FPKM,B30.genes$FPKM,B31.genes$FPKM, B32.genes$FPKM, LO2.genes$FPKM, LO3.genes$FPKM,LO4.genes$FPKM,  LO5.genes$FPKM )
ercc_conc<-cms_095046 <- read.delim("~/Desktop/Siphamia/siph.genes.results/cms_095046.txt", stringsAsFactors=FALSE)
ercc_conc<-ercc_conc[,2:7]
ercc_conc<-ercc_conc[order(ercc_conc$ERCC.ID),]
ERCC<-FPKM[1:92,]

names(Genes)<-c("gene_id","B29","B30","B31", "B32","LO2","LO3","LO4", "LO5", "LO5b")
names(ERCC)<-c("gene_id","B29","B30","B31", "B32","LO2","LO3","LO4", "LO5" )

#Make table for Mix 1 and Mix2
Mix1<- data.frame(ERCC$gene_id, ERCC$LO2, ERCC$LO5, ERCC$B29, ERCC$B32)
Mix2<- data.frame(ERCC$gene_id, ERCC$LO3, ERCC$LO4, ERCC$B30, ERCC$B31)
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

expectobserve<-ggplot(Mix1, aes(x=control,y=logratio)) +geom_point()+ggtitle(expression(paste("Bacteria Log"[2]," mix1/mix2 (FPKM <1 filtered)", sep="")))+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," expected ratio", sep=""))) +
  labs(y=expression(paste("Log"[2]," observed ratio", sep=""))) +
  scale_x_continuous(limits = c(-2, 2)) +
  geom_text(x = -1, y = 3, label = lm_eqn(Mix1), parse = TRUE)


## Make data frames with MIX 1 and MIX 2 FPKM from **ONE sample TYPE** at a time 
Mix1b<- data.frame(ERCC$gene_id, ERCC$B29, ERCC$B32)
Mix2b<- data.frame(ERCC$gene_id, ERCC$B30, ERCC$B31)
Mix1lo<- data.frame(ERCC$gene_id, ERCC$LO2, ERCC$LO5)
Mix2lo<- data.frame(ERCC$gene_id, ERCC$LO3, ERCC$LO4)
# make a column for average of sample FPKM for each row in each mix's dataframs
Mix1lo$Mix1avg<-rowMeans(Mix1lo[2:3])
Mix2lo$Mix2avg<-rowMeans(Mix2lo[2:3])
Mix1b$Mix1avg<-rowMeans(Mix1b[2:3])
Mix2b$Mix2avg<-rowMeans(Mix2b[2:3])
# add the mix two average column to the mix 1 data frame
Mix1lo$Mix2avg<-Mix2lo$Mix2avg
Mix1b$Mix2avg<-Mix2b$Mix2avg
#add column with expected ratios for ech gene from ERCC website
Mix1lo$control<-ercc_conc$log2.Mix.1.Mix.2.
Mix1b$control<-ercc_conc$log2.Mix.1.Mix.2.
# filter any gene with less than 1 FPKM
Mix1b<- Mix1b[Mix1b$Mix1avg>=10,]
Mix1b<-Mix1b[Mix1b$Mix2avg>=10,]
Mix1lo<- Mix1lo[Mix1lo$Mix1avg>=10,]
Mix1lo<-Mix1lo[Mix1lo$Mix2avg>=10,]
# calculate the observed ratio by dividing Mix1avg by Mix2avg  
Mix1b$ratio<- Mix1b$Mix1avg / Mix1b$Mix2avg
Mix1lo$ratio<- Mix1lo$Mix1avg / Mix1lo$Mix2avg
#take log of the ratio
Mix1b$logratio<-log2(Mix1b$ratio)
Mix1lo$logratio<-log2(Mix1lo$ratio)
# filter data frame to remove undefined values
Mix1b1<-Mix1b[is.finite(Mix1b$logratio), ]
Mix1lo1<-Mix1lo[is.finite(Mix1lo$logratio), ]

loANDb<- merge(Mix1lo1, Mix1b1, by ="ERCC.gene_id")
keeps<- c("control.x","logratio.x","logratio.y")
loANDb<- loANDb[keeps]
names(loANDb)<-c("expected", "LO", "B")


library("reshape2", lib.loc="~/Desktop/R-3.2.2/library")

mloANDb <- melt(loANDb, id=c("expected"))

library("ggplot2", lib.loc="~/Desktop/R-3.2.2/library")

cbPalette <- c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73")

ggplot(mloANDb, aes(x=expected, y=value)) + 
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

#b29 mix 1
B29spkindf<-data.frame(ERCC$B29,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(B29spkindf)<-c("FPKM", "conc")
B29spkindf$FPKM<-as.integer(B29spkindf$FPKM)
B29spkindf$conc<-as.integer(B29spkindf$conc)
row_sub = apply(B29spkindf, 1, function(row) all(row !=0 ))
B29spkindf=B29spkindf[row_sub,]
B29spkindf$FPKM<-log2(B29spkindf$FPKM)
B29spkindf$conc<-log2(B29spkindf$conc)

#b30 mix 2
B30spkindf<-data.frame(ERCC$B30,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(B30spkindf)<-c("FPKM", "conc")
B30spkindf$FPKM<-as.integer(B30spkindf$FPKM)
B30spkindf$conc<-as.integer(B30spkindf$conc)
row_sub = apply(B30spkindf, 1, function(row) all(row !=0 ))
B30spkindf=B30spkindf[row_sub,]
B30spkindf$FPKM<-log2(B30spkindf$FPKM)
B30spkindf$conc<-log2(B30spkindf$conc)

#b31 mix 2
B31spkindf<-data.frame(ERCC$B31,ercc_conc$concentration.in.Mix.2..attomoles.ul.)
names(B31spkindf)<-c("FPKM", "conc")
B31spkindf$FPKM<-as.integer(B31spkindf$FPKM)
B31spkindf$conc<-as.integer(B31spkindf$conc)
row_sub = apply(B31spkindf, 1, function(row) all(row !=0 ))
B31spkindf=B31spkindf[row_sub,]
B31spkindf$FPKM<-log2(B31spkindf$FPKM)
B31spkindf$conc<-log2(B31spkindf$conc)

#b32 mix 1
B32spkindf<-data.frame(ERCC$B32,ercc_conc$concentration.in.Mix.1..attomoles.ul.)
names(B32spkindf)<-c("FPKM", "conc")
B32spkindf$FPKM<-as.integer(B32spkindf$FPKM)
B32spkindf$conc<-as.integer(B32spkindf$conc)
row_sub = apply(B32spkindf, 1, function(row) all(row !=0 ))
B32spkindf=B32spkindf[row_sub,]
B32spkindf$FPKM<-log2(B32spkindf$FPKM)
B32spkindf$conc<-log2(B32spkindf$conc)

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
  theme(axis.line = element_line(color = 'black'))+ geom_text(x = 4.5, y = 15, label = lm_eqn(LO2spkindf), parse = TRUE) +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO3
spikeinPlotLO3<-ggplot(LO3spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO3/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO3spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

# LO4
spikeinPlotLO4<-ggplot(LO4spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO4/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO4spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#LO5
spikeinPlotLO5<-ggplot(LO5spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("LO5/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 4.5, y = 15, label = lm_eqn(LO5spkindf), parse = TRUE) +
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B29
spikeinPlotB29<-ggplot(B29spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B29/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B29spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B30
spikeinPlotB30<-ggplot(B30spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B30/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B30spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B31
spikeinPlotB31<-ggplot(B31spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B31/spike-in 2")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B31spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

#B32
spikeinPlotB32<-ggplot(B32spkindf, aes(x=conc,y=FPKM)) +geom_point()+ggtitle("B32/spike-in 1")+
  theme(panel.background =element_blank(), panel.grid.major = element_blank(), panel.grid.minor =element_blank(), panel.border =element_blank()) +
  theme(axis.line = element_line(color = 'black')) +
  geom_text(x = 5, y = 11, label = lm_eqn(B32spkindf), parse = TRUE)+
  geom_smooth(method=lm, se=FALSE, color="black") +
  labs(x=expression(paste("Log"[2]," concentration", sep=""))) +
  labs(y=expression(paste("Log"[2]," FPKM", sep="")))

# put together concentration v. FPKM plots
library("gridExtra", lib.loc="~/Desktop/R-3.2.2/library")
grid.arrange(spikeinPlotLO2, spikeinPlotLO3, spikeinPlotLO4, spikeinPlotLO5, ncol=2, nrow=2)
grid.arrange(spikeinPlotB29, spikeinPlotB30, spikeinPlotB31, spikeinPlotB32, ncol=2, nrow=2)