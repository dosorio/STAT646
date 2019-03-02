##########-----breacher data analysis for pre-explosive samples using minfi--###################
################---Moumita Karmakar----#######################
########--Date 25th February---##################
######################################################################################


## Download the data and save it to a folde. set the working directory to the folder

### Read the phenotype data. You should have 72 samples and 13 columns. remove the samples mentioned in the pdf.column named `sample' contains information regarding whether the sample is pre-explosive(denoted by 1)
# or post-explosive(denoted by 10). Subset the data for pre-explosion. You should have 32 sample for pre-explosion.

options(stringsAsFactors = FALSE)

#### Cleaning phenotype data ######
X<-read.csv("Phenotype_Breacher_HM450.csv",sep=",");rownames(X)<-X[,1];X<-X[,-1]
X<-X[-c(61,62,63,64,68,69,70,71,72),]
## collecting samples only for pre-explosion
baseline<-X[X$sample==1,]
### removing samples "3X01_1", "3X02_2"
baseline<-baseline[-c(26,32),]

#####  Need to correct the format of the categories for column named "X..career..breaches". 
##After correcting the column create the bar plot. re-categorize the column as "<39"  and ">39" 
replace(baseline$X..career..breaches, baseline$X..career..breaches=="1/9/2016", "1-9")
require("car")
career<-recode(baseline$X..career..breaches, "c('40-99','100-199','200-399','400+')='>39'; else='<39'")

### basic bar plot for different categories of career breaches---##
breaches<-factor(baseline$X..career..breaches)
breaches = factor(breaches,levels(breaches)[c(1:3,6,4,5,7)])
levels(breaches)<-c("0" ,"1-9", "10-39", "40-99", "100-199", "200+", "200+")
pdf("breaches.pdf",width= 8,height=4)
barplot(table(breaches), xlab="categories", ylab = "# breaches",main ="# breaches bar plot")
dev.off()
# --- pie chart after recategorizing the column career breaches -------#
pdf("piechart.pdf",width= 8,height=4)
slices<-c(2100/30, 900/30)
lbls <- slices
lbls <- paste(lbls,"%",sep="")
pie(slices,labels = lbls, col=rainbow(length(lbls)),main="Pie Chart", cex=0.8)
legend("topright", c("< 39"," > 39"), cex=0.8, fill= rainbow(length(lbls)))
dev.off()

#  divided bar plot to see the interaction between career breaches and history of TBI-----------------#
pdf("dividedbarplot.pdf",width= 8,height=4)
counts <-table(baseline$hx.of.TBI,breaches)
barplot(counts,  xlab ="# of breachers", col=c("red","black"), legend=c("TBI=no", "TBI=yes"), ylim=c(0, 10))
dev.off()

#####library needed to work with the rgset channel data
library(minfi)
library(limma)

#### reading the rgset channel data.

load("rgset_breacher.Rdata")
rgset<-valid_RGset

# removing the samples and subset the data for pre-explosion. end of this step youshould have 30 samples and all of them are from pre-explosion category.
rgset_reduced<-rgset[,-c(61,62,63,64,68,69,70,71,72)];
rgset_baseline<-rgset_reduced[,-which(c(seq(1,63)) %% 2 == 0)]

#### creating the pheotype data from the rgset data
pheno<-pData(rgset_baseline); dim(pheno);
# changing the column names of rgset channel data to sample names of phenotype
colnames(rgset_baseline)<-pheno$Sample_Name

## preprocessing and normalization. preprocessRaw doesn't perform any kind of normalization. preprocessSWAN
# does Subset-quantile within Array Normalization

mset_baseline<-preprocessRaw(rgset_baseline);
msetillumina_baseline<-preprocessIllumina(rgset_baseline, bg.correct = TRUE, normalize = "controls")
msetsw_baseline<-preprocessSWAN(rgset_baseline, msetillumina_baseline)

####filter out poor quality probes
detP <- detectionP(rgset_baseline)
keep <- rowSums(detP < 0.01) == ncol(rgset_baseline)
mSetSw_baseline <- msetsw_baseline[keep,]

##extracting beta and M values -----------------------
beta_baseline <- getBeta(mSetSw_baseline)
Mval_baseline <-getM(mSetSw_baseline)

###differential methylation using limma for pre-explosion samples only-------------------------
###first matching the order of the samples between file named baseline and file named pheno
baseline<-baseline[order(rownames(baseline),pheno$Sample_Name),]

require("car")
career_new<-recode(baseline$X..career..breaches, "c('40-99','100-199','200-399','400+')='high'; else='low'")

career_breacher <- factor(career_new,levels=c("low","high"))
historyTBI <- factor(baseline$hx.of.TBI,levels=c("yes","no"))
design <- model.matrix(~pheno$age + career_breacher + historyTBI)
fit<- lmFit(Mval_baseline,design)
fit.reduced <- eBayes(fit)

####### top table of coefficient  for high vs low breaching
top_breacher<-topTable(fit.reduced,coef=3,number=dim(Mval_baseline)[1])
fc_breacher<-top_breacher$logFC
fcbreacherplus<-fc_breacher[fc_breacher>0]
fcbreacherminus<-abs(fc_breacher[fc_breacher<=0])

####### top table of coefficient  for history of TBI
top_TBI<-topTable(fit.reduced,coef=4,number=dim(Mval_baseline)[1])
fc_TBI<-top_TBI$logFC
fcTBIplus<-fc_TBI[fc_TBI>0]
fcTBIminus<-abs(fc_TBI[fc_TBI<=0])

### Accessing annotation for Illumina methylation objects
annotation<-getAnnotation(rgset_baseline)


###follow the plots.R for rest of the plots





