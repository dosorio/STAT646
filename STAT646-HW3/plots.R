
##plots---------

###plot1-----------------------------
breaches<-factor(pheno$X..career..breaches)
breaches = factor(breaches,levels(breaches)[c(1:3,6,4,5,7)])
levels(breaches)<-c("0" ,"1-9", "10-39", "40-99", "100-199", "200+", "200+")
pdf("/hpc/users/karmam01/new_plots/breaches.pdf",width= 8,height=4)
barplot(table(breaches), xlab="categories", ylab = "# breaches",main ="# breaches bar plot")
dev.off()

#----------------------------#

pdf("/hpc/users/karmam01/new_plots/detection-p-values.pdf",width= 8,height=4)
barplot(apply(detP,2,mean), main="mean detection p-values", col=as.numeric(factor(targets$group)), xaxt="none")
abline(h=0.01,col="red")
dev.off()

#----------------------------#
pdf("/hpc/users/karmam01/new_plots/piechart.pdf",width= 8,height=4)
slices<-c(2100/32, 1100/32)
lbls <- slices
lbls <- paste(lbls,"%",sep="")
pie(slices,labels = lbls, col=rainbow(length(lbls)),main="Pie Chart", cex=0.8)
legend("topright", c("< 39"," > 39"), cex=0.8, fill= rainbow(length(lbls)))
dev.off()

#----------------------------#
#----------------------------#

pdf("/hpc/users/karmam01/new_plots/dividedbarplot.pdf",width= 8,height=4)
counts <-table(pheno$hx.of.TBI,breaches)
barplot(counts,  xlab ="# of breachers", col=c("red","black"), legend=c("TBI=no", "TBI=yes"), ylim=c(0, 10))
dev.off()

##-------------------------------
snps <- getSnpBeta(rgset); 
colnames(snps); 
match(colnames(snps),paste0(targets$Slide,"_",targets$Array)); 
colnames(snps)<-targets[,1]; colnames(snps); 
x<-dist(t(snps), method = "manhattan");
hc<-hclust(x, method = "single");

pdf(file ="/sc/orga/projects/haghif01a/Mou/Cyril_moumita/dendo.pdf", width= 12,height=6)
par(cex=1,font=3)
plot(hc, main="Dendrogram ")
dev.off()

##
pdf(file = "/hpc/users/karmam01/new_plots/swan.pdf", width=11,height=6)
par(mfrow=c(1,2), cex= 0.75)
plotBetasByType(mset_baseline[,1], main = "Raw")
plotBetasByType(msetsw_baseline[,1], main = "SWAN")
dev.off()

##----------------------------------------------------------------------
pdf(file = "/hpc/users/karmam01/new_plots/hist_breacher.pdf", width=11,height=8.5)
plot.new()
plot.window(ylim = c(-2000, 2000), xlim = range(c(fcbreacherplus,fcbreacherminus)))
p1 <- list(axes = FALSE, xlab = "methylation M-value difference", ylab = "# of CpGs", main = "Breach: low vs high",col = "red", border ="black")
p2 <- list(axes = FALSE, xlab = "methylation M-value difference", ylab = "# of CpGs", main = "Breach: low vs high",col = "blue", border = "black")
par(new=TRUE)
do.call(hist, c(list(x = fcbreacherplus, ylim = c(-2000, 2000)), p1,breaks=25))
par(new=TRUE)
do.call(hist, c(list(x = fcbreacherminus, ylim = c(2000, -2000)), p2,breaks=25))
axis(side = 2, 
      at = pretty(par()$usr[3:4]), 
     labels = abs(pretty(par()$usr[3:4])))
axis(side = 1) 
dev.off()

##----------------------------------------------------------------------
pdf(file = "/hpc/users/karmam01/new_plots/hist_TBI.pdf", width=11,height=8.5)
plot.new()
plot.window(ylim = c(-2000, 2000), xlim = range(c(fcTBIplus,fcTBIminus)))
p1 <- list(axes = FALSE, xlab = "methylation M-value difference", ylab = "# of CpGs", main = "TBI: no vs yes",col = "red", border ="black")
p2 <- list(axes = FALSE, xlab = "methylation M-value difference", ylab = "# of CoGs", main = "TBI: no vs yes",col = "blue", border = "black")
par(new=TRUE)
do.call(hist, c(list(x = fcTBIplus, ylim = c(-2000, 2000)), p1,breaks=25))
par(new=TRUE)
do.call(hist, c(list(x = fcTBIminus, ylim = c(2000, -2000)), p2,breaks=25))
axis(side = 2, 
      at = pretty(par()$usr[3:4]), 
     labels = abs(pretty(par()$usr[3:4])))
axis(side = 1) 
dev.off()
##-------promoter region plot--------------------------------------------------------------------
positioninfo<-annotation$UCSC_RefGene_Group
probename<-annotation$Name
selected_probe<- which(probename  %in%selected_cpgsitebreacher )
position_selected<-positioninfo[selected_probe]
countTSS=sum(grepl("TSS",position_selected))
count5UTR=sum(grepl("5'UTR",position_selected))
countbody=sum(grepl("Body",position_selected))
countexon=sum(grepl("1stExon",position_selected))
cate = c("TSS","body","1stExon","5'UTR");
cnt=c(countTSS,countbody,countexon,count5UTR)
names(cnt)=cate
pdf(file = "/hpc/users/karmam01/new_plots/positionplot.pdf", width=11,height=8.5)
barplot(cnt, xlab="categories of region", ylab="",main="Barplot for regions") 
dev.off()
> cnt
    TSS    body 1stExon   5UTR 
   3672    5444     842    1710
   
##--------------Island info plot-------------------------------------------------------
Islandinfo<-annotation$Relation_to_Island
Islandinfo_selected<-Islandinfo[selected_probe]
countIsland=sum(grepl("Island",Islandinfo_selected))
countNShelf=sum(grepl("N_Shelf",Islandinfo_selected))
countNShore=sum(grepl("N_Shore",Islandinfo_selected))
countSShore=sum(grepl("S_Shore",Islandinfo_selected))
countSShelf=sum(grepl("S_Shelf",Islandinfo_selected))
countSOpenSea=sum(grepl("OpenSea",Islandinfo_selected))
cate_island = c("Island","N_Shelf","N_Shore","OpenSea","S_Shelf", "S_Shore");
cnt_island=c(countIsland,countNShelf,countNShore,countSOpenSea,countSShore,countSShelf)
names(cnt_island)=cate_island
pdf(file = "/hpc/users/karmam01/new_plots/Islandinfoplot.pdf", width=11,height=8.5)
barplot(cnt_island, xlab="categories of Island Name", ylab="",main="Barplot for Island Name")
dev.off()
> cnt_island
 Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
   3316     821    1917    5753    1493     775 
   

##-------------------------------beta value box plot low Vs high breaches------------------------
beta<-beta_baseline
mean_beta<-apply(beta,2,mean)
pdf(file = "/hpc/users/karmam01/new_plots/boxplot_mean.pdf", width=11,height=8.5)
boxplot(mean_beta~career, col=c("gold","darkgreen"),xlab="", ylab="beta values")
dev.off()

##annotation plots for low vs high----------------------------------------------------------
GRset_high<-GRset_baseline[,-as.vector(which(career=="low"))]
GRset_low<-GRset_baseline[,as.vector(which(career=="low"))]
annotation_low<-getAnnotation(GRset_low);
annotation_high<-getAnnotation(GRset_high);
positioninfo_low<-annotation_low$UCSC_RefGene_Group
positioninfo_high<-annotation_high$UCSC_RefGene_Group
probename_low<-annotation_low$Name
probename_high<-annotation_high$Name
selected_probelow<-which(probename_low  %in%selected_cpgsitebreacher )
selected_probehigh<-which(probename_high  %in%selected_cpgsitebreacher )
Lposition_selected<-positioninfo[selected_probelow]
Hposition_selected<-positioninfo[selected_probehigh]
LcountTSS=sum(grepl("TSS",position_selected))
Lcount5UTR=sum(grepl("5'UTR",Lposition_selected))
Lcountbody=sum(grepl("Body",Lposition_selected))
Lcountexon=sum(grepl("1stExon",Lposition_selected))
cate = c("TSS","body","1stExon","5'UTR");
cnt_low=c(LcountTSS,Lcountbody,Lcountexon,Lcount5UTR)
names(cnt_low)=cate
HcountTSS=sum(grepl("TSS",Hposition_selected))
Hcount5UTR=sum(grepl("5'UTR",Hposition_selected))
Hcountbody=sum(grepl("Body",Hposition_selected))
Hcountexon=sum(grepl("1stExon",Hposition_selected))
cate = c("TSS","body","1stExon","5'UTR");
cnt_high=c(HcountTSS,Hcountbody,Hcountexon,Hcount5UTR)
names(cnt_high)=cate
cnt_low
    TSS    body 1stExon   5'UTR 
   3672    5444     842    1710
cnt_high
    TSS    body 1stExon   5'UTR 
   3672    5444     842    1710 
   
 ##----------------------------
 pdf(file = "/hpc/users/karmam01/new_plots/volcano.pdf", width=11,height=8.5)
 par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
plot(top_breacher$logFC, -log10(top_breacher$P.Value), 
     xlab="", ylab="-log10 p-value")
 abline(a=-log10(0.05), b=0, col="red")
 dev.off()


#####-----sexplot----####
RSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)
beta_raw <- getBeta(RSet)
GRset <- mapToGenome(RSet)
xIndex <- which(seqnames(GRset) == "chrX")
yIndex <- which(seqnames(GRset) == "chrY")
xMed <- matrixStats::colMedians(beta_raw[xIndex,], na.rm=TRUE)
yMed <- matrixStats::colMedians(beta_raw[yIndex,], na.rm=TRUE)
x.m <- xMed
y.m <- yMed
require("car")
z<-recode(bL[,4], "c('F')='1'; else='0'")
zz<-as.factor(z)
par(mar=c(7,4,4,2))
plot(c(1,nrow(bL)),c(-1,1),type="n",axes=FALSE,xlab="",ylab="",
main="The median beta values for the XY chromosomes")
points(1:31,x.m,pch=20)
points(1:31,-y.m,pch=20)
axis(1,1:31, bL[,1] ,las=2,cex=0.6)
at<-c(0,0.2,0.4,0.6,0.8,1.0)
axis(2,at=c(at,-at[-1]),c(at,at[-1]))
abline(h=0)
mtext("chrX median beta",2,line=2,at=0.5,adj=0.5)
mtext("ChrY median beta",2,line=2,at=-0.5,adj=0.5)
box()
segments(c(1:48),-y.m,y1=x.m,col=c("forestgreen","red")[zz])

