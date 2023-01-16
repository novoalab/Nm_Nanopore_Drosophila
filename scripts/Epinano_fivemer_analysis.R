library(plyr)
library(ggplot2)
library(ggrepel)
library(MASS)
library(reshape2)

pred <- read.table("Nm_Nanopore_Drosophila/Nm_pos.bed", sep="\t", header=F, stringsAsFactors = F)
pred$chr_pos <- paste(pred$V1, pred$V3, sep=":")
pred$V4 <- as.character(pred$V4)
pred <- pred[grep("[A/G/C/U/Ψ]m", pred$V4),]


#Cleanup
cleanup <- function(input, label) {
  #Filter low coverage reads
  input <- subset(input, cov>30)
  #Filter read starts
  input <- subset(input, pos>20)
  #For Stress
  #Add a column with position
  input$position<- paste(input$X.Ref,input$pos)
  #Change column names
  input <- input[, c("X.Ref","pos","position", feature)]
  colnames(input)<- c("Chr","Position","chr_pos",feature )
  data_melted<- melt(data = input, id.vars = c("Chr", "Position", "chr_pos"))
  colnames(data_melted)[which(names(data_melted) == "value")] <- paste(label, "value", sep="_")
  return(data_melted)
}

###ANALYSIS FOR REPLICATE 1 

#Arguments 
fivemer1 <- read.delim("NT.tsv.per.site.var.per_site_var.5mer.csv",sep=",")  #1st variable
label1 <- as.character("NT")  #1st label
fivemer2 <-read.delim("FBL.tsv.per.site.var.per_site_var.5mer.csv", sep=",") #2nd variable
label2 <- as.character("FBL") #2nd label
feature<- as.character("sum") #Feature

##fivemers1 -- getting the sum of the score for 5 consecutive positions

fivemer1 <- separate(fivemer1, col = Window, into = c("1", "2", "3", "4", "5"), sep = ":")
fivemer1 <- separate(fivemer1, col = Coverage, into = c("cov1", "cov2", "cov3", "cov4", "cov5"), sep = ":")
fivemer1$chr_pos <- paste(fivemer1$Ref, fivemer1$`3`, sep=":")

fivemer1$cov3 <- as.numeric(fivemer1$cov3)
fivemer1$`3` <- as.numeric(fivemer1$`3`)
kmer_rep1 <- fivemer1[,c("Ref", "3","cov3")]
kmer_rep1$mis <- fivemer1$mis1+fivemer1$mis2+fivemer1$mis3+fivemer1$mis4+fivemer1$mis5
kmer_rep1$del <- fivemer1$del1+fivemer1$del2+fivemer1$del3+fivemer1$del4+fivemer1$del5
kmer_rep1$ins <- fivemer1$ins1+fivemer1$ins2+fivemer1$ins3+fivemer1$ins4+fivemer1$ins5
kmer_rep1$sum <- kmer_rep1$del+kmer_rep1$ins+kmer_rep1$mis
input1 <- kmer_rep1
colnames(input1) <- c("X.Ref", "pos", "cov", "mis", "del", "ins", "sum")

##fivemers2 -- getting the sum of the score for 5 consecutive positions
fivemer2 <- separate(fivemer2, col = Window, into = c("1", "2", "3", "4", "5"), sep = ":")
fivemer2 <- separate(fivemer2, col = Coverage, into = c("cov1", "cov2", "cov3", "cov4", "cov5"), sep = ":")
fivemer2$chr_pos <- paste(fivemer2$Ref, fivemer2$`3`, sep=":")

fivemer2$cov3 <- as.numeric(fivemer2$cov3)
fivemer2$`3` <- as.numeric(fivemer2$`3`)
kmer_rep2 <- fivemer2[,c("Ref", "3","cov3")]
kmer_rep2$mis <- fivemer2$mis1+fivemer2$mis2+fivemer2$mis3+fivemer2$mis4+fivemer2$mis5
kmer_rep2$del <- fivemer2$del1+fivemer2$del2+fivemer2$del3+fivemer2$del4+fivemer2$del5
kmer_rep2$ins <- fivemer2$ins1+fivemer2$ins2+fivemer2$ins3+fivemer2$ins4+fivemer2$ins5
kmer_rep2$sum <- kmer_rep2$del+kmer_rep2$ins+kmer_rep2$mis
input2 <- kmer_rep2
colnames(input2) <- c("X.Ref", "pos", "cov", "mis", "del", "ins", "sum")

#Cleanup and process the data
data1 <- cleanup(input1, label1)
data2 <- cleanup(input2, label2)

data1$NT_value <- rescale(data1$NT_value, to = c(0, 1))
data2$FBL_value <- rescale(data2$FBL_value, to = c(0, 1))

merged <- merge(data1,data2[, c(3,5)], by="chr_pos")
merged$chr_pos <- paste(merged$Chr, merged$Position, sep = ":")
merged$Chr <- NULL
merged$Position <- NULL
merged$base <- NULL
merged$variable <- NULL


merged$score<- abs(merged[,c(paste(label1, "value", sep="_"))] - merged[,c(paste(label2, "value", sep="_"))])
merged <- separate(merged, col = chr_pos, into =c("chr", "pos"), sep=":", remove=F)
merged$pos <- as.numeric(merged$pos)
#merged <- merged[merged$pos > 50,] 

merged1 <- merged


###ANALYSIS FOR REPLICATE 2

#Arguments
fivemer1 <- read.delim("NT.tsv.per.site.var.per_site_var.5mer.csv",sep=",")  #1st variable
label1 <- as.character("NT")  #1st label
fivemer2 <-read.delim("FBL.tsv.per.site.var.per_site_var.5mer.csv", sep=",") #2nd variable
label2 <- as.character("FBL") #2nd label
feature<- as.character("sum") #Feature

##fivemers1 -- getting the sum of the score for 5 consecutive positions

fivemer1 <- separate(fivemer1, col = Window, into = c("1", "2", "3", "4", "5"), sep = ":")
fivemer1 <- separate(fivemer1, col = Coverage, into = c("cov1", "cov2", "cov3", "cov4", "cov5"), sep = ":")
fivemer1$chr_pos <- paste(fivemer1$Ref, fivemer1$`3`, sep=":")

fivemer1$cov3 <- as.numeric(fivemer1$cov3)
fivemer1$`3` <- as.numeric(fivemer1$`3`)
kmer_rep1 <- fivemer1[,c("Ref", "3","cov3")]
kmer_rep1$mis <- fivemer1$mis1+fivemer1$mis2+fivemer1$mis3+fivemer1$mis4+fivemer1$mis5
kmer_rep1$del <- fivemer1$del1+fivemer1$del2+fivemer1$del3+fivemer1$del4+fivemer1$del5
kmer_rep1$ins <- fivemer1$ins1+fivemer1$ins2+fivemer1$ins3+fivemer1$ins4+fivemer1$ins5
kmer_rep1$sum <- kmer_rep1$del+kmer_rep1$ins+kmer_rep1$mis
input1 <- kmer_rep1
colnames(input1) <- c("X.Ref", "pos", "cov", "mis", "del", "ins", "sum")



##fivemers2  -- getting the sum of the score for 5 consecutive positions

fivemer2 <- separate(fivemer2, col = Window, into = c("1", "2", "3", "4", "5"), sep = ":")
fivemer2 <- separate(fivemer2, col = Coverage, into = c("cov1", "cov2", "cov3", "cov4", "cov5"), sep = ":")
fivemer2$chr_pos <- paste(fivemer2$Ref, fivemer2$`3`, sep=":")

fivemer2$cov3 <- as.numeric(fivemer2$cov3)
fivemer2$`3` <- as.numeric(fivemer2$`3`)
kmer_rep2 <- fivemer2[,c("Ref", "3","cov3")]
kmer_rep2$mis <- fivemer2$mis1+fivemer2$mis2+fivemer2$mis3+fivemer2$mis4+fivemer2$mis5
kmer_rep2$del <- fivemer2$del1+fivemer2$del2+fivemer2$del3+fivemer2$del4+fivemer2$del5
kmer_rep2$ins <- fivemer2$ins1+fivemer2$ins2+fivemer2$ins3+fivemer2$ins4+fivemer2$ins5
kmer_rep2$sum <- kmer_rep2$del+kmer_rep2$ins+kmer_rep2$mis
input2 <- kmer_rep2
colnames(input2) <- c("X.Ref", "pos", "cov", "mis", "del", "ins", "sum")



#Cleanup and process the data
data1 <- cleanup(input1, label1)
data2 <- cleanup(input2, label2)

merged <- merge(data1,data2[, c(3,5)], by="chr_pos")
merged$chr_pos <- paste(merged$Chr, merged$Position, sep = ":")
merged$Chr <- NULL
merged$Position <- NULL
merged$base <- NULL
merged$variable <- NULL


merged$score<- abs(merged[,c(paste(label1, "value", sep="_"))] - merged[,c(paste(label2, "value", sep="_"))])



merged <- separate(merged, col = chr_pos, into =c("chr", "pos"), sep=":", remove=F)
merged$pos <- as.numeric(merged$pos)
#merged <- merged[merged$pos > 50,]


##barplot
######
barplot<- function(subs){
  subs$score<- abs(subs[,c(paste(label1, "value", sep="_"))] - subs[,c(paste(label2, "value", sep="_"))])
  sites <- pred[pred$chr_pos %in% subs[subs$score > 5*median(subs$score) & subs$score>0.03,]$chr_pos,]
  pos1 <- position_jitter(width = 0, seed = 0)
  redbars <- subs[subs$chr_pos %in% sites[sites$V1 == ref,]$chr_pos,]
  #pdf(file=paste("18S","sum", label1, label2, "barplot.pdf", sep="_"),height=5,width=15,onefile=FALSE)
  print(ggplot(subs, aes(x=pos, y=score)) +
          geom_bar(stat = "identity", width=1, fill="deepskyblue4") +
          geom_bar(data=subs[subs$score > 3*median(subs$score),], aes(x=pos, y=score),stat = "identity", width=1, fill="tan1")+
          geom_bar(data=subs[subs$score > 5*median(subs$score),], aes(x=pos, y=score),stat = "identity", width=1, fill="tomato")+
          ggtitle(paste("summed_errors_per_5mer", label1, label2, sep="_"))+
          geom_vline(xintercept = pred[pred$V1 == ref,]$V3, linetype = "dashed", color = "slategray3")+
          #geom_label_repel(data=redbars, aes(x=pos, y=max(redbars$score)-0.2, label= pos),size=5, color="black", segment.size  = 0.2,segment.color = "black", fill="white", position = pos1, min.segment.length = 0)+
          geom_label_repel(data=redbars, aes(x=pos, y=0.09, label= pos),size=6, color="black", segment.size  = 0.2,segment.color = "black", fill="white", min.segment.length = 0)+
          #geom_text_repel(data=subset(subs, score > 4*median(subs$score)), aes(pos, score, label=pos),size=3, color="red", segment.size  = 1,segment.color = "black")+
          #geom_hline(yintercept = median(subs$score), linetype = "dashed", color= "dimgray")+
          xlab("Position")+
          ylab("Delta summed errors") +
          ylim(0, 0.1)+ ## here specify a subset of coordinates of the transcript if you need to
          theme_bw()+
          theme(axis.text.x = element_text(face="bold", color="black",size=20),
                axis.text.y = element_text(face="bold", color="black", size=20),
                plot.title = element_text(color="black", size=24, face="bold",hjust = 0.5),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                panel.background = element_blank(),
                legend.position = "none",
                axis.line = element_line(colour = "black", size=0.5)))
  #dev.off()
}

for (ref in levels(merged$chr)){
  pdf(file=paste(ref,"sum", "barplot_rep2.pdf", sep="_"),height=4,width=12,onefile=FALSE) 
  barplot(merged[merged$chr==ref,])
  dev.off()
  pdf(file=paste(ref,"sum", "barplot_rep1.pdf", sep="_"),height=4,width=12,onefile=FALSE) 
  barplot(merged1[merged1$chr==ref,])
  dev.off()
}

### final scores and sites in common between the two reps (4th and 5th column containind the score in rep1 and rep2 respectively):

write.table(merged1, "rep1_FBL_kmer_scores.txt", sep ="\t", col.names=F, row.names=F, quote=F)
write.table(merged, "rep2_FBL_kmer_scores.txt", sep ="\t", col.names=F, row.names=F, quote=F)

subs1 <- merged[merged$chr == "18S",]
subs2 <- merged1[merged1$chr =="18S",]

common18 <- merge(subs1[subs1$score > 3* median(subs1$score),c(1,2,3,6)], subs2[subs2$score > 3* median(subs2$score),c(1,6)], by = "chr_pos", all =F, suffixes= c("_rep1", "_rep2"))

subs1 <- merged[merged$chr == "28S",]
subs2 <- merged1[merged1$chr =="28S",]

common28 <- merge(subs1[subs1$score > 3* median(subs1$score),c(1,2,3,6)], subs2[subs2$score > 3* median(subs2$score),c(1,6)], by = "chr_pos", all =F, suffixes= c("_rep1", "_rep2"))

write.table(common18, "FBL_kmer_replicable_positions_18S.txt", sep ="\t", col.names=F, row.names=F, quote=F)
write.table(common28, "FBL_kmer_replicable_positions_28S.txt", sep ="\t", col.names=F, row.names=F, quote=F)
