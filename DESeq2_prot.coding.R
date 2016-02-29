rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_YFPs/HtSeq_reverse/")
library("DESeq2")
library("edgeR") # u need this for the RPKM calculations

####################################################################################################################################
# Run DESeq2   --> Treated =YFP_POS
####################################################################################################################################


# Run only once to get tables, then add up tech reps. 

#files <- c(
#        "treatedB1T1.txt",
 #       "treatedB1T2.txt",
  #      "treatedB2T1.txt",
   #     "treatedB2T2.txt",
    #    "treatedB3T1.txt",
     #   "treatedB3T2.txt",
      #  "untreatedB1T1.txt",
       # "untreatedB1T2.txt",
        #"untreatedB2T1.txt",
        #"untreatedB2T2.txt",
        #"untreatedB3T1.txt",
        #"untreatedB3T2.txt",
        #"untreatedB4T1.txt",
        #"untreatedB4T2.txt"
#)

#i = 1
#while(i<14){
 #       file = read.table(files[i],sep="\t")
  #      tmp = read.table(files[i+1],sep="\t")
   #     treatB1 = merge(file,tmp,by=c(1))
    #    treatB1$count = treatB1[,2]+treatB1[,3]
     #   treatB1 = subset(treatB1,select=c(V1,count))
      #  write.table(treatB1,file=paste("treat",i,sep="."),sep="\t",quote=F,row.names=F,col.names=F)
       # i = i+2
#}

#generate protein coding gene tables

files <- c(
  "treatedB1.txt",
  "treatedB2.txt",
  "treatedB3.txt",
  "untreatedB1.txt",
  "untreatedB2.txt",
  "untreatedB3.txt",
  "untreatedB4.txt"
)

anno = read.table("TAIR10.gene.list.txt",sep="\t")

for (filename in files) {
  AllCountsFinal = read.table(filename,sep="\t")
  AllCountsFinal = merge(AllCountsFinal,anno,by=c(1))
  AllCountsFinal = AllCountsFinal[AllCountsFinal[,3]=="protein_coding_gene",]
  AllCountsFinal = AllCountsFinal[,-3]
  write.table(AllCountsFinal,file=paste(filename,".pro.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}

############### 
#Import files

sampleFiles <- grep("treated.*pro",list.files("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_YFPs/HtSeq_reverse/"),value=TRUE)
sampleFiles

sampleCondition<-c("treated","treated","treated",
                   "untreated","untreated","untreated","untreated")
sampleCondition

sampleTable<-data.frame(sampleName=sampleFiles,fileName=sampleFiles, condition=sampleCondition)
sampleTable

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory="~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_YFPs/HtSeq_reverse/", design=~condition)
ddsHTSeq

#Levels in colData are important bc they're used in log calculations; set untreated/control 1st
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
ddsHTSeq

#Standard differential expression analysis steps are wrapped into a single function, DESeq.
dds<-DESeq(ddsHTSeq)
str(dds) 
res<-results(dds)
res<-res[order(res$pvalue),] # order results table by the smallest adjusted p value:


#Information about variables and tests were used can be found by calling the function mcols on the results object.
mcols(res)$description
mcols(res,use.names=TRUE)

head(results(dds, addMLE=TRUE),10)

# Results which pass an adjusted p value threshold
resSig <- subset(res, pvalue  < 1) # genes with p-val SMALLER than 0.1
resSig 
dim(resSig)

#write.table(x=resSig, file="BL_YFPs_DESeq2_p-val_0.01")

#####################################################

###caculate rpkm for genes using edgeR###
count = counts(ddsHTSeq)
count = as.data.frame(count)
count$ID = rownames(count)

len = read.table("tair10.whole.genelength.txt",sep="\t")
names(len) = c("ID","length")

count = merge(count,len,by=c("ID"))
rpkm = rpkm(count[,2:8],count$length,normalized.lib.sizes=TRUE)
rpkm = as.data.frame(rpkm)
rpkm$ID = count$ID

####get rpkm of diff exp genes ####
df_exp = as.data.frame(resSig)
df_exp$ID = rownames(df_exp)
df_exp = merge(df_exp,rpkm,by=c("ID"))  ###final table contains rpkm of each trt
###order according to logfoldchange
df_exp = df_exp[order(df_exp[,3],decreasing = TRUE),]

###

write.table(df_exp, file="BL_Uni_YFPs_DESeq2_DEGs_TechRepsAdded_reversed")




