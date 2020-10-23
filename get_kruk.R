library(ggplot2)
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################
filename="lipid_data"

my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)	
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=log10(A[1:d[1], 5:d[2]])
rownames(B)=A[,1]
Cname=colnames(A)[5:d[2]]
#colnames(A)[1:4]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)
ID=rep(NA,Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    ID[mm]=splitname[[mm]][4]
}
Group=c(Group1,Group2, Group3)
GSSID=c(ID,ID,ID)
nn=sort(unique(Group))
nn<-nn[nn!='NotUse']
mmm=paste0(nn, collapse=',')
print(paste("probe","syn", "GeneName", "EZID","KP",
            "wilcox.Adult_Baby","wilcox.20_60", "wilcox.FT_LPT","wilcox.FT3_PT3", "wilcox.FT3_PT5", "wilcox.PT3_PT5", 
            "truefc.Adult_Baby","truefc.20_60", "truefc.FT_LPT","truefc.FT3_PT3", "truefc.FT3_PT5", "truefc.PT3_PT5",
            "AverageSignal",mmm, sep=","))

for (i in 1:d[1]){
    probe=A[i,1]
    syn=A[i,2]
    long=A[i,3]
    ez=A[i,4]
    gene0<-as.numeric(A[i, 5:d[2]])
    gene<-as.numeric(B[i,])
    log10Gene=c(gene, gene, gene)
    Gene=c(gene0,gene0,gene0)
    mydata0=data.frame(log10Gene,Gene,Group, GSSID )
    mydata<-mydata0[Group!="NotUse",]
    KP=kruskal.test(Gene ~ Group, data = mydata)$p.value
    pAB=my.wilcox.p.value(Gene[Group=="Adult"], Gene[Group=="Baby"])
    pOY=my.wilcox.p.value(Gene[Group=="A20s"], Gene[Group=="A60s"])
    pFL=my.wilcox.p.value(Gene[Group=="FT"], Gene[Group=="LPT"])
    pF3P3=my.wilcox.p.value(Gene[Group=="FT3"], Gene[Group=="PT3"])
    pF3P5=my.wilcox.p.value(Gene[Group=="FT3"], Gene[Group=="PT5"])
    pP3P5=my.wilcox.p.value(Gene[Group=="PT3"], Gene[Group=="PT5"])
    
    fAB=truefc(mean(Gene[Group=="Adult"])/mean(Gene[Group=="Baby"]))
    fOY=truefc(mean(Gene[Group=="A20s"])/mean(Gene[Group=="A60s"]))
    fFL=truefc(mean(Gene[Group=="FT"])/mean(Gene[Group=="LPT"]))
    fF3P3=truefc(mean(Gene[Group=="FT3"])/mean(Gene[Group=="PT3"]))
    fF3P5=truefc(mean(Gene[Group=="FT3"])/mean(Gene[Group=="PT5"]))
    fP3P5=truefc(mean(Gene[Group=="PT3"])/mean(Gene[Group=="PT5"]))
    
    #mean=aggregate(x = mydata$Gene,by = list(mydata$Group),function(x) c(mean=mean(x), sd=sd(x)))
    MMM=aggregate(x = mydata$Gene,by = list(mydata$Group),function(x) (mean=mean(x)))
    if(KP<=0.05){
    jpeg(paste0(syn,".",probe,".jpg"), width = 300, height = 300)
    p1<-ggplot(mydata, aes(x=Group, y=Gene)) + 
      geom_boxplot()+
      labs(title=paste(probe, syn, long) ,x="Group", y = "Expression")+
      theme_classic()+ geom_jitter(shape=16, position=position_jitter(0.2))
    print(p1)
    dev.off()
    }
    
    XXX=paste(MMM$x, collapse  =",")
    print(paste(probe,syn, long, ez,KP,pAB, pOY, pFL, pF3P3, pF3P5, pP3P5, 
                fAB, fOY, fFL, fF3P3, fF3P5, fP3P5,  
                mean(Gene), XXX, sep=","))
    
} 
