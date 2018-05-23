#Rscript to make MullerPlot
library(MullerPlot)
#prepare data and load some packages and functions
if (TRUE) source("Rscripts/PrepareForDataViz.R")


for (patname in c("P00026", 
"P00072",
"P00077",
"P00084",
"P00089",
"P00132",
"P00140")){
PatData<-PatientOverview[PatientOverview$patient==patname,MutColumns]

#OTUs: WT, 67+70, 70 , 41
if (patname=="P00026")listRelevantMuts<-c("NR41","NR67","NR70")
if (patname=="P00072")listRelevantMuts<-c("PI82","PI88","PI54", "PI46")
if (patname=="P00077")listRelevantMuts<-c("PI46","PI54","PI24")
if (patname=="P00084")listRelevantMuts<-c("NN103","NN188")
if (patname=="P00089")listRelevantMuts<-c("NN101", "NN103","NN190")
if (patname=="P00132")listRelevantMuts<-c("NN100", "NN101")
if (patname=="P00140")listRelevantMuts<-c("NR41", "NR67","NR70","NR215")

PatData<-PatData[,which(names(PatData)%in%listRelevantMuts)]
Timepoints<-PatientOverview$day[PatientOverview$patient==patname]

DataTable<-as.data.frame(t(PatData))
x<-nrow(DataTable)+1
DataTable[x,]<-rep(0,ncol(DataTable))
row.names(DataTable)[x]<-"WT"
names(DataTable)<-Timepoints
if (patname == "P00026") {DataTable[3,]<-DataTable[3,]-DataTable[2,]
 rownames(DataTable)[which(rownames(DataTable)=="NR67")]<-"NR67/70"
}
if (patname=="P00072") {
    rownames(DataTable)[1]<-"PI46/54/82"
    rownames(DataTable)[2]<-"PI54/82"
    DataTable[3,]<-DataTable[3,]-DataTable[2,]
    DataTable[2,]<-DataTable[2,]-DataTable[1,]
    DataTable<-DataTable[c(4,3,1,2,5),]
}
if (patname=="P00077") {
    rownames(DataTable)[1]<-"PI24/46"
    DataTable[5,]<-c(0,0,0,0.25,0,0,0,0) #add row for 46/54
    DataTable[2,]<-DataTable[2,]-DataTable[5,]
    rownames(DataTable)[5]<-"PI46/54"
    DataTable[3,]<-DataTable[3,]-DataTable[5,]
    DataTable[2,]<-DataTable[2,]-DataTable[1,]
    DataTable<-DataTable[c(5,3,2,1,4),]
}
if (patname=="P00084"){
    DataTable[4,]<-c(0,0,0,0,0,0,0.1666667,0,0) #add row for 46/54
    rownames(DataTable)[4]<-"NN103/188"
    DataTable[1,7]<-0
    DataTable[2,7]<-1-0.1666667
    DataTable<-DataTable[c(1,4,2,3),]
}
if (patname=="P00089"){
    DataTable[5,]<-rep(0,ncol(DataTable))
    DataTable[6,]<-rep(0,ncol(DataTable))
    rownames(DataTable)[5]<-"NN101/190"
    rownames(DataTable)[6]<-"NN103/190"
    DataTable<-DataTable[c(1,2,3,5,6,4),]
    DataTable[,7]<-c(0,2/6,0,2/6,1/6,1/6)
    DataTable[,8]<-c(0,4/7,0,2/7,1/7,0)
    DataTable[,9]<-c(0,5/7,0,2/7,0,0)
    DataTable[,10]<-c(1/6,1/6,0,4/6,0,0)
    DataTable[,11]<-c(0,1/7,0,6/7,0,0)
    DataTable<-DataTable[c(2,5,1,4,3,6),] #put 103 as the first one
}
if (patname =="P00132"){
    DataTable[4,]<-rep(0,ncol(DataTable))
    rownames(DataTable)[4]<-"NN100/101"
    DataTable[,5]<-c(5/8,2/8,0,1/8)
    DataTable<-DataTable[c(2,4,1,3),]
}
if(patname=="P00140"){
    DataTable[6,]<-rep(0,ncol(DataTable))
    rownames(DataTable)[6]<-"NR67/70"
    DataTable[7,]<-rep(0,ncol(DataTable))
    rownames(DataTable)[7]<-"N41/67/215"
    DataTable[8,]<-rep(0,ncol(DataTable))
    rownames(DataTable)[8]<-"N41/67/70/215"
    DataTable<-DataTable[c(1,2,3,4,6,7,8,5),]
    DataTable[,1]<-c(0,0,0,0,5/8,3/8,0,0)
    DataTable[,2]<-c(0,0,0,0,0,0,1,0)
    DataTable[,3]<-c(0,0,0,0,0,1,0,0) 
}

x=which(rownames(DataTable)=="WT")
for (i in 1:ncol(DataTable)){
    DataTable[x,i]<-max(0,1-sum(DataTable[1:(x-1),i]))
}

DataTable<-DataTable[which(apply(DataTable,MARGIN = 1,sum)!=0),]

SummedDataTable<-DataTable
for (i in (nrow(DataTable)-1):1){
    SummedDataTable[i,] = SummedDataTable[i,] + SummedDataTable[i+1,]
}

png(paste(patname,"Muller.png",sep=""),width=500,height=400,units="px",pointsize=12,bg="white")

layout(matrix(c(1,1)))
par(mar=c(3, 3, 2, 0))
plot(1,1,col=0,xlim=c(0,as.numeric(names(DataTable)[ncol(DataTable)])),ylim=c(0,1),main = paste("patient", substr(patname, 4,6)))
xvals<-as.numeric(names(DataTable))
xvals<-c(xvals,xvals[length(xvals):1])
for (i in 1:nrow(DataTable)){
    polygon(x=xvals,y=1-c(t(SummedDataTable[i,]),rep(0,ncol(SummedDataTable))),col=brewer.pal(7,"Set3")[i+1])
}
legend("right", legend = rownames(DataTable) ,col = brewer.pal(7,"Set3")[2:10],pch = 19,cex=1.3)

dev.off()
}
