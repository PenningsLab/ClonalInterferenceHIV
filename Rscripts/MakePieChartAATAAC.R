
#prepare data and load some packages and functions
if (TRUE) source("Rscripts/PrepareForDataViz.R")

Counting606 <- c(0,0,0,0)

for (patname in List99Pats){
    #set filename and read fasta file into patfasta
    filename=paste("OriginalData/FASTAfiles/",patname,".fasta",sep="")
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE)
    
    #what is there at pos 606?
    numA<-length(which(patfasta[,606]=="a"))
    numT<-length(which(patfasta[,606]=="t"))
    numC<-length(which(patfasta[,606]=="c"))

    if (numC==0&numT==0) Counting606[1]=Counting606[1]+1 #neither C nor T
    if (numC==0&numT>0) Counting606[2]=Counting606[2]+1 #T but not C
    if (numC>0&numT==0) Counting606[3]=Counting606[3]+1 #C but not T
    if (numC>0&numT>0) Counting606[4]=Counting606[4]+1 #T and C
}

png("AACAATPiechart.png",width=500,height=400,units="px",pointsize=12,bg="white")
par(mar=c(1,1,2,1))
pie(Counting606,labels = c("neither AAC nor AAT", "only AAT", "only AAC", "both AAC and AAT"), col = brewer.pal(8,"Set3")[c(3,6,4,5)], main = "Occurrence of AAC and AAT codons at position RT 103")    
dev.off()
    

Counting41 <- c(0,0,0,0) #RT M41L mutation at nucleotide 418

for (patname in List99Pats){
    #set filename and read fasta file into patfasta
    filename=paste("OriginalData/FASTAfiles/",patname,".fasta",sep="")
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE)
    
    #what is there at pos 606?
    numA<-length(which(patfasta[,418]=="a"))
    numT<-length(which(patfasta[,418]=="t"))
    numC<-length(which(patfasta[,418]=="c"))
    
    if (numC==0&numT==0) Counting41[1]=Counting41[1]+1 #neither C nor T
    if (numC>0&numT==0) Counting41[2]=Counting41[2]+1 #C but not T
    if (numC==0&numT>0) Counting41[3]=Counting41[3]+1 #T but not C
    if (numC>0&numT>0) Counting41[4]=Counting41[4]+1 #T and C
}

png("CTCTTC_M41L_Piechart.png",width=500,height=400,units="px",pointsize=12,bg="white")
par(mar=c(1,1,2,1))
pie(Counting41,labels = c("neither CTC nor TTC", "only CTC", "only TTC", "both CTC and TTC"), col = brewer.pal(8,"Set3")[3:8], main = "Occurrence of CTC and TTC codons at position RT 41")    
dev.off()

