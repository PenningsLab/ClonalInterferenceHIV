
patstoexclude<-c() #filter out pats with less than 5 seqs, less than two time points or less than 2 seqs at the second most sample rich time point.
for (patname in List99Pats){
    filename=paste("OriginalData/FASTAfiles/",patname,".fasta",sep="")
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE)
    days<-sort(unique(substr(names(patfasta[,1]),5,7)),decreasing=TRUE)#last day first
    if (length(days)<2) patstoexclude<-c(patstoexclude, patname)
    else if (dim(patfasta)[1]<5) patstoexclude<-c(patstoexclude, patname)
    else if (sort(table(substr(names(patfasta[,1]),5,7)), decreasing = TRUE)[2]<2) {
        #print(patname)
        patstoexclude<-c(patstoexclude, patname)} #make sure the second most sample rich day has at least 2 sequences
} 
List99Pats<-List99Pats[-which(List99Pats%in%patstoexclude)]

#Count median number of days and sequences
numdays<-c(); numseqs<-c()
for (patname in List99Pats){
    filename=paste("OriginalData/FASTAfiles/",patname,".fasta",sep="")
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE)
    days<-sort(unique(substr(names(patfasta[,1]),5,7)),decreasing=TRUE)#last day first
    numdays<-c(numdays, length(days))
    numseqs<-c(numseqs,dim(patfasta)[1])
} 
