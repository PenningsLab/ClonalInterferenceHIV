#Make Barplot

# for each mutation , I want to count in how many of the 118 patients it is seen. 

#prepare data and load some packages and functions
if (TRUE) source("Rscripts/PrepareForDataViz.R")
if (TRUE) source("Rscripts/FilteringPatients.R")


MutCounts<-c()
for (i in 1:length(MutColumns)){
    Mut<-names(PatientOverview)[MutColumns[i]]
    count=0
    for (p in 1:length(List99Pats)){
        Pat<-List99Pats[p]
        if (max(PatientOverview[which(PatientOverview$patient==Pat),MutColumns[i]])>0) count = count+1
    }
    print(paste(Mut, count))
    MutCounts<-c(MutCounts,count)
}

png("Output/MutOccurrence.png") 
par(mar = c(7, 4, 2, 2) + 0.2)
barplot(MutCounts/length(List99Pats), main="Occurrence of resistance mutations in Bacheler dataset", 
        names.arg = names(PatientOverview)[MutColumns], 
        ylim=c(0,1), 
        #xlab = "Mutations", 
        ylab="Fraction of patients with mut. observed", 
        yaxt="n",
        las=2, 
        col = 0)

abline(h=(1:5)/5)
brewer.pal(8,"Set3")[3]

barplot(MutCounts/length(List99Pats), add=TRUE, las=1, col=c(rep(brewer.pal(6,"Set3")[3], 7), rep(brewer.pal(6,"Set3")[4], 11), rep(brewer.pal(6,"Set3")[2], 12)))
dev.off()
        