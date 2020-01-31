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

NamesMuts<-names(PatientOverview)[MutColumns]
x<-which(MutCounts==0)
MutCounts<-MutCounts[-x]
NamesMuts<-NamesMuts[-x]

NamesMuts<-c("L100I", "K101EP", "K103NS", "V106AMI", "Y181CIV", "Y188LCH", "G190ASE",
             "M41L",  "K65R" , "D67N" , "K70R" , "L74VI" , "M184VI", "L210W", "T215FY", "K219QE",
             "L24I",  "V32I" , "M46IL"  ,"I47VA" , "I50LV" , "I54VTALM","L76V", "V82ATFS" , "I84V",   "N88DS" , "L90M" )

png("Output/MutOccurrence.png") 
par(mar = c(7, 4, 2, 2) + 0.2)
barplot(MutCounts/length(List99Pats), main="Occurrence of resistance mutations in Bacheler dataset", 
        names.arg = NamesMuts, 
        ylim=c(0,1), 
        #xlab = "Mutations", 
        ylab="Fraction of patients with mut. observed", 
        yaxt="n",
        las=2, 
        col = 0)

abline(h=(1:5)/5)
brewer.pal(8,"Set3")[3]

barplot(MutCounts/length(List99Pats), add=TRUE, las=1, col=c(rep(brewer.pal(6,"Set3")[3], 7), rep(brewer.pal(6,"Set3")[4], 9), rep(brewer.pal(6,"Set3")[2], 11)))
dev.off()
        