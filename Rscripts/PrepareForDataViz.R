if (TRUE){
    #  setwd("~/Dropbox/ExplDataAnalysisR/Kadie/ClonalInteference/RScripts")
    source("Rscripts/RfunctionsVisualizeData.r")	
    source("Rscripts/GetConsensusB.r")	
    read.csv("OriginalData/PatientOverview.csv")->PatientOverview  #data on 170 patients
    library(ape)
    library(seqinr)
    library(pegas)	
    library(RColorBrewer)
    NNRTIcolumns<-which(substr(names(PatientOverview),1,2)=="NN")
    NRTIcolumns<-which(substr(names(PatientOverview),1,2)=="NR")
    PIcolumns<-which(substr(names(PatientOverview),1,2)=="PI")
    MutColumns<-c(NNRTIcolumns,NRTIcolumns,PIcolumns)	
    #For EFV treatment, NN179, NN227 and NN138 are not relevant. I am going to remove them.
    C<-which(names(PatientOverview)=="NN138"|names(PatientOverview)=="NN179"|names(PatientOverview)=="NN227")
    MutColumns<-MutColumns[-which(MutColumns==12|MutColumns==13|MutColumns==17)]
    List99Pats<-unique(PatientOverview$patient)
    
}
