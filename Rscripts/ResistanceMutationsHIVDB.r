#Major HIV-1 Drug Resistance Mutations 
#Updated March 9, 2015 
#Updated summary from the HIV Drug Resistance Database. This document can be downloaded from the http://hivdb.stanford.edu home page. Detailed and referenced versions of each drug class summary can be found at http://hivdb.stanford.edu/pages/drugSummaries.html 


NNRTImuts<-data.frame("pos"=numeric(10),"wt"=numeric(10),"mut"=numeric(10))
NNRTImuts$pos<-c(100,101,103,106,138,179,181,188,190,227)
NNRTImuts$mut<-c("I","PEH","NS","MA","KAGQ","DEF","CIV","LCH","SAEQ","LC")
NRTImuts<-data.frame("pos"=numeric(11),"wt"=numeric(11),"mut"=numeric(11))
NRTImuts$pos<-c(41,65,67,70,74,115,151,184,210,215,219)
NRTImuts$mut<-c("L","R","N","ER","IV","F","M","IV","W","YF","QE")
PImuts<-data.frame("pos"=numeric(12),"wt"=numeric(12),"mut"=numeric(12))
PImuts$pos<-c(24,32,46,47,48,50,54,76,82,84,88,90)
PImuts$mut<-c("I","I","IL","VA","VM","LV","VTALM","V","AFTS","V","S","M")
NumImportantMuts=10+11+12
PositionsRT<-sort(c(NNRTImuts$pos,NRTImuts$pos))
PositionsPRO<-PImuts$pos

AllNtPositionsInvolvedInResistance<-sort(c(c((PositionsRT+99),PositionsPRO)*3,c((PositionsRT+99),PositionsPRO)*3-1,c((PositionsRT+99),PositionsPRO)*3-2))

RTImuts<-rbind(NNRTImuts,NRTImuts)
RTImuts$mut[which(RTImuts$pos==67)]<-"NR"



