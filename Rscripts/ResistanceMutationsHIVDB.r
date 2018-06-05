#Major HIV-1 Drug Resistance Mutations 
#Updated March 9, 2015 
#Updated June 2018 https://hivdb.stanford.edu/assets/media/resistance-mutation-handout-Dec2017.b8f72e32.pdf


NNRTImuts<-data.frame("pos"=numeric(8),"wt"=numeric(8),"mut"=numeric(8))
NNRTImuts$pos<-c(100,101,103,106,181,188,190,230)
NNRTImuts$mut<-c("I","PE","NS","MA","CIV","LCH","SAEQ","L")
NRTImuts<-data.frame("pos"=numeric(11),"wt"=numeric(11),"mut"=numeric(11))
NRTImuts$pos<-c(41,65,67,70,74,115,151,184,210,215,219)
NRTImuts$mut<-c("L","R","N","E","IV","F","M","IV","W","YF","QE")
PImuts<-data.frame("pos"=numeric(11),"wt"=numeric(11),"mut"=numeric(11))
PImuts$pos<-c(32,46,47,48,50,54,76,82,84,88,90)
PImuts$mut<-c("I","IL","VA","VM","LV","VTALM","V","AFTS","V","S","M")
NumImportantMuts=length(c(NNRTImuts,NRTImuts,PImuts))
PositionsRT<-sort(c(NNRTImuts$pos,NRTImuts$pos))
PositionsPRO<-PImuts$pos

AllNtPositionsInvolvedInResistance<-sort(c(c((PositionsRT+99),PositionsPRO)*3,c((PositionsRT+99),PositionsPRO)*3-1,c((PositionsRT+99),PositionsPRO)*3-2))

RTImuts<-rbind(NNRTImuts,NRTImuts)
RTImuts$mut[which(RTImuts$pos==67)]<-"NR"



