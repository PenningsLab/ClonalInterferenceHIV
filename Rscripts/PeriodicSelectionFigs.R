#prepare data and load some packages and functions
if (TRUE) source("Rscripts/PrepareForDataViz.R")

    for (patname in  c("P00094", "P00021", "P00077")){
    PatData<-PatientOverview[PatientOverview$patient==patname,MutColumns]
    
    Timepoints<-PatientOverview$day[PatientOverview$patient==patname]
    
    DataTable<-as.data.frame(t(PatData))
    names(DataTable)<-Timepoints
    
    x<-apply(DataTable, MARGIN = 2, sum)
    
    png(paste("Output/", patname,"PeriodicSelection.png",sep=""),width=600,height=400,units="px",pointsize=12,bg="white")
    
    plot(names(x), x, pch = 21, cex=3, bg = brewer.pal(7,"Set3"),
         ylab="Mean number of res mutations", 
         main = "Average number of resistance mutations per sequence",
         ylim=c(0,max(x)+1), 
         xlab="days in trial"
         )
    mtext(text=paste("Viral sequences from patient",substr(patname,4,6)), side=3)
    abline(h=1:5, lty=2)
    points(names(x), x, pch = 21, cex=3, bg = brewer.pal(7,"Set3"))
    
    dev.off()
}
