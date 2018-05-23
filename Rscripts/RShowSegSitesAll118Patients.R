#prepare data and load some packages and functions
if (TRUE) source("Rscripts/PrepareForDataViz.R")

shorterpat13=FALSE
ShowSingletons=FALSE
ShowOnlyResSites=TRUE

for (patname in List99Pats[89]){
    #gather information about the sequences in this patient before the sweep and at the timepoint where the sweep is detected. 
    if (TRUE){	
        #set filename and read fasta file into patfasta
        filename=paste("OriginalData/FASTAfiles/",patname,".fasta",sep="")
        patfasta<-read.dna(filename, format = "fasta",as.character=TRUE)
        patfastafasta<-read.dna(filename, format = "fasta")
        #seqlabels is the list of sequence labels
        seqlabels<-sort(names(patfasta[,1]),decreasing=TRUE)#last day first
        
        days<-sort(unique(substr(names(patfasta[,1]),5,7)),decreasing=TRUE)#last day first
        alldays<-sort(substr(names(patfasta[,1]),5,7),decreasing=TRUE)#last day first	
        if (shorterpat13) {days <-days[3:4] ; alldays <- alldays [11:24]; seqlabels<-seqlabels[11:24]}
        #which of the sequences are from the first day of treatment
        DayTipsDay0<-which(substr(names(patfasta[,1]),5,7)==days[length(days)])
        #get.ListOfSegSites is a function I wrote to get information about each site in the sequence (e.g. whether it is polymorphic, silent etc.) 
        AAseqs<-read.table("OriginalData/AAseqs.csv",header=TRUE,sep=",")
        L<-get.ListOfSegSites(patfasta,DayTipsDay0,AAseqs)
        numseqs=length(seqlabels)
    }	
    
    #prepare for making a figure (esp which sites to include in the figure)
    if (TRUE){
        #how many sites to draw	
        #find which resistance sites to include in the figure  	
        PatSpecResistanceCodons<-vector()
        for (i in 1:length(MutColumns)){
            if (sum(PatientOverview[PatientOverview$patient==patname,MutColumns[i]])>0){
                print (paste(i,names(PatientOverview)[MutColumns[i]]))  
                PatSpecResistanceCodons<-c(PatSpecResistanceCodons,names(PatientOverview)[MutColumns[i]])
            }}
        #which of the mutations in the patient are RT mutations? 
        PositionsRTImp<-as.numeric(substr(PatSpecResistanceCodons[which(substr(PatSpecResistanceCodons,1,2)!="PI")],3,6))
        RTIpositions<-vector(); for (c in PositionsRTImp){RTIpositions<-c(RTIpositions,which(L$codon==c+99))}
        #which are PRO mutations? 
        PositionsPROImp<-as.numeric(substr(PatSpecResistanceCodons[which(substr(PatSpecResistanceCodons,1,2)=="PI")],3,6))
        PROpositions<-vector(); for (c in PositionsPROImp){PROpositions<-c(PROpositions,which(L$codon==c))}
        #find out which mutations are polymorphic OR relevant for resistance, these should be included in the figure	
        if (ShowSingletons)SitesToDraw<-sort(unique(c(which(L$status!="N"),PROpositions,RTIpositions))) #N = not polymorphic
        if (!ShowSingletons)SitesToDraw<-sort(unique(c(which(L$status=="P"),PROpositions,RTIpositions))) #P = polymorphic
        if (ShowOnlyResSites)SitesToDraw<-sort(unique(c(PROpositions,RTIpositions)))
        PatSpecResistancePositions=c(PROpositions,RTIpositions)
    }	
    
    
    #make a figure
    if (length(days)>1 & length(PatSpecResistanceCodons)>=1
        #Add sth here to make sure we only plot is length(days)>1
        #Add sth here to make sure we only plot id num res mut >1
    ){
        cexpos=1; wl=0.4; wr=0.4; he=0.4; HE=0.4
        #open a png file to make a figure
        widthpng=300+15*length(SitesToDraw)
        heightpng=300+15*numseqs	
        figurefilename=paste("Output/Jan2018Graphs/", patname,"Full.png",sep="");
        if (ShowOnlyResSites) figurefilename=paste("Output/Jan2018Graphs/", patname,"Short.png",sep="");
        png(figurefilename,width=widthpng,height=heightpng,units="px",pointsize=12,bg="white")
        par(mar=c(1,1,2.,0))
        #make empty plot
        plot(1:2,1:2,col="white",ylim=c(-4,(length(seqlabels))+length(days))+2,xlim=c(-1,length(SitesToDraw)+6),ylab="",xlab="sites",yaxt="n",xaxt="n",frame.plot=FALSE)
        title(main=paste("Viral sequences from patient",substr(patname,4,6), "                    "),cex.main = 2,line = -1.)	
        
        if (length(SitesToDraw[SitesToDraw<298])>=3) rect(xleft = 0.,ybottom = -3,xright = length(SitesToDraw[SitesToDraw<298])+0.5,ytop = -1.3,col="black")
        if (length(SitesToDraw[SitesToDraw>=298])>=3) rect(xleft = length(SitesToDraw[SitesToDraw<298]) +0.5,ybottom = -3,xright =length(SitesToDraw)+1  ,ytop = -1.3,col="darkgrey")
        extrasize=0.6; protext="PROTEASE"; RTtext="REVERSE TRANSCRIPTASE"
        if ((length(SitesToDraw[SitesToDraw<298])<=6)|(length(SitesToDraw[SitesToDraw>=298])<=6)){extrasize=0.6;protext="PRO"; RTtext="RT"}
        if (length(SitesToDraw[SitesToDraw<298])>=3) text(1+length(SitesToDraw)/8,-2,protext,cex=cexpos+extrasize,col="white")
        if (length(SitesToDraw[SitesToDraw>=298])>=3) text((min(which(SitesToDraw>297))+length(SitesToDraw))/2,-2,RTtext,cex=cexpos+extrasize,col="white")
        
        #indicate the locations of important resistance mutations (codon number and darkgrey color)
        for (p in PatSpecResistanceCodons){
            if (substr(p,1,2)!="PI"){ # non PI sites
                rect(xleft= which(SitesToDraw == 3*(as.numeric(substr(p,3,6))+99))-2.5, ybottom = length(seqlabels)+length(days)+0.2-0.3, xright = which(SitesToDraw == 3*(as.numeric(substr(p,3,6))+99))+.5, ytop =length(seqlabels)+length(days)+1.8-0.5 , col="darkgrey")
                text(which(SitesToDraw == 3*(as.numeric(substr(p,3,6))+99))-1,length(seqlabels)+length(days)+1-0.3,substr(p,3,6),cex=1.7, col="white")}
            if (substr(p,1,2)=="PI"){ # PI sites
                rect(xleft= which(SitesToDraw == 3*(as.numeric(substr(p,3,6))))-2.5, ybottom = length(seqlabels)+length(days)+0.2-0.3, xright = which(SitesToDraw == 3*(as.numeric(substr(p,3,6))))+.5, ytop =length(seqlabels)+length(days)+1.8-0.5 , col="black")
                text(which(SitesToDraw == 3*(as.numeric(substr(p,3,6))))-1,length(seqlabels)+length(days)+1-0.3,substr(p,3,6),cex=1.7, col="white")}}
        for (p in PatSpecResistancePositions){ #Not sure what this does...
            for (d in 1:length(days)){
                height= range(d - 1 + which(substr(seqlabels,5,7)==days[d]))
                rect(which(SitesToDraw==p)-wl-0.1,height[1]-.5,which(SitesToDraw==p)+wr+0.1,height[2]+0.5,density=-1,col="white",lwd=0, border="darkgrey")}}
        
        #indicate resistance status of codon in light blue
        if (TRUE){
            for (gene in c("PI","RTI")){
                if (gene=="PI") respos= PositionsPROImp
                if (gene=="RTI") respos= PositionsRTImp
                for (cod in respos){	
                    for (seq in seqlabels){
                        num=which(seqlabels==seq)+which(days==substr(seq,5,7))-1
                        num2=which(names(patfasta[,1])==seq)
                        if (gene=="PI") {
                            y<-which(L$codon==cod); x<-which(PImuts$pos==cod)
                            if (length(grep(translate(patfasta[num2,y]),PImuts$mut[x]))>0){co="darkgrey"
                            rect(which(SitesToDraw==y[1])-HE,num-.5,which(SitesToDraw==y[3])+HE,num+.5,density=-1,col=co,lwd=0, border=co)}
                        }
                        
                        if (gene=="RTI"){
                            y<-which(L$codon==cod+99); x<-which(RTImuts$pos==cod)
                            if (length(grep(translate(patfasta[num2,y]),RTImuts$mut[x]))>0){
                                co="darkgrey"
                                rect(which(SitesToDraw==y[1])-HE,num-.5,which(SitesToDraw==y[3])+HE,num+.5,density=-1,col=co,lwd=0, border=co)}}
                        
                    }}
            }
        }
        #add lines between the codons
        #for (x in seq(3.5, length(PatSpecResistancePositions),by=3)){
        #    lines(c(x,x),c(.5,-.5+length(seqlabels)+length(days)),col="white",lwd=3)}
        
        #draw colored rectangles for all other resistance relevant sites	
        for (s in PatSpecResistancePositions){
            #if (s<604|s>606){
                for (seq in seqlabels){
                    num=which(seqlabels==seq)+which(days==substr(seq,5,7))-1
                    num2=which(names(patfasta[,1])==seq)
                    col=brewer.pal(9,"Set3")[s%%9]
                    col=brewer.pal(6,"Set3")[2+which(c("a","c","g","t")==patfasta[num2,s])]
                        if (L$status[s]!="N"&patfasta[num2,s]!=as.character(consensusB)[s]){ # compare to majority day 0 for that patient
                            rect(which(SitesToDraw==s)-HE,num-HE,which(SitesToDraw==s)+HE,num+HE,density=-1,col=col,lwd=0, border=col)
                            }
                        if (L$status[s]!="N"&as.character(consensusB)[s] != L$majoritydayzero[s] & patfasta[num2,s]!=L$majoritydayzero[s]){
                            rect(which(SitesToDraw==s)-HE,num-HE,which(SitesToDraw==s)+HE,num+HE,density=-1,col=col,lwd=0, border=col)
                        }
                    text(which(SitesToDraw==s),num,toupper(patfasta[num2,s]),col="black",cex=cexpos*1)
                   }
        #}
        }
        
        #show segregating sites, both syn and non syn
        for (s in SitesToDraw){#for all multipletons
            if (s>=297) text(which(SitesToDraw==s),-0.5,s-297,srt=90)
            if (s<297)	text(which(SitesToDraw==s),-0.5,s,srt=90)	
            if (length(which(PatSpecResistancePositions==s))==0){
                for (seq in seqlabels){
                    num=which(seqlabels==seq)+which(days==substr(seq,5,7))-1
                    num2=which(names(patfasta[,1])==seq)
                    if (patfasta[num2,s]!=L$majoritydayzero[s]){
                        if(L$synnonsyn[s]=="non")co=brewer.pal(8,"Set3")[8] #"red"
                        if(L$synnonsyn[s]!="non")co=brewer.pal(6,"Set3")[1] 
                        if(L$fourfold[s]==1)co=brewer.pal(6,"Set3")[1] 
                        rect(which(SitesToDraw==s)-HE,num-HE,which(SitesToDraw==s)+HE,num+HE,density=-1,col=co,lwd=0, border=co)}}}
        }
        
        
        #indicate days and samples	
        for (d in 1:length(days)){
            height= range(d - 1 + which(substr(seqlabels,5,7)==days[d]))
            lines(x=c(0,length(SitesToDraw)+1),y=c(height[1]-0.5,height[1]-0.5),col=1)	
            lines(x=c(0,length(SitesToDraw)+1),y=c(height[2]+0.5,height[2]+0.5),col=1)
            lines(x=c(0,0),y=c(height[1]-0.5,height[2]+0.5))
            lines(x=c(length(SitesToDraw)+1,length(SitesToDraw)+1),y=c(height[1]-0.5,height[2]+0.5))
            text(pos=4, length(SitesToDraw)+1.2, mean(height),paste("day",as.numeric(days[d])),cex=cexpos*2,srt=0)}
        
        #mtext("Data from Bacheler et al 2000, Figure related to and by: Pennings, Kryazhimskiy, Wakeley 2013", 4, line=-2)
        
        #close the file
        dev.off() 
    }
    
}

