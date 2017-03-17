#FUNCTION checkdata
checkDNAbcd<-function(seq,DistModel="K80") {   
 mysplit<-c();
 mylabels<-c();
 listsp<-data.frame();
 DNAlength<-c();
 mysplit<-strsplit(labels(seq),"_");
 mygen<-c(); for (i in 1:length(mysplit)){ mygen<-c(mygen,mysplit[[i]][1]);}
 mysp<-c(); for (i in 1:length(mysplit)){ mysp<-c(mysp,mysplit[[i]][2]);}
 mylabels<-data.frame(genus=mygen,species=mysp,id=labels(seq));
 listsp<-data.frame(species=sort(unique(paste(mylabels$genus,mylabels$species,sep="_"))),Nseq=NA,Nhap=NA);
 for(i in 1:length(listsp$species)) {
  listsp$Nseq[[i]]<-length(grep(listsp$species[[i]],labels(seq)));  
  if (listsp$Nseq[[i]]>1){
   listsp$Nhap[[i]]<-dim(haplotype(seq[grep(listsp$species[[i]],labels(seq)), ]))[[1]];  
  }
  else {
   listsp$Nhap[[i]]<- listsp$Nseq[[i]];
  }
 }
 mytbl<-c(); 
 for (j in 1:dim(seq)[1]){
  mytbl<-table(as.character(seq[j,]));
  DNAlength[j]<-dim(seq)[2]-sum(mytbl[which(names(mytbl)=="-" |names(mytbl)=="?"| names(mytbl)=="N"| names(mytbl)=="N" |names(mytbl)=="_")]);
 }
 
#INTRA- AND INTERSPECIFIC DISTANCES
 dist<-dist.dna(seq,DistModel,as.matrix=TRUE,pairwise.deletion=TRUE);
 diag(dist)<-NA
 disttmp<-matrix();
 disttmp<-dist;
 disttmp[upper.tri(disttmp,diag=TRUE)]<-NA;
 spdist<-list(); #intra and inter
 for (i in 1:length(listsp$species)){   #intraspecific dist
  if (length(grep(as.character(listsp$species[i]),labels(seq)))>0) {
   x<-matrix();
   x<-disttmp[grep(as.character(listsp$species[i]),labels(seq)),grep(as.character(listsp$species[i]),labels(seq))];
   spdist$intra<-c(spdist$intra,as.vector(as.numeric(x)));
   spdist$intra<-na.omit(spdist$intra);
  }
 }
 for (i in 1:length(listsp$species)){  #interspecific dist
  x<-matrix();
  x<-disttmp[grep(as.character(listsp$species[i]),labels(seq)),-grep(as.character(listsp$species[i]),labels(seq))];
  spdist$inter<-c(spdist$inter,as.vector(as.numeric(x))); 
  spdist$inter<-na.omit(spdist$inter);
 }
 colnames(dist)<-paste(mylabels$genus,mylabels$species,sep="_");
 write.csv(mylabels,"mylabels.csv");
 write.csv(listsp,"listsp.csv");
return(list(mylabels=mylabels, listsp=listsp, DNAlength=DNAlength, dist=dist, spdist=spdist, seq=seq));
}

