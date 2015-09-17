#########################################################################################################
#                                                                                                       #
#     This code is part of the software enclosed to the paper entitled "A semi-supervised approach for  #
#     refining transcriptional signatures of drug response and repositioning predictions",              #
#     by Francesco Iorio et al, submitted as research paper to PLoS ONE.                                #
#                                                                                                       #
#     Copyright (c) 2014 - 2019, EMBL - European Bioinformatics Institute                               #
#     Author: Francesco Iorio (iorio@ebi.ac.uk)                                                         #
#     Distributed under the GPLv3 License.                                                              #
#     See accompanying file LICENSE.txt or copy at http://www.gnu.org/licenses/gpl-3.0.html             #
#                                                                                                       #
#     Paper website: http://www.ebi.ac.uk/~iorio/PLoS_ONE_Submission                                    #
#                                                                                                       #
#########################################################################################################


library(pheatmap)

load('DATA/DRUG_DISTANCES.ro')
load('DATA/DRUG_COMMUNITIES.ro')
load('DATA/enrichedMoas.ro')
load('DATA/DRUG_PRLs.ro')
load('DATA/affy_ps_annotation.ro')


DNquery<-function(seed='paclitaxel',distTh=0.8065,printToFile=FALSE){
  neighborDistances<-sort(DRUG_DISTANCES[seed,])
  neighborDistances<-neighborDistances[which(neighborDistances<distTh)]
  
  seedId<-which(names(neighborDistances)==seed)
  neighborDistances<-neighborDistances[-seedId]
  
  neighbors<-names(neighborDistances)
  
  WholeDistanceSet<-c(DRUG_DISTANCES)
  WholeDistanceSet<-WholeDistanceSet[which(WholeDistanceSet>0)]
  
  quantiles<-rep(NA,length(neighbors))
  names(quantiles)<-neighbors
  
  drugCommunities<-DRUG_COMMUNITIES[neighbors,1]
  names(drugCommunities)<-neighbors
  communityOccurrence<-rep(NA,length(drugCommunities))
  communityCardinality<-rep(NA,length(drugCommunities))
  
  flag<-0
  for (i in neighbors){
    flag<-flag+1
    quantiles[i]<-length(which(WholeDistanceSet<=neighborDistances[i]))/length(WholeDistanceSet)*100
    communityOccurrence[flag]<-length(which(drugCommunities[1:flag]==drugCommunities[i]))
    communityCardinality[flag]<-length(which(DRUG_COMMUNITIES[,1]==drugCommunities[i]))  
  }
  
  totalNdrugs<-rep(nrow(DRUG_COMMUNITIES),length(neighbors))
  
  subNeighOrder<-1:length(neighbors)
  PVALS<-phyper(communityOccurrence-1,communityCardinality,totalNdrugs-communityCardinality,subNeighOrder,lower.tail=FALSE)
  
  PVALS[which(communityOccurrence<2)]<-NA
  
  a.pval<-p.adjust(PVALS[!is.na(PVALS)],'fdr')
  ADJ.PVAL<-rep(NA,length(drugCommunities))
  ADJ.PVAL[!is.na(PVALS)]<-a.pval
  
  MOAs<-enrichedMOAs[as.character(drugCommunities)]
  
  neighborhood<-as.data.frame(cbind(neighborDistances,quantiles,neighbors,drugCommunities,subNeighOrder,communityOccurrence,communityCardinality,totalNdrugs,PVALS,
                                    ADJ.PVAL,MOAs))
  
  
  
  colnames(neighborhood)<-c('D','quantile %','Drug','C id','order','C Occ','C card','Total #drugs','C Overrep p-val','Adj p-val','MOAs')
  
  if(printToFile){
    write.table(neighborhood,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_DN_neighborhood.txt',sep=''))
  }
  
  return(neighborhood)
}

DeriveSingleSignature<-function(seed='paclitaxel'){
  
  seedUP<-DRUG_PRLs[1:250,seed]
  
  nprobes<-nrow(DRUG_PRLs)

  current_percentile<-100*match(seedUP,DRUG_PRLs[,seed])/nprobes  
  percentiles<-current_percentile
  
  Signature<-cbind(seedUP,percentiles)
  colnames(Signature)[1]<-'ProbeSets'
  rownames(Signature)<-Signature[,1]
  UP<-as.data.frame(Signature)
  
  seedDOWN<-DRUG_PRLs[nprobes:(nprobes-250+1),seed]
  
  current_percentile<-100*match(seedDOWN,DRUG_PRLs[,seed])/nprobes  
  percentiles<-current_percentile
  
  Signature<-cbind(seedDOWN,percentiles)
  colnames(Signature)[1]<-'ProbeSets'
  rownames(Signature)<-Signature[,1]
  
  DOWN<-as.data.frame(Signature)
  return(list(seedUPreg=UP,seedDOWNreg=DOWN))
  
  
  
}
DeriveInConsistentSignature<-function(seed='paclitaxel',otherCompounds=c('MG-132','celastrol','5224221'),PTH=30,FUZZYNESS=2,printToFile=FALSE){
  
  seedUP<-DRUG_PRLs[1:250,seed]
  
  nprobes<-nrow(DRUG_PRLs)
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedUP,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedUP
  
  inconsistency<-(rowSums(percentiles>=(100-PTH)))>=FUZZYNESS
  inconsistentSignature<-which(inconsistency)
  
  inconsistentSignature<-percentiles[inconsistentSignature,]
  
  seedPercentiles<-100*match(rownames(inconsistentSignature),DRUG_PRLs[,seed])/nprobes
  
  inconsistentSignature<-cbind(seedPercentiles,inconsistentSignature)
  colnames(inconsistentSignature)[1]<-seed
  
  inconsistentSignature<-cbind(rownames(inconsistentSignature),inconsistentSignature)
  colnames(inconsistentSignature)[1]<-'ProbeSets'
  
  UP<-as.data.frame(inconsistentSignature)
  
  if(printToFile){
    
    write.table(inconsistentSignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                                      '_inconsistentSignatureUP.txt',sep=''))
  }
  
  
  seedDOWN<-DRUG_PRLs[nprobes:(nprobes-250+1),seed]
  
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedDOWN,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedDOWN
  
  inconsistency<-(rowSums(percentiles<=PTH))>=FUZZYNESS
  inconsistentSignature<-which(inconsistency)
  
  inconsistentSignature<-percentiles[inconsistentSignature,]
  
  seedPercentiles<-100*match(rownames(inconsistentSignature),DRUG_PRLs[,seed])/nprobes
  
  inconsistentSignature<-cbind(seedPercentiles,inconsistentSignature)
  colnames(inconsistentSignature)[1]<-seed
  
  inconsistentSignature<-cbind(rownames(inconsistentSignature),inconsistentSignature)
  colnames(inconsistentSignature)[1]<-'ProbeSets'
  
  DOWN<-as.data.frame(inconsistentSignature)
  
  if(printToFile){
    
    write.table(inconsistentSignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                                      '_inconsistentSignatureDOWN.txt',sep=''))
  }
  
  return(list(seedUPreg=UP,seedDOWNreg=DOWN))
  
  
  
}
DeriveConsistentSignature<-function(seed='paclitaxel',otherCompounds=c('MG-132','celastrol','5224221'),PTH=30,FUZZYNESS=2,printToFile=FALSE){
  
  seedUP<-DRUG_PRLs[1:250,seed]
  
  nprobes<-nrow(DRUG_PRLs)
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedUP,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedUP
  
  inconsistency<-(rowSums(percentiles>=(100-PTH)))>=FUZZYNESS
  consistentSignature<-which(!inconsistency)
  
  consistentSignature<-percentiles[consistentSignature,]
  
  seedPercentiles<-100*match(rownames(consistentSignature),DRUG_PRLs[,seed])/nprobes
  
  consistentSignature<-cbind(seedPercentiles,consistentSignature)
  colnames(consistentSignature)[1]<-seed
  
  consistentSignature<-cbind(rownames(consistentSignature),consistentSignature)
  colnames(consistentSignature)[1]<-'ProbeSets'
  
  UP<-as.data.frame(consistentSignature)
  
  if(printToFile){
    
    write.table(consistentSignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                                    '_consistentSignatureUP.txt',sep=''))
  }
  
  
  seedDOWN<-DRUG_PRLs[nprobes:(nprobes-250+1),seed]
  
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedDOWN,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedDOWN
  
  inconsistency<-(rowSums(percentiles<=PTH))>=FUZZYNESS
  consistentSignature<-which(!inconsistency)
  
  consistentSignature<-percentiles[consistentSignature,]
  
  seedPercentiles<-100*match(rownames(consistentSignature),DRUG_PRLs[,seed])/nprobes
  
  consistentSignature<-cbind(seedPercentiles,consistentSignature)
  colnames(consistentSignature)[1]<-seed
  
  consistentSignature<-cbind(rownames(consistentSignature),consistentSignature)
  colnames(consistentSignature)[1]<-'ProbeSets'
  
  DOWN<-as.data.frame(consistentSignature)
  
  if(printToFile){
    
    write.table(consistentSignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                                    '_consistentSignatureDOWN.txt',sep=''))
  }
  
  return(list(seedUPreg=UP,seedDOWNreg=DOWN))
  
  
  
  
}
DeriveMSTSignature<-function(seed='paclitaxel',otherCompounds=c('albendazole','fenbendazole','nocodazole','parbendazole'),PTH=25,FUZZYNESS=4,printToFile=FALSE){
  
  nprobes<-nrow(DRUG_PRLs)
  
  UpLim<-round(nprobes*PTH/100)
  
  seedUP<-DRUG_PRLs[1:UpLim,seed]
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedUP,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedUP
  
  inconsistency<-(rowSums(percentiles>=(100-PTH)))>=FUZZYNESS
  MSTsignature<-which(inconsistency)
  
  MSTsignature<-percentiles[MSTsignature,]
  
  seedPercentiles<-100*match(rownames(MSTsignature),DRUG_PRLs[,seed])/nprobes
  
  MSTsignature<-cbind(seedPercentiles,MSTsignature)
  colnames(MSTsignature)[1]<-seed
  
  MSTsignature<-cbind(rownames(MSTsignature),MSTsignature)
  colnames(MSTsignature)[1]<-'ProbeSets'
  
  UP<-as.data.frame(MSTsignature)
  
  if(printToFile){
    
    write.table(MSTsignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                             '_MST_UP.txt',sep=''))
  }
  
  
  seedDOWN<-DRUG_PRLs[nprobes:(nprobes-UpLim+1),seed]
  
  
  for (i in 1:length(otherCompounds)){
    
    current_percentile<-100*match(seedDOWN,DRUG_PRLs[,otherCompounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-otherCompounds
  rownames(percentiles)<-seedDOWN
  
  inconsistency<-(rowSums(percentiles<=PTH))>=FUZZYNESS
  MTDSignature<-which(inconsistency)
  
  MTDSignature<-percentiles[MTDSignature,]
  
  seedPercentiles<-100*match(rownames(MTDSignature),DRUG_PRLs[,seed])/nprobes
  
  MTDSignature<-cbind(seedPercentiles,MTDSignature)
  colnames(MTDSignature)[1]<-seed
  
  MTDSignature<-cbind(rownames(MTDSignature),MTDSignature)
  colnames(MTDSignature)[1]<-'ProbeSets'
  
  DOWN<-as.data.frame(MTDSignature)
  
  if(printToFile){
    
    write.table(MTDSignature,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                                                             '_MST_DOWN.txt',sep=''))
  }
  
  return(list(seedUPreg=UP,seedDOWNreg=DOWN))
  
  
  
  
}
percHeatMaps<-function(probes,seed,otherCompounds,printToFile=FALSE){
  
  compounds<-c(seed,otherCompounds)
  nprobes<-nrow(DRUG_PRLs)
  
  for (i in 1:length(compounds)){
    
    current_percentile<-100*match(probes,DRUG_PRLs[,compounds[i]])/nprobes
    
    if (i == 1){
      percentiles<-current_percentile
    }else{
      percentiles<-cbind(percentiles,current_percentile)
    }
  }
  
  colnames(percentiles)<-compounds
  
  
  
  if (printToFile){
    pheatmap(filename=paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                            '_percHeatMap.png',sep=''),border_color=NA,
             percentiles,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(colors=c('darkred','white','darkblue'))(100))
  }else{
    pheatmap(percentiles,border_color=NA,cluster_rows=FALSE,cluster_cols=FALSE,color=colorRampPalette(colors=c('darkred','white','darkblue'))(100))  
  }
  
}
plotRunningSums<-function(consistentSigTable,inconsistentSigTable,
                          seed='paclitaxel',otherCompounds=c('MG-132','celastrol','5224221'),
                          printToFile=FALSE){
  
  if(printToFile){
    png(width=1024,height=300,paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                    '_RS.png',sep=''))
  }
  
  Consistent_signature1<-list(UP=consistentSigTable$seedUPreg$ProbeSets,
                              DOWN=consistentSigTable$seedDOWNreg$ProbeSets)
  
  Inconsistent_signature1<-list(UP=inconsistentSigTable$seedUPreg$ProbeSets,
                                DOWN=inconsistentSigTable$seedDOWNreg$ProbeSets)
  
  compounds<-c(seed,otherCompounds)
  ncompounds<-length(compounds)
  
  
  layout(matrix(c(1:(2*ncompounds),rep(2*ncompounds+1,ncompounds)), ncol=ncompounds, byrow=TRUE), heights=c(4, 4, 2))
  
  par(mar=c(2,5,4,4))
  
  for (i in 1:ncompounds){
    if(i==1) {ylab='ES RS consistent S'} else {ylab=''}
    ConsSig_CS<-cMap_CS(DRUG_PRLs[,compounds[i]],Consistent_signature1,returnRS=TRUE)
    plot(ConsSig_CS$RSUP,type='l',ylab=ylab,xlab='',col='darkred',ylim=c(-1,1),lwd=3,main=compounds[i],cex.main=2,cex.lab=1.5)
    peakUP<-which(abs(ConsSig_CS$RSUP)==max(abs(ConsSig_CS$RSUP)))[1]
    peakDOWN<-which(abs(ConsSig_CS$RSDOWN)==max(abs(ConsSig_CS$RSDOWN)))[1]
    
    par(new=TRUE)
    plot(ConsSig_CS$RSDOWN,type='l',ylab='',xlab='',col='darkblue',ylim=c(-1,1),lwd=2.5,axes=FALSE)
    abline(h=0,lty=2,col='darkgray')
    abline(v=peakUP,lty=3,col='darkred',lwd=3)
    abline(v=peakDOWN,lty=3,col='darkblue',lwd=3)
    
  }
  
  for (i in 1:ncompounds){
    if(i==1) {ylab='ES RS inconsistent S'} else {ylab=''}
    InconsSig_CS<-cMap_CS(DRUG_PRLs[,compounds[i]],Inconsistent_signature1,returnRS=TRUE)
    plot(InconsSig_CS$RSUP,type='l',ylab=ylab,xlab='rank position',col='darkred',ylim=c(-1,1),lwd=3,main=compounds[i],cex.main=2,cex.lab=1.5)
    peakUP<-which(abs(InconsSig_CS$RSUP)==max(abs(InconsSig_CS$RSUP)))[1]
    peakDOWN<-which(abs(InconsSig_CS$RSDOWN)==max(abs(InconsSig_CS$RSDOWN)))[1]
    
    par(new=TRUE)
    plot(InconsSig_CS$RSDOWN,type='l',ylab='',xlab='',col='darkblue',ylim=c(-1,1),lwd=3,axes=FALSE)
    abline(h=0,lty=2,col='darkgray')
    abline(v=peakUP,lty=3,col='darkred',lwd=3)
    abline(v=peakDOWN,lty=3,col='darkblue',lwd=3)
  }
  
  plot(0,0,col='white',axes=FALSE,xlab='',ylab='')
  legend('left',col=c('darkred','darkblue'),legend=c('up-regulated part','down-regulated part'),bty='n',cex=2,horiz=TRUE,lty=1,lwd=3)
  legend('right',col=c('darkred','darkblue'),legend=c('positive peak','negative peak'),bty='n',cex=2,horiz=TRUE,lty=3,lwd=3)
  
  if(printToFile){
    dev.off()
  }
  
}
plotRunningSumsMST<-function(MSTsignatureTable,
                             seed='paclitaxel',
                             otherCompounds=c('albendazole','fenbendazole','nocodazole','parbendazole'),
                             printToFile=FALSE){
  
  if(printToFile){
    png(width=1024,height=300,paste('OUTPUT/',seed,'_',paste(otherCompounds,collapse=','),
                                    '_RS.png',sep=''))
  }
  
  MST_signature<-list(UP=MSTsignatureTable$seedUPreg$ProbeSets,
                      DOWN=MSTsignatureTable$seedDOWNreg$ProbeSets)
  
  
  compounds<-c(seed,otherCompounds)
  ncompounds<-length(compounds)
  
  
  layout(matrix(c(1:ncompounds,rep(ncompounds+1,ncompounds)), ncol=ncompounds, byrow=TRUE), heights=c(4, 2))
  
  par(mar=c(2,5,4,4))
  
  for (i in 1:ncompounds){
    if(i==1) {ylab='Microtubule Stabiliser Signature'} else {ylab=''}
    MSTsig_CS<-cMap_CS(DRUG_PRLs[,compounds[i]],MST_signature,returnRS=TRUE)
    plot(MSTsig_CS$RSUP,type='l',ylab=ylab,xlab='',col='darkred',ylim=c(-1,1),lwd=3,main=compounds[i],cex.main=2,cex.lab=1.5)
    peakUP<-which(abs(MSTsig_CS$RSUP)==max(abs(MSTsig_CS$RSUP)))[1]
    peakDOWN<-which(abs(MSTsig_CS$RSDOWN)==max(abs(MSTsig_CS$RSDOWN)))[1]
    
    par(new=TRUE)
    plot(MSTsig_CS$RSDOWN,type='l',ylab='',xlab='',col='darkblue',ylim=c(-1,1),lwd=2.5,axes=FALSE)
    abline(h=0,lty=2,col='darkgray')
    abline(v=peakUP,lty=3,col='darkred',lwd=3)
    abline(v=peakDOWN,lty=3,col='darkblue',lwd=3)
    
  }
  
  plot(0,0,col='white',axes=FALSE,xlab='',ylab='')
  legend('left',col=c('darkred','darkblue'),legend=c('up-regulated part','down-regulated part'),bty='n',cex=2,horiz=TRUE,lty=1,lwd=3)
  legend('right',col=c('darkred','darkblue'),legend=c('positive peak','negative peak'),bty='n',cex=2,horiz=TRUE,lty=3,lwd=3)
  
  if(printToFile){
    dev.off()
  }
  
}