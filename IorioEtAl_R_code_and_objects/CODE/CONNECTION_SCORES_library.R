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

library(mixtools)

combine_2CS<-function(CS1,CS2,printToFile=FALSE,fn=''){
  
  connectedDRUGS<-colnames(DRUG_PRLs)
  
  RESULTS<-cbind(CS1$CS[connectedDRUGS],CS1$Pval[connectedDRUGS],100*CS1$adjP[connectedDRUGS],CS1$NCS[connectedDRUGS],
                 CS2$CS[connectedDRUGS],CS2$Pval[connectedDRUGS],100*CS2$adjP[connectedDRUGS],CS2$NCS[connectedDRUGS],
                 rowMeans(cbind(CS1$NCS[connectedDRUGS],CS2$NCS[connectedDRUGS])))
  RESULTS<-RESULTS[order(RESULTS[,9],decreasing=TRUE),]
  
  colnames(RESULTS)<-c('cons S CS','cons S pvalue','cons S fdr %','cons S NCS',
                       'incons S CS','incons S pvalue','incons S fdr %','incons S NCS',
                       'avg NCS')
  
  if (printToFile){
    RESULTS<-cbind(rownames(RESULTS),RESULTS)
    colnames(RESULTS)[1]<-'DRUG'
    write.table(RESULTS,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',fn,'_refined_neighborhood.txt',sep=''))
  }
  return(RESULTS)
}
combine_3CS<-function(CS1,CS2,CS3,previousNeighBr='',printToFile=FALSE,fn=''){
  connectedDRUGS<-previousNeighBr
  
  RESULTS<-cbind(CS1$NCS[connectedDRUGS],
                 CS2$NCS[connectedDRUGS],
                 CS3$NCS[connectedDRUGS],
                 rowMeans(cbind(CS1$NCS[connectedDRUGS],CS2$NCS[connectedDRUGS],CS3$NCS[connectedDRUGS])))
  
  RESULTS<-RESULTS[order(RESULTS[,4],decreasing=TRUE),]
  
  colnames(RESULTS)<-c('P/PI cons S NCS',
                       'P/PI incons S NCS',
                       'MST S NCS',
                       'avg NCS')
  
  if (printToFile){
    RESULTS<-cbind(rownames(RESULTS),RESULTS)
    colnames(RESULTS)[1]<-'DRUG'
    write.table(RESULTS,quote=FALSE,sep='\t',row.names=FALSE,file=paste('OUTPUT/',fn,'_neighborhood.txt',sep=''))
  }
  return(RESULTS)
}
cMap_CS<-function(ranked_list,opsig1,returnRS=FALSE){
  ESUP1<-qES(ranked_list,opsig1$UP,display=FALSE,returnRS=returnRS)
  ESDOWN1<-qES(ranked_list,opsig1$DOWN,display=FALSE,returnRS=returnRS)
  
  if (returnRS){
    RSUP<-ESUP1$RS
    RSDOWN<-ESDOWN1$RS
    
    ESUP1<-ESUP1$ES
    ESDOWN1<-ESDOWN1$ES
    
    TES1<-(ESUP1-ESDOWN1)/2
    
    return(list(TES=TES1,ESUP=ESUP1,ESDOWN=ESDOWN1,RSUP=RSUP,RSDOWN=RSDOWN))
  }
  
  TES1<-(ESUP1-ESDOWN1)/2
  
  return(TES1)
}

CS<-function(signature,RANKED_LISTS,show_progress=TRUE){
  
  signature<-list(UP=signature$seedUPreg$ProbeSets,            
                  DOWN=signature$seedDOWNreg$ProbeSets)
  
  ns<-ncol(RANKED_LISTS)
  
  CS<-rep(NA,ns)
  Pvals<-rep(NA,ns)
  
  names(CS)<-colnames(RANKED_LISTS)
  names(Pvals)<-colnames(RANKED_LISTS)
  
  cat('simulating null model\n')
  mixmdl<-est_emp_Cs(signature,10000,RANKED_LISTS,show_progress=show_progress)
  cat('done!\n')
  
  cat('computing connectivity scores\n')
  
  if(show_progress){
    pb <- txtProgressBar(min=1,max=ns,style=3)
  }
  
  for (i in 1:ns){
    CS[i]<-cMap_CS(RANKED_LISTS[,i],signature)
    Pvals[i]<-pnormmix(CS[i], mixmdl)
    if(show_progress){
      setTxtProgressBar(pb, i)
    }
  }
 
  if(show_progress){
    Sys.sleep(1)
    close(pb)
  }
  
  NCS<-rep(NA,length(CS))
  NCS[CS>=0]<-CS[CS>=0]/max(mixmdl$mu)
  NCS[CS<0]<--CS[CS<0]/min(mixmdl$mu)
  
  names(NCS)<-names(CS)
  
  res<-list(CS=CS,Pval=Pvals,adjP=p.adjust(Pvals,method='fdr'),NCS=NCS)
  
  cat('Done!\n')
  
  return(res)
}
qES<-function(RANKEDLIST,REGULON,display=TRUE,returnRS=FALSE){
  
  REGULON<-intersect(as.character(REGULON),RANKEDLIST)
  
  HITS<-is.element(RANKEDLIST,REGULON)+0
  
  hitCases<-cumsum(HITS)
  missCases<-cumsum(1-HITS)
  
  N<-length(RANKEDLIST)
  NR<-length(REGULON)
  
  Phit<-hitCases/NR
  Pmiss<-missCases/(N-NR)
  
  m<-max(abs(Phit-Pmiss))
  t<-which(abs(Phit-Pmiss)==m)
  
  if (length(t)>1){t<-t[1]}
  peak<-t
  ES<-Phit[t]-Pmiss[t]
  RS<-Phit-Pmiss
  
  if (display){
    if (ES>=0){c<-"red"}else{c<-"green"}
    plot(0:N,c(0,Phit-Pmiss),col=c,type="l",xlim=c(0,N),ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))),xaxs="i",bty="l",axes=FALSE,
         xlab="Gene Rank Position",ylab="Running Sum")
    par(new=TRUE)
    plot(0:N,rep(0,N+1),col='gray',type="l",new=FALSE,xlab="",ylab="",ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))))
    axis(side=2)
    
  }
  
  if (returnRS){
    POSITIONS<-which(HITS==1)
    names(POSITIONS)<-RANKEDLIST[which(HITS==1)]
    
    POSITIONS<-POSITIONS[order(names(POSITIONS))]
    names(POSITIONS)<-names(POSITIONS)[order(names(POSITIONS))]
    
    return(list(ES=ES,RS=RS,POSITIONS=POSITIONS,PEAK=t))
  } else {return(ES)}
}

est_emp_Cs<-function(signature,nt,RANKED_LISTS,show_progress=TRUE){
  EMP_CS<-rep(NA,nt)
  
  ng<-nrow(RANKED_LISTS)
  
  if (show_progress){
    pb <- txtProgressBar(min=1,max=nt,style=3)
  }
  
  for (i in 1:nt){
    EMP_CS[i]<-cMap_CS(RANKED_LISTS[sample(1:ng,ng),1],signature)
    if (show_progress){
      setTxtProgressBar(pb, i)
    }
  }
  if(show_progress){
    Sys.sleep(1)
    close(pb)
  }
  mixmdl = normalmixEM(EMP_CS,k=3,verb=FALSE)
  
  return(mixmdl)
}

pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    if (x>=0){
      lambda[component]*pnorm(-x,mean=-mixture$mu[component],sd=mixture$sigma[component],lower.tail=TRUE)
    }else {
      lambda[component]*pnorm(x,mean=mixture$mu[component],sd=mixture$sigma[component],lower.tail=TRUE)
    }
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(sum(pnorms))
}

getDrugName<-function(id){
  return(as.character(DRUG_PROPS[id,'DRUG_NAME']))
}
getDrugTarget<-function(id){
  return(as.character(DRUG_PROPS[id,'PUTATIVE_TARGET']))
}
