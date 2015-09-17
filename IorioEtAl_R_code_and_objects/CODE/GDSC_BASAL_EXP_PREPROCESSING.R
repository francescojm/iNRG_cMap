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

library(sROC)
load('DATA/GDSC_basalEXP.ro')

EL_statistic<-function(expression_pattern,ret.pvals=FALSE){
  ep<-expression_pattern
  
  
  CDF<-kCDF(ep,xgrid=ep,adjust=1)
  nep<-CDF$Fhat[match(ep,CDF$x)]
  
  el<-log(nep/(1-nep))
  
  if (ret.pvals){
    pvals<-rep(NA,length(nep))
    pvals[el>=0]<-(1-nep[el>=0])
    pvals[el<0]<-nep[el<0]
    
    return(list(el=el,pvals=pvals))
  }
  
  return(el)
}

basalRanked_lists<-function(medNorm_basalExp){
  cat('Generating rankings...\n')
  ranks<-apply(medNorm_basalExp,MARGIN=2,FUN='order',decreasing=TRUE)
  
  
  nc<-ncol(ranks)
  basal_ranked_lists<-matrix(NA,nrow(medNorm_basalExp),nc)
  
  
  for (i in 1:nc){
    
    basal_ranked_lists[,i]<-rownames(medNorm_basalExp)[ranks[,i]]
  }
  cat('Done!\n')
  
  colnames(basal_ranked_lists)<-colnames(medNorm_basalExp)
  return(basal_ranked_lists)
}

EL_statistics<-function(expression_data,show_progress=TRUE){
  ed<-expression_data
  
  
  nn<-dim(ed)
  EL<-matrix(NA,nn[1],nn[2])
  
  cat('Computing Expression Level statistics...\n')
  if(show_progress){pb <- txtProgressBar(min=1,max=nn[1],style=3)}
  
  for (i in 1:nn[1]){
    EL[i,]<-EL_statistic(ed[i,])
    if(show_progress){setTxtProgressBar(pb, i)}
  }
  
  if(show_progress){
    Sys.sleep(1)
    close(pb)
  }
  cat('Done!\n')
  
  rownames(EL)<-rownames(ed)
  colnames(EL)<-colnames(ed)
  return(EL)
}

