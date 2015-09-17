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


library(beeswarm)

test_pred_ability<-function(multiple_sig_cs,DRUGS,th=0.30,mainTitle,display=TRUE){

  ndrugs<-length(DRUGS)
  
  connected_cell_lines<-cell_lines_connected_to_mult_sig(multiple_sig_cs)
  
  if (display){
    layout(matrix(c(1,1,1,2,3,4), ncol=3, byrow=TRUE),heights=c(1,4))
    par(xpd=NA)
    plot(0,0,col=NA,axes=FALSE,xlab='',ylab='')
    text(0,0,mainTitle,cex=1.5)
  }
  
  for (i in 1:length(DRUGS)){
    res<-genericTtest(SCREENING$IC50s,ALLscreenedCellLines=rownames(SCREENING$IC50s),display=display,
                 specific_cell_lines=connected_cell_lines$NEG,drug=DRUGS[i],
                 labels=c(paste('predicted\nsensitive'),'others'))
    currentLine<-c(res$PVAL,res$deltaMEAN,res$effectSize,res$N1,res$N2)
    if (i == 1){
      totres<-currentLine
   
    }else{
      totres<-rbind(totres,currentLine)
    }
  }
  totres<-cbind(DRUGS,as.character(DRUG_PROPS[DRUGS,1]),as.character(DRUG_PROPS[DRUGS,4]),totres)
  totres<-cbind(rep(mainTitle,length(DRUGS)),totres)
  
  totres<-as.data.frame(totres,row.names=NA)
  colnames(totres)<-c('used signature(s)','drug id','drug','target','p-val','deltaMean','effectSize','N1','N2')
  return(totres)
}
cell_lines_connected_to_mult_sig<-function(multiple_sig_cs,th=0.30){
  
  nsig<-length(multiple_sig_cs)
  
  
  for (i in 1:nsig){
    if (i == 1){
      idxsNEG<-names(which(multiple_sig_cs[[i]]$NCS < 0 & multiple_sig_cs[[i]]$adjP < th))
      idxsPOS<-names(which(multiple_sig_cs[[i]]$NCS > 0 & multiple_sig_cs[[i]]$adjP < th))
    }
    else{
      idxsNEG<-intersect(idxsNEG,names(which(multiple_sig_cs[[i]]$NCS < 0 & multiple_sig_cs[[i]]$adjP < th)))
      idxsPOS<-intersect(idxsNEG,names(which(multiple_sig_cs[[i]]$NCS > 0 & multiple_sig_cs[[i]]$adjP < th)))
    }
  }
  
  NEG<-intersect(idxsNEG,rownames(SCREENING$IC50s))
  POS<-intersect(idxsPOS,rownames(SCREENING$IC50s))

  return(list(NEG=NEG,POS=POS))
}
genericTtest<-function(IC50s,ALLscreenedCellLines,specific_cell_lines,drug=drug_id,display=TRUE,labels=NULL){
  
  ALLscreenedCellLines<-intersect(ALLscreenedCellLines,rownames(IC50s))
  specific_cell_lines<-intersect(specific_cell_lines,rownames(IC50s))
  
  ALLscreenedCellLines<-ALLscreenedCellLines[!is.na(IC50s[ALLscreenedCellLines,as.character(drug)])]
  specific_cell_lines<-specific_cell_lines[!is.na(IC50s[specific_cell_lines,as.character(drug)])]
  
  
  IC50pattern<-IC50s[ALLscreenedCellLines,as.character(drug)]
  names(IC50pattern)<-ALLscreenedCellLines
  
  if (length(labels)==0){
    labels=c('g1','g2')
  }
  XLAB='GDSC cell lines'
  
  N1<-length(which(!is.na(IC50pattern[setdiff(ALLscreenedCellLines,specific_cell_lines)])))
  N2<-length(which(!is.na(IC50pattern[specific_cell_lines])))
  
  if (N1>=2  & N2>=2){
    
    
    TT<-t.test(IC50pattern~(is.element(ALLscreenedCellLines,specific_cell_lines)))
    P<-TT$p.value
    DM<-TT$estimate[2]-TT$estimate[1]
    DM<-DM[[1]]
    
    
    lx <- N1 - 1
    ly <- N2 - 1
    md  <- abs(DM)        ## mean difference (numerator)
    
    
    x<-IC50pattern[setdiff(ALLscreenedCellLines,specific_cell_lines)]
    x<-x[which(!is.na(x))]
    
    y<-IC50pattern[specific_cell_lines]
    y<-y[which(!is.na(y))]
    
    
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation
    
    cd  <- md/csd                        ## cohen's d
    
    if(display){  
      
      
      MAIN<-paste(getDrugName(drug),'\n')
      
      
      
      beeswarm(IC50pattern~(!is.element(ALLscreenedCellLines,specific_cell_lines)),labels=labels,xlab=XLAB,ylab='log(IC50)',
               col=c(rgb(150,0,255,150,maxColorValue=255),
                     rgb(100,100,100,100,maxColorValue=255)),
                     pch = 16,cex=2,
               main=paste(MAIN,
                          'DM = ',format(DM,digits=2),', Eff.S = ',format(cd,digits=2),
                          ', P = ',format(P,digits=2,scientific=TRUE),
                          sep=''),cex.main=1,corral='wrap',cex.lab=1,cex.names=0.5)
      
      
      boxplot(IC50pattern~(!is.element(ALLscreenedCellLines,specific_cell_lines)),add=TRUE,outline=FALSE,col=NA,boxwex=0.6,names=c('',''),lwd=2)
      
      Malt<-mean(IC50pattern[which(is.element(ALLscreenedCellLines,specific_cell_lines))])
      Mwt<-mean(IC50pattern[which(!is.element(ALLscreenedCellLines,specific_cell_lines))])
      
      SDalt<-sd(IC50pattern[which(is.element(ALLscreenedCellLines,specific_cell_lines))])
      SDwt<-sd(IC50pattern[which(!is.element(ALLscreenedCellLines,specific_cell_lines))])
      
      lines(x=c(0.85,1.15),y=c(Malt,Malt),col='red',lwd=5)
      lines(x=c(1.85,2.15),y=c(Mwt,Mwt),col='red',lwd=5)
      
      lines(x=c(0.70,1.30),y=c(Malt+SDalt,Malt+SDalt),col='red',lwd=2)
      lines(x=c(0.70,1.30),y=c(Malt-SDalt,Malt-SDalt),col='red',lwd=2)
      
      lines(x=c(1.70,2.30),y=c(Mwt+SDwt,Mwt+SDwt),col='red',lwd=2)
      lines(x=c(1.70,2.30),y=c(Mwt-SDwt,Mwt-SDwt),col='red',lwd=2)
      
    }
    
  }else{
    P<-NA
    DM<-NA
    cd<-NA
    effectSize<-NA
  }
  
  
  return(list(PVAL=P,deltaMEAN=DM,effectSize=cd,N1=N1,N2=N2))
}