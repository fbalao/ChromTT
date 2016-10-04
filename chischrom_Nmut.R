#Analisis para los stem con número de mutaciones por nodo

chischrom_Nmut<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  if (dim(preMedt)==2){
    preMedt[1]<-sum(rowSums(preMed[,5:9], na.rm = T))
  }  else {
    preMedt<-c(sum(rowSums(preMed[,5:9], na.rm = T)), preMedt)
    }  
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  if (dim(Medt)==2){
    preMedt[1]<-sum(rowSums(Med[,5:9], na.rm = T))
  }  else {
    Medt<-c(sum(rowSums(Med[,5:9], na.rm = T)), Medt)
    }
  chidata<-cbind(preMedt, Medt)
  dm<-prod(dim(chidata))
  if(dm!=4){
    restodo<-c(NA,NA,NA) #Check 2x2 table
  } else {
    restodo<-(chisq.test(chidata))
  }
  
  #Chisqtest for tips 
  preMedtips<-tips[tips$agestem>3.4,]
  Medtips<-tips[tips$agestem<3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  if (dim(preMedtips)==2){
    preMedtipst[1]<-sum(rowSums(preMedtips[,5:9], na.rm = T))
  }  else {preMedtipst<-c(sum(rowSums(preMedtips[,5:9], na.rm = T)), preMedtipst)
  }
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)   
  if (dim(Medtips)==2){
    Medtipst[1]<-sum(rowSums(Medtips[,5:9], na.rm = T))
  }  else {
    Medtipst<-c(sum(rowSums(Medtips[,5:9], na.rm = T)), Medtipst)
    }
  chidatatips<-cbind(preMedtipst, Medtipst)
  dm<-prod(dim(chidatatips))
  if(dm!=4){
    restips<-c(NA,NA,NA)
  } else {
      restips<-(chisq.test(chidatatips))
      }

  
  #Chisqtest for nodes 
  
  preMednodes<-nodes[nodes$agestem>3.4,]
  Mednodes<-nodes[nodes$agestem<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  if (dim(preMednodest)==2){
    preMednodest[1]<-sum(rowSums(preMednodes[,5:9], na.rm = T))
  } else {
    preMednodest<-c(sum(rowSums(preMednodes[,5:9], na.rm = T)), preMednodest)}
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)  
  if (dim(Mednodest)==2){Mednodest[1]<-sum(rowSums(Mednodes[,5:9], na.rm = T))
  } else {Mednodest<-c(sum(rowSums(Mednodes[,5:9], na.rm = T)), Mednodest)}
  chidatanodes<-cbind(preMednodest, Mednodest)
  dm<-prod(dim(chidatanodes))
  if(dm!=4){
    resnodes<-c(NA,NA,NA)
  } else {resnodes<-(chisq.test(chidatanodes))
    }
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}


setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()
resultschrom<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom[i,]<-chischrom_Nmut(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom)<-listfile
resultschrom
write.table(resultschrom, file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_stem_Nmutransitions.txt", sep="\t")
