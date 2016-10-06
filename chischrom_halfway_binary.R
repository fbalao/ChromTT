chischrom.halfway<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[((x$agestem - x$age)/2 + x$age) >3.4,]
  Med<-x[((x$agestem - x$age)/2 + x$age)<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  dm<-prod(dim(chidata))
  if(dm!=4){
    restodo<-c(NA,NA,NA)
  }   else {
    restodo<-(chisq.test(chidata))
  }
  
  
  # New method uniform
  Nnodestotal <- dim(x)[1]
  NnodespreMed <- dim(preMed)[1]
  NnodesMed <- dim(Med)[1]
  NmutpreMed <-  NnodespreMed- sum(rowSums(preMed[,5:9], na.rm = T)==0)
  NmutMed <-  NnodesMed- sum(rowSums(Med[,5:9], na.rm = T)==0)
  restodo <-chisq.test(c(NmutpreMed,NmutMed), p=c(NnodespreMed/Nnodestotal,NnodesMed/Nnodestotal ))

  
# Chisqtest for tips HALFWAY
  preMedtips<-tips[((tips$agestem - tips$age)/2 + tips$age) >3.4,]
  Medtips<-tips[((tips$agestem - tips$age)/2 + tips$age) <3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)    
  chidatatips<-cbind(preMedtipst, Medtipst)
  dm<-prod(dim(chidatatips))
  if(dm!=4){
    restips<-c(NA,NA,NA)
  }   else {
    restips<-(chisq.test(chidatatips))
  }
  
  # New method uniform
  Ntipstotal <- dim(x)[1]
  NtipspreMed <- dim(preMed)[1]
  NtipsMed <- dim(Med)[1]
  NmuttipspreMed <-  NtipspreMed- sum(rowSums(preMedtips[,5:9], na.rm = T)==0)
  NmuttipsMed <-  NtipsMed- sum(rowSums(Medtips[,5:9], na.rm = T)==0)
  restips <-chisq.test(c(NmuttipspreMed,NmuttipsMed), p=c(NtipspreMed/Ntipstotal,NtipsMed/Ntipstotal ))

  
  
  #Chisqtest for nodes HALFWAY NODE 
  
  preMednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age) >3.4,]
  Mednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age)<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)    
  chidatanodes<-cbind(preMednodest, Mednodest)
  dm<-prod(dim(chidatanodes))
  if(dm!=4){
    resnodes<-c(NA,NA,NA)
  }   else {
    resnodes<-(chisq.test(chidatanodes))
  }
  
  # New method uniform
  Nnodestotal <- dim(x)[1]
  NnodespreMed <- dim(preMed)[1]
  NnodesMed <- dim(Med)[1]
  NmutnodespreMed <-  NnodespreMed- sum(rowSums(preMednodes[,5:9], na.rm = T)==0)
  NmutnodesMed <-  NnodesMed- sum(rowSums(Mednodes[,5:9], na.rm = T)==0)
  resnodes <-chisq.test(c(NmutnodespreMed,NmutnodesMed), p=c(NnodespreMed/Nnodestotal,NnodesMed/Nnodestotal ))

  
  
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}

#Analisis para los stem con mutaciones 1/0 HALFWAY
setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()

resultschrom.halfway<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom.halfway[i,]<-chischrom.halfway(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom.halfway)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom.halfway)<-listfile
resultschrom.halfway

write.table(round(resultschrom.halfway, 5), file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_halfway_binarytransitions.txt", sep="\t")
