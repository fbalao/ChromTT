chischrom.halfway.totmut_rates<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[((x$agestem - x$age)/2 + x$age) >3.4,]
  Med<-x[((x$agestem - x$age)/2 + x$age)<3.4,]
  Nnodestotal <- dim(x)[1]
  NnodespreMed <- dim(preMed)[1]
  NnodesMed <- dim(Med)[1]
  NmutpreMed <-  sum(rowSums(preMed[,5:9], na.rm=T))
  NmutMed <-  sum(rowSums(Med[,5:9], na.rm=T))
  if (sum(c(NmutpreMed,NmutMed))!=0){
    restodo <-chisq.test(c(NmutpreMed,NmutMed), p=c(NnodespreMed/Nnodestotal,NnodesMed/Nnodestotal ), simulate.p.value = T, B=999)
  } else {
    restodo<-c(NA,NA,NA)
  }
  
  # Chisqtest for tips HALFWAY
  preMedtips<-tips[((tips$agestem - tips$age)/2 + tips$age) >3.4,]
  Medtips<-tips[((tips$agestem - tips$age)/2 + tips$age) <3.4,]
  Ntipstotal <- dim(tips)[1]
  NtipspreMed <- dim(preMedtips)[1]
  NtipsMed <- dim(Medtips)[1]
  NmuttipspreMed <-  sum(rowSums(preMedtips[,5:9], na.rm=T))
  NmuttipsMed <-  sum(rowSums(Medtips[,5:9], na.rm=T))
  if (sum(c(NmuttipspreMed,NmuttipspreMed))!=0){
    restips <-chisq.test(c(NmuttipspreMed,NmuttipspreMed), p=c(NtipspreMed/Ntipstotal,NtipsMed/Ntipstotal ), simulate.p.value = T, B=999)
  } else {
    restips<-c(NA,NA,NA)
  }
  
  
  #Chisqtest for nodes HALFWAY NODE 
  
  preMednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age) >3.4,]
  Mednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age)<3.4,]
  Nnodestotal <- dim(nodes)[1]
  NnodespreMed <- dim(preMednodes)[1]
  NnodesMed <- dim(Mednodes)[1]
  NmutnodespreMed <-  sum(rowSums(preMednodes[,5:9], na.rm=T))
  NmutnodesMed <-  sum(rowSums(Mednodes[,5:9], na.rm=T))
  if (sum(c(NmutnodespreMed,NmutnodesMed))!=0){
    resnodes <-chisq.test(c(NmutnodespreMed,NmutnodesMed), p=c(NnodespreMed/Nnodestotal,NnodesMed/Nnodestotal ), simulate.p.value = T, B=999)
  } else {
    resnodes<-c(NA,NA,NA)
  }
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}

#Analisis para los stem con mutaciones 1/0 HALFWAY
setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()

resultschrom.halfway_totmut.rates<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom.halfway_totmut.rates[i,]<-chischrom.halfway.totmut_rates(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom.halfway_totmut.rates)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom.halfway_totmut.rates)<-listfile
resultschrom.halfway_totmut.rates

write.table(round(resultschrom.halfway_totmut.rates, 5), file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_halfway_totmut_rates.txt", sep="\t")
