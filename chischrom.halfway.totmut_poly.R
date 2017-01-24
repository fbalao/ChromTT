chischrom.halfway.totmut_poly<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  nodes<-nodes[nodes$name!='N1',]
  tips<-x[-(grep("N[1-9]", x$n)),]
  tips<-tips[tips$name!='N1',]
  preMed<-x[((x$agestem - x$age)/2 + x$age) >3.4,]
  preMed<-preMed[preMed$name!='N1',]
  Med<-x[((x$agestem - x$age)/2 + x$age)<3.4,]
  Med<-Med[Med$name!='N1',]
  Ntotal <- dim(x)[1]-1
  NpreMed <- dim(preMed)[1]
  NMed <- dim(Med)[1]
  NmutpreMed <-  sum(rowSums(preMed[,7:9], na.rm=T))
  NmutMed <-  sum(rowSums(Med[,7:9], na.rm=T))
  if (sum(c(NmutpreMed,NmutMed))!=0){
    restodo <-chisq.test(c(NmutpreMed,NmutMed), p=c(NpreMed/Ntotal,NMed/Ntotal ), simulate.p.value = T, B=999)
  } else {
    restodo<-c(NA,NA,NA)
  }
  
  # Chisqtest for tips HALFWAY
  preMedtips<-tips[((tips$agestem - tips$age)/2 + tips$age) >3.4,]
  Medtips<-tips[((tips$agestem - tips$age)/2 + tips$age) <3.4,]
  Ntipstotal <- dim(tips)[1]
  NtipspreMed <- dim(preMedtips)[1]
  NtipsMed <- dim(Medtips)[1]
  NmuttipspreMed <-  sum(rowSums(preMedtips[,7:9], na.rm=T))
  NmuttipsMed <-  sum(rowSums(Medtips[,7:9], na.rm=T))
  if (sum(c(NmuttipspreMed,NmuttipsMed))!=0){
    restips <-chisq.test(c(NmuttipspreMed,NmuttipsMed), p=c(NtipspreMed/Ntipstotal,NtipsMed/Ntipstotal ), simulate.p.value = T, B=999)
      } else {
    restips<-c(NA,NA,NA)
         }
  
  
  #Chisqtest for nodes HALFWAY NODE 
  
  preMednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age) >3.4,]
  Mednodes<-nodes[((nodes$agestem - nodes$age)/2 + nodes$age)<3.4,]
  Nnodestotal <- dim(nodes)[1]
  NnodespreMed <- dim(preMednodes)[1]
  NnodesMed <- dim(Mednodes)[1]
  NmutnodespreMed <-  sum(rowSums(preMednodes[,7:9], na.rm=T))
  NmutnodesMed <-  sum(rowSums(Mednodes[,7:9], na.rm=T))
  if (sum(c(NmutnodespreMed,NmutnodesMed))!=0){
    resnodes <-chisq.test(c(NmutnodespreMed,NmutnodesMed), p=c(NnodespreMed/Nnodestotal,NnodesMed/Nnodestotal ), simulate.p.value = T, B=999)
  } else {
    resnodes<-c(NA,NA,NA)
  }
  
  chisq<-c(c(NpreMed, NMed, NmutpreMed, NmutMed,restodo[[1]],restodo[[3]]),c(NtipspreMed, NtipsMed,NmuttipspreMed, NmuttipsMed,restips[[1]],restips[[3]]), c(NnodespreMed, NnodesMed,NmutnodespreMed, NmutnodesMed,resnodes[[1]],resnodes[[3]]))
  chisq
}

#Analisis para los stem con mutaciones 1/0 HALFWAY
setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()

resultschrom.halfway_totmut.poly<-matrix(nrow=length(listfile), ncol=18)
for (i in 1:length(listfile)){
  tryCatch(resultschrom.halfway_totmut.poly[i,]<-chischrom.halfway.totmut_poly(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom.halfway_totmut.poly)<-c("NpreMed", "NMed" ,"NmutpreMed", "NmutMed","Chisq_all","p-value_all","NpreMed", "NMed" ,"NmutpreMed", "NmutMed","Chisq_tips","p-value_tips","NpreMed", "NMed" ,"NmutpreMed", "NmutMed","Chisq_nodes","p-value_nodes")
row.names(resultschrom.halfway_totmut.poly)<-listfile
resultschrom.halfway_totmut.poly

write.table(round(resultschrom.halfway_totmut.poly, 5), file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_halfway_totmut_poly.txt", sep="\t")
