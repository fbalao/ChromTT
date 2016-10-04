chischrom<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  dm<-prod(dim(chidata))
  if(dm!=4){res<-c(NA,NA,NA)
  } 
  else 
  res<-(chisq.test(chidata))
  res
}


chischrom_disploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:6], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:6], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  dm<-prod(dim(chidata))
  if(dm!=4){res<-c(NA,NA,NA)
  } 
  else 
  res<-(chisq.test(chidata))
  res
}


chischrom_polyploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,7:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,7:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  dm<-prod(dim(chidata))
  if(dm!=4){res<-c(NA,NA,NA)
  } 
  else 
  res<-(chisq.test(chidata))
  res
}

#Loop para correr el test en todos los archivo y hacer una tabla
setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()
results<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom(file=listfile[i])
  results[i,1]<-res[[1]]
  results[i,2]<-res[[3]]
}


results<-as.data.frame(results)
rownames(results)<- listfile
colnames(results)<-c("Chisq","p-value")
results

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con disploidía
resultsdisplo<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom_disploidy(file=listfile[i])
  resultsdisplo[i,1]<-res[[1]]
  resultsdisplo[i,2]<-res[[3]]
}

resultsdisplo<-as.data.frame(resultsdisplo)
rownames(resultsdisplo)<- listfile
colnames(resultsdisplo)<-c("Chisq","p-value")
resultsdisplo

#Loop para correr el test en todos los archivo y hacer una tabla DATOS con poliploidía
resultspoly<-matrix(nrow=length(listfile), ncol=2)
for (i in 1:length(listfile)){
  res<-chischrom_polyploidy(file=listfile[i])
  resultspoly[i,1]<-res[[1]]
  resultspoly[i,2]<-res[[3]]
}

resultspoly<-as.data.frame(resultspoly)
rownames(resultspoly)<- listfile
colnames(resultspoly)<-c("Chisq","p-value")
resultspoly

write.table(results, file ="ChisqResults_stem_todo.txt", sep="\t", row.names=T)
write.table(resultspoly, file ="ChisqResults_stem_todo_polyploidy.txt", sep="\t", row.names=T)
write.table(resultsdisplo, file ="ChisqResults_stem_todo_disploidy.txt", sep="\t", row.names=T)

