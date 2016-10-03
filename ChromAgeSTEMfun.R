# Esta función toma los resultados de Chromevol y los resume en una tabla
# Se necesita el arbol original (Nexus tree) y los archivos *_mlAncestors.tree y *_expectations.txt 
# del mejor modelo estimado de ChromEvol
# La estructura de los archivos debe ser como estan en Drive. En la carpeta principal el tree
# en la carpeta OUT/1_Best el resto de archivos
# Tiene dos argumentos name: Nombre general del archivo (género) y path: la dirección 
# a la carpeta con los datos



#ChromR(name="Reseda",path="/home/fbalao/Datos/R/Rpackages/ChromTT/Reseda", best="1_Best_CONST_RATE_DEMI_EST")

ChromR<- function(name, path, bestmodel){
#Cargar paquetes necesarios. si no los tienes instalar con install.packages
  library(phytools)
  library(data.table)
  library(plyr)
#Fija el directorio
setwd(path)
#lee el arbol original para obtener las edades de los nodos
tree<-paste(name,".tree", sep="")
arbol<-read.nexus(tree)
arbol$tip.label<-gsub("\\'|\\]", "", arbol$tip.label)
ntip<-length(arbol$tip.label)
nnodes<-length(node.height(arbol))
ninodes<-nnodes-ntip
nodenames<-paste("N", 1:ninodes, sep="")
nodeage<- branching.times(arbol)


#Obtener la edad del stem asociado a cada nodo
stem1<-arbol$edge
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
agestem<-cbind(as.numeric.factor(revalue(factor(stem1[,1]), nodeage)), stem1[,2])
rootnode<-ntip+1
rootage<-nodeage[as.character(rootnode)]
agestem<-rbind(agestem, c(rootage,rootnode))
agestem<-agestem[order(agestem[,2]),]

#
names(nodeage)<-nodenames
#plot(arbol)
#nodelabels(nodenames)
time<-c(rep(0,ntip),nodeage)

#Obtener los números de cromosomas en cada nodo/tip *_mlAncestors.tree
path2<- paste(path,"/OUT/", bestmodel,"/", name, "_mlAncestors.tree", sep="")
chromtree<-readLines(path2)
chromtree2<-gsub("\\[|\\]", "", chromtree)
chromnumtree<-read.tree(text=chromtree2)
chromnum<-data.frame(name=c(chromnumtree$tip.label, chromnumtree$node.label))
chromnum$chrom<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 2)
chromnum$name<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 1)

chromtable<-data.frame(name=c(arbol$tip.label, nodenames), age=time, agestem=agestem[,1])

chromtable2<-merge(chromtable, chromnum, by="name")
chromtable2<-chromtable2[with(chromtable2, order(-age)), ]

#scatter.smooth(-chromtable2$age, chromtable2$chrom, pch=16, lwd=2, col=1)

# Leer la tabla de mutaciones *_expectations.txt con paquete data.table
path3<- paste(path,"/OUT/", bestmodel,"/", name, "_expectations.txt", sep="")
mutations <- fread(path3, 
                   skip="#ALL EVENTS EXPECTATIONS PER NODE", nrows=nnodes-1,  header=T) # No importa el error
colnames(mutations)[1]<- "name"
D1<-data.frame("N1",0,0,0,0,0)
colnames(D1)<- colnames(mutations)
mutations<-rbind(D1, mutations, fill=T)
mutations<-mutations[with(mutations, order(name)), ]


chromtable3<-merge(chromtable2, mutations, by="name")

#Redondeo de las mutaciones
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


chromtable3[,5:dim(chromtable3)[2]]<-round2(chromtable3[,5:dim(chromtable3)[2]],0)
chromtable3

#Exportar la tabla a archivo tabulado
name2<-paste(name, "_summary.txt", sep="")
write.table(as.matrix(chromtable3), file=name2, sep="\t", fileEncoding = "utf8", row.names = F)
}


#plot(-chromtable2$age, chromtable2$chrom, type="p")
#cbind(-chromtable2$age, as.numeric(chromtable2$chrom))

#Chi2 test
chischrom<-function(file){
x<-read.table(file, header=T, strip.white = T)
preMed<-x[x$agestem>3.4,]
Med<-x[x$agestem<3.4,]

preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
chidata<-cbind(preMedt, Medt)
print(chisq.test(chidata))
}


chischrom_disploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,5:6], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:6], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
}


chischrom_polyploidy<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  
  preMedt<-table(rowSums(preMed[,7:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,7:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
  print(chisq.test(chidata))
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



#Función para hacer los Chisq test en todos los nodos, solo tips y solo nodos internos
chischrom2<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  chidata<-cbind(preMedt, Medt)
   restodo<-(chisq.test(chidata))
  
 #Chisqtest for tips 
  preMedtips<-tips[tips$agestem>3.4,]
  Medtips<-tips[tips$agestem<3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)    
  chidatatips<-cbind(preMedtipst, Medtipst)
  restips<-(chisq.test(chidatatips))
  
  #Chisqtest for tips 
  
  preMednodes<-nodes[nodes$agestem>3.4,]
  Mednodes<-nodes[nodes$agestem<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)    
  chidatanodes<-cbind(preMednodest, Mednodest)
  resnodes<-(chisq.test(chidatanodes))
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
    }

#Analisis para los stem con mutaciones 1/0
setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
listfile<-dir()
resultschrom<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
    tryCatch(resultschrom[i,]<-chischrom2(file=listfile[i]), error=function(e) {
      print('Error')    })
 }

colnames(resultschrom)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom)<-listfile

write.table(resultschrom, file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_stem_binarytransitions.txt", sep="\t")


#Analisis para los stem con número de mutaciones por nodo

setwd("/home/fbalao/Datos/R/Rpackages/ChromTT/summary/")
chischrom2_Nmut<-function(file){
  x<-read.table(file, header=T, strip.white = T)
  x$n<- substr(x$name, start=1, stop=2)
  nodes<-x[grep("N[1-9]", x$n),]
  tips<-x[-(grep("N[1-9]", x$n)),]
  preMed<-x[x$agestem>3.4,]
  Med<-x[x$agestem<3.4,]
  preMedt<-table(rowSums(preMed[,5:9], na.rm = T)==0)
  if (dim(preMedt)==2){preMedt[1]<-sum(rowSums(preMed[,5:9], na.rm = T))
  }  else {preMedt<-c(sum(rowSums(preMed[,5:9], na.rm = T)), preMedt)}  
  Medt<-table(rowSums(Med[,5:9], na.rm = T)==0)    
  if (dim(Medt)==2){preMedt[1]<-sum(rowSums(Med[,5:9], na.rm = T))
  }  else {Medt<-c(sum(rowSums(Med[,5:9], na.rm = T)), Medt)}
  chidata<-cbind(preMedt, Medt)
  restodo<-(chisq.test(chidata))
  
  #Chisqtest for tips 
  preMedtips<-tips[tips$agestem>3.4,]
  Medtips<-tips[tips$agestem<3.4,]
  preMedtipst<-table(rowSums(preMedtips[,5:9], na.rm = T)==0)
  if (dim(preMedtips)==2){preMedtipst[1]<-sum(rowSums(preMedtips[,5:9], na.rm = T))
  }
  else {preMedtipst<-c(sum(rowSums(preMedtips[,5:9], na.rm = T)), preMedtipst)}
  Medtipst<-table(rowSums(Medtips[,5:9], na.rm = T)==0)   
  if (dim(Medtips)==2){Medtipst[1]<-sum(rowSums(Medtips[,5:9], na.rm = T))
  }
  else {Medtipst<-c(sum(rowSums(Medtips[,5:9], na.rm = T)), Medtipst)}
  chidatatips<-cbind(preMedtipst, Medtipst)
  restips<-(chisq.test(chidatatips))
  
  #Chisqtest for tips 
  
  preMednodes<-nodes[nodes$agestem>3.4,]
  Mednodes<-nodes[nodes$agestem<3.4,]
  preMednodest<-table(rowSums(preMednodes[,5:9], na.rm = T)==0)
  if (dim(preMednodest)==2){preMednodest[1]<-sum(rowSums(preMednodes[,5:9], na.rm = T))
  } else {preMednodest<-c(sum(rowSums(preMednodes[,5:9], na.rm = T)), preMednodest)}
  Mednodest<-table(rowSums(Mednodes[,5:9], na.rm = T)==0)  
  if (dim(Mednodest)==2){Mednodest[1]<-sum(rowSums(Mednodes[,5:9], na.rm = T))
  } else {Mednodest<-c(sum(rowSums(Mednodes[,5:9], na.rm = T)), Mednodest)}
  chidatanodes<-cbind(preMednodest, Mednodest)
  resnodes<-(chisq.test(chidatanodes))
  
  chisq<-c(c(restodo[[1]],restodo[[3]]),c(restips[[1]],restips[[3]]), c(resnodes[[1]],resnodes[[3]]))
  chisq
}

listfile<-dir()
resultschrom<-matrix(nrow=length(listfile), ncol=6)
for (i in 1:length(listfile)){
  tryCatch(resultschrom[i,]<-chischrom2_Nmut(file=listfile[i]), error=function(e) {
    print('Error')    })
}

colnames(resultschrom)<-c("Chisq_all","p-value_all","Chisq_tips","p-value_tips","Chisq_nodes","p-value_nodes")
row.names(resultschrom)<-listfile

write.table(resultschrom, file="/home/fbalao/Datos/R/Rpackages/ChromTT/results/Analysis_stem_Nmutransitions.txt", sep="\t")
