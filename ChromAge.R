# Este script toma los resultados de Chromevol y los resume en una tabla
# Se necesita el arbol original (Nexus tree) y los archivos *_mlAncestors.tree y *_expectations.txt 
# del mejor modelo estimado de ChromEvol

#Cargar paquetes necesarios. si no los tienes instalar con install.packages
library(phytools)
library(data.table)

#Fija el directorio
setwd("./Hedera/")
#lee el arbol original para obtener las edades de los nodos
tree<-"Hedera.tree"
arbol<-read.nexus(tree)
ntip<-length(arbol$tip.label)
nnodes<-length(node.height(arbol))
ninodes<-nnodes-ntip
nodenames<-paste("N", 1:ninodes, sep="")
nodeage<- branching.times(arbol)
names(nodeage)<-nodenames
plot(arbol)
nodelabels(nodenames)
time<-c(rep(0,ntip),nodeage)

#Obtener la edad del stem asociado a cada nodo
stem1<-arbol$edge
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
agestem<-cbind(as.numeric.factor(revalue(factor(stem1[,1]), nodeage)), stem1[,2])
rootnode<-ntip+1
rootage<-nodeage[as.character(rootnode)]
agestem<-rbind(agestem, c(rootage,rootnode))
agestem<-agestem[order(agestem[,2]),]


#Obtener los números de cromosomas en cada nodo/tip *_mlAncestors.tree
chromtree<-readLines("./OUT/1_Best_BASE_NUM/Hedera_mlAncestors.tree")
chromtree2<-gsub("\\[|\\]", "", chromtree)
chromnumtree<-read.tree(text=chromtree2)
chromnum<-data.frame(name=c(chromnumtree$tip.label, chromnumtree$node.label))
chromnum$chrom<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 2)
chromnum$name<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 1)

chromtable<-data.frame(name=c(arbol$tip.label, nodenames), age=time, agestem=agestem[,1])

chromtable2<-merge(chromtable, chromnum, by="name")
chromtable2<-chromtable2[with(chromtable2, order(-age)), ]

plot(-chromtable2$age, chromtable2$chrom, type="p")

# Leer la tabla de mutaciones *_expectations.txt con paquete data.table
mutations <- fread("./OUT/1_Best_BASE_NUM/Hedera_expectations.txt", 
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

chromtable3[,4:8]<-round2(chromtable3[,4:8],0)
chromtable3

#Exportar la tabla a archivo tabulado
write.table(as.matrix(chromtable3), file="Hedera_summary.csv",sep="\t", fileEncoding = "utf8", row.names = F)
