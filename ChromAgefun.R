# Esta función toma los resultados de Chromevol y los resume en una tabla
# Se necesita el arbol original (Nexus tree) y los archivos *_mlAncestors.tree y *_expectations.txt 
# del mejor modelo estimado de ChromEvol
# La estructura de los archivos debe ser como estan en Drive. En la carpeta principal el tree
# en la carpeta OUT/1_Best el resto de archivos
# Tiene dos argumentos name: Nombre general del archivo (género) y path: la dirección 
# a la carpeta con los datos



#ChromR(name="Reseda",path="/home/fbalao/Datos/R/Rpackages/ChromTT/Reseda/")

ChromR<- function(name, path){
#Cargar paquetes necesarios. si no los tienes instalar con install.packages
  library(phytools)
  library(data.table)
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
names(nodeage)<-nodenames
plot(arbol)
nodelabels(nodenames)
time<-c(rep(0,ntip),nodeage)

#Obtener los números de cromosomas en cada nodo/tip *_mlAncestors.tree
path2<- paste("./OUT/1_Best_CONST_RATE_DEMI_EST/", name, "_mlAncestors.tree", sep="")
chromtree<-readLines(path2)
chromtree2<-gsub("\\[|\\]", "", chromtree)
chromnumtree<-read.tree(text=chromtree2)
chromnum<-data.frame(name=c(chromnumtree$tip.label, chromnumtree$node.label))
chromnum$chrom<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 2)
chromnum$name<-lapply(strsplit(as.character(chromnum$name), "\\-"), "[", 1)

chromtable<-data.frame(name=c(arbol$tip.label, nodenames), age=time)

chromtable2<-merge(chromtable, chromnum, by="name")
chromtable2<-chromtable2[with(chromtable2, order(-age)), ]

scatter.smooth(-chromtable2$age, chromtable2$chrom, pch=16, lwd=2, col=1)

# Leer la tabla de mutaciones *_expectations.txt con paquete data.table
path3<- paste("./OUT/1_Best_CONST_RATE_DEMI_EST/", name, "_expectations.txt", sep="")
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

chromtable3[,4:8]<-round2(chromtable3[,4:8],0)
chromtable3

#Exportar la tabla a archivo tabulado
name2<-paste(name, "_summary.txt", sep="")
write.table(as.matrix(chromtable3), file=name2, sep="\t", fileEncoding = "utf8", row.names = F)}


plot(-chromtable2$age, chromtable2$chrom, type="p")
cbind(-chromtable2$age, as.numeric(chromtable2$chrom))

#Hacer una columna si es del onset o del stem