# Script de seleçao dos preditores. Seleciona-se o de maior score e em caso de empate escolhe o de menor números de preditorese.
# input = arquivos de preditores para k = 1, 2, 3 e 4
# output = 1 arquivo com os preditores selecionados.

library(igraph)
net = "_completo"
tipo = "IM_pen_zero_SDB_126"
file_read = paste("Edges/",tipo,"/Network",net,"/edges.txt", sep = "")
file_read = paste("gold_standard/gold_standard_126", sep = "")
preditors <- as.data.frame(read.table(file_read, header=FALSE, sep = "\t",as.is=TRUE))
preditors.network<-graph.data.frame(preditors, directed=T)
V(preditors.network)$color <- "orange"
tkplot(preditors.network)

g2  <-graph.data.frame(preditors, directed=T)
g <- simplify(g2)
V(g)$color <- "orange"

tkplot(g)
## V(preditors.network) #prints the list of vertices (people)
## E(preditors.network) #prints the list of edges (relationships)
## degree(preditors.network) #print the number of edges per vertex (relationships per people)

## bad.vs<-V(preditors.network)[degree(preditors.network)<3] #identify those vertices part of less than three edges
## preditors.network<-delete.vertices(preditors.network, bad.vs) #exclude them from the graph

#write.csv(preditors, file = "Área de Trabalho/Projeto R/preditors_selected/preditors.csv")

