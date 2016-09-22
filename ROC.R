library("pROC")
library("caTools")


ngenes = 126;
ngenes_2 = ngenes*ngenes;

file_read = paste("Edges/",tipo,"/Network",net,"/edges.txt", sep = "")
gold_standard = paste("gold_standard/gold_standard_126",sep="")
gold_network = as.matrix(read.table(gold_standard, header = FALSE, sep = "\t"))
test_network = as.matrix(read.table(file_read, header = FALSE, sep = "\t"))

verdadeiro_positivo = 0;
verdadeiro_negativo = 0;
falso_positivo = 0;
falso_negativo = 0;
test_positivo =  0;
test_negativo =  0;
gold_positivo =  0;
gold_negativo =  0;

teste.matriz <- matrix(data=0,nrow=ngenes,ncol=ngenes)
gold.matriz <- matrix(data=0,nrow=ngenes,ncol=ngenes)
gold.matriz[126,126]
#######################333
resposta = vector(mode =  "integer", length =  ngenes_2)
teste = vector(mode =  "integer", length =  ngenes_2)

aupr_X_precision =vector(mode =  "double", length =  ngenes_2)
aupr_Y_recall = vector(mode =  "double", length =  ngenes_2)
auroc_X_tru = vector(mode =  "double", length =  ngenes_2)
auroc_Y_false = vector(mode =  "double", length =  ngenes_2)

length(resposta) 
length(teste)
i = 1;
while(i <= nrow(gold_network) || i <= nrow(test_network))
{
  if(i <= nrow(gold_network))
  {
    #if((gold_network[i,1]) != (gold_network[i,2]))
     # gold_positivo = gold_positivo + 1
    resposta[((gold_network[i,1])*ngenes)+gold_network[i,2]] = 1
    resposta[((gold_network[i,2])*ngenes)+gold_network[i,1]] = 1
    gold.matriz[(gold_network[i,1])+1,(gold_network[i,2])+1] = 1
    gold.matriz[(gold_network[i,2])+1,(gold_network[i,1])+1] = 1
  }
  if(i <= nrow(test_network)){
    teste[((test_network[i,1])*ngenes)+test_network[i,2]] = 1
    teste[((test_network[i,2])*ngenes)+test_network[i,1]] = 1
    teste.matriz[(test_network[i,1])+1,(test_network[i,2])+1] = 1
    teste.matriz[(test_network[i,2])+1,(test_network[i,1])+1] = 1
  }
  i = i + 1
}


for(i in 1:ngenes)
{
  for(j in 1:ngenes)
  {
    if(gold.matriz[i,j] == 0)
    {
      if(teste.matriz[i,j] == 0)
        verdadeiro_negativo = verdadeiro_negativo + 1
      
      else
        falso_positivo = falso_positivo + 1
      
    }
    else
    {
      if(teste.matriz[i,j] == 1)
        verdadeiro_positivo = verdadeiro_positivo + 1
      else
        falso_negativo = falso_negativo + 1
      
      
    }
  }
  
 # taxa_verdadeiro_positivo = as.double(verdadeiro_positivo/gold_positivo)
  #taxa_falso_positivo = as.double(falso_positivo/gold_negativo)
  #precision = as.double(verdadeiro_positivo/(verdadeiro_positivo+falso_positivo))
  #recall = as.double(taxa_verdadeiro_positivo)
  #aupr_X_precision[i] = precision
  #aupr_Y_recall[i]= recall
  #auroc_X_tru[i] = taxa_verdadeiro_positivo
  #auroc_Y_false[i] = taxa_falso_positivo
}

test_positivo =  verdadeiro_positivo+falso_positivo;
test_negativo =  falso_negativo+verdadeiro_negativo;
gold_positivo =  verdadeiro_positivo+falso_negativo;
gold_negativo =  verdadeiro_negativo+falso_positivo;

verdadeiro_positivo = 0;
verdadeiro_negativo = 0;
falso_positivo = 0;
falso_negativo = 0;

for(i in 1:ngenes_2)
{
  if(resposta[i] == 1)
  {
    if(teste[i] == 1)
      verdadeiro_positivo = verdadeiro_positivo + 1
    else
      falso_negativo = falso_negativo + 1
  }
  else
  {
    if(teste[i] == 0)
      verdadeiro_negativo = verdadeiro_negativo + 1
    else
      falso_positivo = falso_positivo + 1
  }
  taxa_verdadeiro_positivo = verdadeiro_positivo/gold_positivo
  taxa_falso_positivo = falso_positivo/gold_negativo
  precision = verdadeiro_positivo/(verdadeiro_positivo+falso_positivo)
  recall = taxa_verdadeiro_positivo
  aupr_X_precision[i] = precision
  aupr_Y_recall[i]= recall
  auroc_X_tru[i] = taxa_verdadeiro_positivo
  auroc_Y_false[i] = taxa_falso_positivo
  
}

especificity = verdadeiro_negativo/(verdadeiro_negativo+falso_positivo)
PPV = verdadeiro_positivo/(verdadeiro_positivo+falso_positivo)
similarity = sqrt(PPV * especificity)
print("verdadeiro_negativo")
print(verdadeiro_negativo)
print("verdadeiro_positivo")
print(verdadeiro_positivo)
print("falso_negativo")
print(falso_negativo)
print("falso_positivo")
print(falso_positivo)
print("especificity")
print(especificity)
print("PPV")
print(PPV)
print("similarity")
print(similarity)


AUROC = trapz(auroc_Y_false, auroc_X_tru)
AUPR = trapz(aupr_Y_recall,aupr_X_precision) 


print("verdadeiro_positivo+falso_negativo")
print(verdadeiro_positivo+falso_negativo)
print("gold_positivo")
print(gold_positivo)
print("falso_positivo+verdadeiro_negativo")
print(falso_positivo+verdadeiro_negativo)
print("gold_negativo")
print(gold_negativo)
print("AUROC")
print(AUROC)
print("AUPR")
print(AUPR)

roc1 = roc(response = resposta, predictor = teste)
area = auc(roc1)
relatorio = paste("AUC de", tipo, "=", area, sep = " ")
cria_diretorio("AUC")
file_write = paste("AUC/Relatorio_Network",net,".txt", sep = "")
write(relatorio, file = file_write, append = TRUE)
#plot.roc(roc1)
