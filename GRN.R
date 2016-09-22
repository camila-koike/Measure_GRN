library("Rcpp")
sourceCpp("IFFS.cpp")

cria_diretorio<-function (diretorio)
{
  if(dir.exists(diretorio) == FALSE)
  {
    dir.create(diretorio)
  }
}

InputSDB = paste("Dados_String/input_string_intersec_", tam, "_genes_normalizado.txt", sep = "")
InputExpr = paste("Dados_Microarray/input_filtro_intersec_spell_kegg_",tam, net, ".txt", sep = "")
InputCS = paste("chip_set/chip_set", net, ".txt", sep = "")


cria_diretorio("Relatorios de tempo")
pasta = paste("Relatorios de tempo/Network", net,sep = "")
cria_diretorio(pasta)

time_report = paste("Relatorios de tempo/Network", net,"/Relatorio de tempo ", tipo, sep = "")
unlink(time_report)

tempo_total_ini = Sys.time()
for(i in 1:4)
{
  a = Sys.time()
  max_k = i
  p = paste("Network",net, " K = ", i, sep = "")
  print(p)
  ################### criando as pastas #############
  cria_diretorio("Resultado_com_Score")
  pasta = paste("Resultado_com_Score/",tipo,sep = "")
  cria_diretorio(pasta)
  pasta = paste(pasta,"/Network",net,sep = "")
  cria_diretorio(pasta)
  cria_diretorio("Preditores_dos_Alvos")
  pasta = paste("Preditores_dos_Alvos/",tipo,sep = "")
  cria_diretorio(pasta)
  pasta = paste(pasta,"/Network",net,sep = "")
  cria_diretorio(pasta)
  pasta = paste(pasta,"/k",max_k,sep = "")
  cria_diretorio(pasta)
  #########################################
  resultado_com_score = paste("Resultado_com_Score/", tipo,"/Network",net,"/Network",net,"_preditors",max_k,".txt", sep = "")
  output = paste("Preditores_dos_Alvos/", tipo,"/Network",net,"/k",max_k,"/alvo", sep = "")
  RunIffs(InputExpr, InputSDB, InputCS, output, resultado_com_score, max_k, delta, w_cod, w_cod_rara, w_cod_zero, w_im,
         w_im_rara, w_im_zero, w_tsallis, w_sdb)
  b = Sys.time()
  text =as.character.POSIXt(difftime(b,a))
  print(text)
  relatorio = paste("Network", net, "tipo de função J", tipo, "K =", i, "tempo de", text, sep = " ")
  write(relatorio,time_report, append = TRUE)
}
tempo_total_fim = Sys.time()
text =as.character.POSIXt(difftime(tempo_total_fim,tempo_total_ini))
relatorio = paste("TOTAL = ", text, sep = " ")
write(relatorio,time_report, append = TRUE)
