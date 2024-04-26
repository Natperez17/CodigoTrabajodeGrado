#Histogramas 

data <- list(E_acp_alto_global_n4[[1]],E_acp_medio_global_n4[[1]],E_acp_bajo_global_n4[[1]],
             eacpdim_alto_n4[[1]],eacpdim_medio_n4[[1]],eacpdim_bajo_n4[[1]],
             E_afm_alto_global_n4[[1]],E_afm_medio_global_n4[[1]], E_afm_bajo_global_n4[[1]],
             eafmdim_alto_n4[[1]], eafmdim_medio_n4[[1]], eafmdim_bajo_n4[[1]])

data <- unlist(data)
Grupo <- rep(1:12, each = 1000)

data_final <- data.frame(data, Grupo)
nombres_grupos <- c("Método1 Escenario1", "Método1 Escenario2", "Método1 Escenario3", 
                    "Método2 Escenario1", "Método2 Escenario2", "Método2 Escenario3",
                    "Método3 Escenario1", "Método3 Escenario2","Método3 Escenario3", 
                    "Método4 Escenario1", "Método4 Escenario2", "Método4 Escenario3")

# Asignar nombres a los niveles del factor Grupo
data_final$Grupo <- factor(data_final$Grupo, levels = 1:12, labels = nombres_grupos)
colnames(data_final) <- c("Valores", "Grupo")



# Crear el histograma con facet_wrap
ggplot(data_final, aes(x = Valores)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
  facet_wrap(~ Grupo, ncol = 3) +
  labs(title = "Histograma de los puntajes de las escalas", x = "Puntajes", y = "Frecuencia")

#-------------------------------------------------------------------------------
#Gráfico del círculo de correlaciones

colores_grupo <- c("#1874CD","#1874CD","#1874CD","#1874CD","#1874CD","#1874CD","#1874CD","#1874CD", 
                   "#FF7F24", "#FF7F24","#FF7F24","#FF7F24","#FF7F24", "#FF7F24","#FF7F24", 
                   "darkgoldenrod1", "darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1", "darkgoldenrod1", "darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1")
nombres_variables <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25")

matriz <-as.data.frame(r.3.l1.alto.25.n4[[523]])
colnames(matriz) <- nombres_variables

resultACP <- PCA(matriz, scale.unit = TRUE, graph = FALSE)

fviz_pca_var(resultACP, col.var = colores_grupo,
             title="Círculo de correlaciones",
             col.ind = "red", legend.title=("Grupos"))


#-------------------------------------------------------------------------------
#Análisis factorial confirmatorio

colnames(Datos_Carolina_Engegement) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25")
modelo <- 
  " 
   Vigor =~ X1+ X2+ X3+ X4+ X5+ X6
   Dedicacion=~ X7+ X8+ X9+ X10+ X11
   Absorcion =~ X12+ X13+ X14+ X15+ X16+ X17
"
afc <- lavaan::cfa(modelo, data = Datos_Carolina_Engegement)
summary(afc)
semPlot::semPaths(afc, nCharNodes = 0,intercepts = FALSE, edge.label.cex=1.3, optimizeLatRes = T, groups = "lat",pastel = T, sizeInt=5,edge.color ="black",esize = 5, label.prop=0,sizeLat = 11,"std",layout="circle3", exoVar = F)
