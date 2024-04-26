
#SEMILLA
set.seed(1726)
options(scipen = 999)

#Librerias
pacman::p_load(FactoMineR, stats,factoextra, ggplot2, readxl, dplyr, lavaan, ade4, scales, psych, 
               parallel, doSNOW, arules, homals, Gifi, Metrics, tidyr, tidyverse, hrbrthemes, 
               viridis, stringr, envalysis, ggthemes, tictoc, semPlot)


#-------------------------------------------------------------------------------
#ESTUDIO CON DATOS REALES
#-------------------------------------------------------------------------------

#Base de datos 

Datos_Carolina_Engegement <- read_excel("Datos_Carolina_Engegement.xlsx")



#Función para estandarizar la matriz de datos
estandarizar = function(z){
  
  media = mean(z)
  desviacion = sd(z)
  
  (z - media) / desviacion
  
}

#-------------------------------------------------------------------------------
#Construcción de Escalas con el MÉTODO 1 

i.acp.pesos <- function(x) {
  
  #dat <- as.data.frame(lapply(x ,as.numeric))
  dat <- as.data.frame(x)
  
  acp <- FactoMineR::PCA(dat, graph = FALSE)
  
  # Extraer vectores propios cuantitativos 
  VectorPropio<-round(acp$svd$V,2)  # Primer vector propio (Pesos de las variables observables)
  I<- VectorPropio[,1] 
  
  I
}


i.acp.brutos <- function(x, I) {
  
  z = as.data.frame(x)
  z = as.matrix(z)
  
  z = estandarizar(z) #Se estandariza la matriz de datos original
  
  pacp <- diag(plyr::ldply(I)[,1]) # Diagonal de pesos
  
  matriz <- as.matrix(z) # Matriz de datos estandarizada
  
  xp_acp <- matriz%*%pacp # multiplicación matriz datos * diagonal pesos
  
  iacp <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2) #Suma por filas xp_acp

  iacp #Escalas MÉTODO 1
}

#-------------------------------------------------------------------------------
#Pesos de las variables de los datos según Método 1 
acp1000 = i.acp.pesos(Datos_Carolina_Engegement) 

#Escalas construidas con el Método 1
escalas.ACP_reales_items = i.acp.brutos(Datos_Carolina_Engegement, acp1000) 

#-------------------------------------------------------------------------------

#Construcción de Escalas con el MÉTODO 3

i.afm.pesos <- function(x) {

  dat <- as.data.frame(x)
  
  afm <- FactoMineR::MFA(dat,  group = c(6,5,6), type = c("s", "s", "s") , graph = FALSE) #AFM matriz de datos
  
  # Extraer vectores propios cuantitativos 
  I <- round(afm$global.pca$svd$V[,1], 2)  # Primer vector propio (Pesos de las variables observables)  
  I
}


i.afm.brutos <- function(x, I) {
  
  
  z = as.data.frame(x)
  z = as.matrix(z)
  
  z = estandarizar(z)  #Se estandariza la matriz de datos original
  
  pacp <- diag(plyr::ldply(I)[,1]) # Diagonal de pesos
  
  matriz <- as.data.frame(z)
  matriz <- as.matrix(matriz) # Matriz de datos estandarizada
  
  xp_acp <- matriz%*%pacp # multiplicación matriz datos * diagonal pesos
  
  iafm <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2) #Suma por filas xp_acp
  
  iafm  #Escalas MÉTODO 3
}


#-------------------------------------------------------------------------------
#Pesos de las variables de los datos según Método 1 
afm1000 = i.afm.pesos(Datos_Carolina_Engegement)

#Escalas construidas con el Método 1
escalas.AFM_reales_items = i.afm.brutos(Datos_Carolina_Engegement, afm1000) 

#-------------------------------------------------------------------------------

#Construcción de Escalas con el MÉTODO 2

EscalasACP.dim <-function(datos){
  
  matriz <- as.data.frame(datos)
  matriz1 <- as.matrix(matriz[ ,1:6])      #Matriz de datos grupo 1
  matriz2 <- as.matrix(matriz[ ,7:11])     #Matriz de datos grupo 2
  matriz3 <- as.matrix(matriz[ ,12:17])   #Matriz de datos grupo 3
  result1 <- PCA(matriz1, scale.unit = TRUE, graph = FALSE)     #ACP grupo 1
  result2 <- PCA(matriz2, scale.unit = TRUE, graph = FALSE)     #ACP grupo 2
  result3 <- PCA(matriz3, scale.unit = TRUE, graph = FALSE)     #ACP grupo 3
  Vectorpropio1<- c(result1$svd$V[,1], result2$svd$V[,1], result3$svd$V[,1])    #Primer vector propio de cada grupo 
  Vectorpropio1dim <- round(Vectorpropio1, 2) #Pesos de las variables observables
  
  
  
  
  nueva_data_ACP<- list()
  
  matriz <- as.data.frame(datos)
  matriz1 <- as.matrix(as.data.frame(matriz[ ,1:6]))    #Matriz de datos grupo 1
  matriz2 <- as.matrix(as.data.frame(matriz[ ,7:11]))     #Matriz de datos grupo 2
  matriz3 <- as.matrix(as.data.frame(matriz[ ,12:17]))    #Matriz de datos grupo 3
  
  
  z1 = estandarizar(matriz1) #Estandarización grupo 1
  z2 = estandarizar(matriz2) #Estandarización grupo 2
  z3 = estandarizar(matriz3) #Estandarización grupo 3
  
  
  pesosacp1 <- diag(Vectorpropio1dim[1:6])# Diagonal de pesos # Diagonal de pesos grupo 1 
  pesosacp2 <- diag(Vectorpropio1dim[7:11])# Diagonal de pesos # Diagonal de pesos grupo 2
  pesosacp3 <- diag(Vectorpropio1dim[12:17])# Diagonal de pesos # Diagonal de pesos grupo 3 
  
  matriz1 <- as.matrix(z1)      #Matriz de datos grupo 1
  matriz2 <- as.matrix(z2)     #Matriz de datos grupo 2
  matriz3 <- as.matrix(z3)    #Matriz de datos grupo 3
  
  
  
  xp_acp1 <- matriz1%*%pesosacp1 # multiplicación matriz grupo 1 * diagonal pesos grupo 1
  xp_acp2 <- matriz2%*%pesosacp2 # multiplicación matriz grupo 2 * diagonal pesos grupo 2
  xp_acp3 <- matriz3%*%pesosacp3 # multiplicación matriz grupo 3 * diagonal pesos grupo 3
  
  xp_acp1 <- round(as.numeric(apply(xp_acp1, 1, sum)), digits = 2) #Suma por filas xp_acp1
  xp_acp2 <- round(as.numeric(apply(xp_acp2, 1, sum)), digits = 2) #Suma por filas xp_acp2
  xp_acp3 <- round(as.numeric(apply(xp_acp3, 1, sum)), digits = 2) #Suma por filas xp_acp3
  
  matriz_dim1 <- as.matrix(as.data.frame(xp_acp1))  #Escalas grupo 1
  matriz_dim2 <- as.matrix(as.data.frame(xp_acp2))  #Escalas grupo 2
  matriz_dim3 <- as.matrix(as.data.frame(xp_acp3))  #Escalas grupo 3
  
  nueva_data_ACP <- cbind(matriz_dim1, matriz_dim2,matriz_dim3) #Matriz con las 3 nuevas variables 
  #conformadas por escalas
  
  
  matriz <- nueva_data_ACP
  result <- PCA(matriz, scale.unit = TRUE, graph = FALSE) #ACP nueva base de datos
  VectorPropio1_data<- round(result$svd$V[,1],2) #Primer vector propio (Pesos)
  
  
  eacpdim<- list()
  
  
  z = as.data.frame(nueva_data_ACP) 
  z = as.matrix(z) #Nueva base de datos
  
  z = estandarizar(z) #Nueva base de datos estandarizada
  
  pacp <- diag(plyr::ldply(VectorPropio1_data)[,1])# Diagonal de pesos
  
  
  matriz <- as.matrix(z) # Matriz de datos
  
  xp_acp <- matriz%*%pacp # producto entre matriz datos * diagonal pesos
  
  eacpdim[[i]] <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2) # escalas por dimensiones
  
  
  eacpdim #Escalas MÉTODO 2
}


#Escalas construidas con el Método 2
escalas.ACP_reales_dim <- EscalasACP.dim(Datos_Carolina_Engegement)

#-------------------------------------------------------------------------------

#Construcción de Escalas con el MÉTODO 4

EscalasAFM.dim <- function(datos){
  
  matriz <- as.data.frame(datos)
  Afm <-MFA(matriz, c(6,5,6), c("s", "s", "s"), graph = FALSE)  #Realización AFM
  VectorpropioAFM_dim<- c(Afm$separate.analyses$Gr1$svd$V[,1], Afm$separate.analyses$Gr2$svd$V[,1], Afm$separate.analyses$Gr3$svd$V[,1])
  #Primer vector propio (pesos unidos)
  
  
  z = as.data.frame(datos)
  z = as.matrix(z)
  
  z = estandarizar(z)  #Matriz de datos estandarizada
  
  pesosafm <- diag(plyr::ldply(VectorpropioAFM_dim)[,1])# Diagonal de pesos
  
  
  matriz <- as.matrix(z) # Matriz de datos estandarizada
  # Matriz de datos
  
  xp_afm <- matriz%*%pesosafm #Producto entre matriz datos * diagonal pesos
  
  eafm_Separado <- round(as.numeric(apply(xp_afm, 1, sum)), digits = 2) #Suma xp_afm
  
  eafm_Separado #Escalas MÉTODO 4
}



#Escalas construidas con el Método 2
escalas.AFM_reales_dim <- EscalasAFM.dim(Datos_Carolina_Engegement)

#-------------------------------------------------------------------------------
#TRANSFORMACIÓN DE ESCALAS
#-------------------------------------------------------------------------------

#Normalizar 0 a 1000

normalizar <- function(datos) {
  data <- as.data.frame(datos)
  EscalaMinima <- min(data) #Puntaje mínimo de la escala
  EscalaMaxima <- max(data) #Puntaje máximo de la escala
  Escala <- round(((data - EscalaMinima) / (EscalaMaxima - EscalaMinima)) * 100, 2)
  Escala <- as.vector(Escala)
  return(Escala)
}

#Escalas normalizadas para cada método
escalas.ACP_reales_items_NORM = normalizar(escalas.ACP_reales_items) 
escalas.ACP_reales_dim_NORM = normalizar(escalas.ACP_reales_dim) 
escalas.AFM_reales_items_NORM = normalizar(escalas.AFM_reales_items) 
escalas.AFM_reales_dim_NORM = normalizar(escalas.AFM_reales_dim) 

#-------------------------------------------------------------------------------
#PSEUDO ERROR CUADRÁTICO MEDIO
#-------------------------------------------------------------------------------

#Se replica el mismo procedimiento utilizado para el estudio de datos Simulados
set.seed(1726)
dataReal <- as.vector(escalas.AFM_reales_dim) #Este cambia según el método a evaluar


#Se realiza remuestreo a las escalas para obtener:


#Bootstrap con una muestra de 50 individuos replicado 1000 veces
muestras50 <-list()
for (i in 1:1000) {
  muestra_bootstrap <- sample(dataReal, size = 50, replace = TRUE)
  muestras50[[i]] <- muestra_bootstrap
}

#Bootstrap con una muestra de 250 individuos replicado 1000 veces
muestras250 <-list()
for (i in 1:1000) {
  muestra_bootstrap <- sample(dataReal, size = 250, replace = TRUE)
  muestras250[[i]] <- muestra_bootstrap
}

#Bootstrap con una muestra de 500 individuos replicado 1000 veces
muestras500 <-list()
for (i in 1:1000) {
  muestra_bootstrap <- sample(dataReal, size = 500, replace = TRUE)
  muestras500[[i]] <- muestra_bootstrap
}

#Bootstrap con una muestra de 700 individuos replicado 1000 veces
muestras700 <-list()
for (i in 1:1000) {
  muestra_bootstrap <- sample(dataReal, size = 700, replace = TRUE)
  muestras700[[i]] <- muestra_bootstrap
}


# Función para calcular PECM usando bootstrap
PECM <- function(muestras50,  muestras250, muestras500, muestras700){
  #Población
  datos1<- unlist(muestras50) #Escalas n=700
  datos2<- unlist(muestras250) #Escalas n=500
  datos3<- unlist(muestras500) #Escalas n=250
  datos4<- unlist(muestras700) #Escalas n=50
  
  datos <-c(datos1, datos2, datos3, datos4)  #Unión de las escalas construidas con
  # el mismo método
  mediaPoblacional <-median(datos) #Media de todos los puntajes de las escalas unidas
  
  datosEstimador<- data.frame(muestras700) #Escalas n=700
  
  MediaEstimador= as.numeric(as.vector(apply(datosEstimador, 2, median))) #Media de cada lista de escalas
  
  sesgo_cuadrado <- (median(MediaEstimador)-mediaPoblacional)^2  #Calculo del sesgo al cuadrado n=700
  
  ECM<-0
  for (i in 1:1000) {
    ECM<-ECM+((MediaEstimador[i]-mediaPoblacional)^2/1000) + sesgo_cuadrado
  }
  ECM #PECM
}

#Resultados PECM según Método de construcción de las Escalas
PECM_Reales_ACP.items<- PECM(muestras50,  muestras250, muestras500, muestras700)
PECM_Reales_ACP.dim <- PECM(muestras50,  muestras250, muestras500, muestras700)
PECM_Reales_AFM.global <- PECM(muestras50,  muestras250, muestras500, muestras700)
PECM_Reales_AFM.sep <- PECM(muestras50,  muestras250, muestras500, muestras700)



#-------------------------------------------------------------------------------
#PUNTOS DE CORTE
#-------------------------------------------------------------------------------

# Crear una función para ordenar los resultados por mediana
ordenar_intervalo <- function(intervalo) {
  mediana <- sapply(intervalo, function(x) x["Median"])
  orden <- order(mediana)
  return(intervalo[orden])
}
PuntosCorte <- function(datos){
  
  resumen <- list()
  
  # Iterar sobre las 1000 matrices del boostrap realizado para 700 individuos
  for (i in 1:1000) {
    EscalaMinima <- min(unlist(datos[[i]])) #Puntaje mínimo de las Escalas
    EscalaMaxima <- max(unlist(datos[[i]])) #Puntaje máximo de las Escalas
    Escala <- round(((unlist(datos[[i]]) - EscalaMinima) / (EscalaMaxima - EscalaMinima)) * 100, 2) #Escala transformada
    C <- as.data.frame(Escala) #Nuevo data frame
    kmeans <- kmeans(C, 3, iter.max = 1000, nstart = 10) #tres clusters, máximo 1000 interacciones y 
    #10 intentos de inicialización aleatoria
    C$cluster <- kmeans$cluster #Añade columna cluster a C
    cluster_summary <- tapply(C$Escala, C$cluster, summary) #Resumen de los datos
    
    resumen[[i]] <- cluster_summary
  }
  
  
  # Organizar los resultados por intervalo (bajo, medio, alto)
  resultados_organizados <- lapply(resumen, ordenar_intervalo)
  mediasMinimo <-list()
  mediasMaximo <-list()
  for(i in 1:3){
    valores <-sapply(resultados_organizados, function(x) x[[i]])
    mediasMinimo[[i]]<-mean(valores["Min.",])
    mediasMaximo[[i]]<-mean(valores["Max.",])
  }
  Medias <-cbind(mediasMinimo, mediasMaximo)  # Media del mínimo y máximo de cada cluster
}

# Resumen de la media del mínimo y máximo de cada cluster de los puntajes de las escalas
#de todos los métodos
resumen_REALES_ACP_items <- PuntosCorte(muestras700)
resumen_REALES_ACP_dim <- PuntosCorte(muestras700)
resumen_REALES_AFM_items <- PuntosCorte(muestras700)
resumen_REALES_AFM_dim <- PuntosCorte(muestras700)
