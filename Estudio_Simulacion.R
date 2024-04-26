
#SEMILLA
set.seed(1726)
options(scipen = 999)

#Librerias
pacman::p_load(FactoMineR, stats,factoextra, ggplot2, readxl, dplyr, lavaan, ade4, scales, psych, 
               parallel, doSNOW, arules, homals, Gifi, Metrics, tidyr, tidyverse, hrbrthemes, 
               viridis, stringr, envalysis, ggthemes, tictoc, semPlot)


#-------------------------------------------------------------------------------
#ESTUDIO CON DATOS SIMULADOS
#-------------------------------------------------------------------------------

#Base de datos
r.e.l1.alto.25.n1 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.alto.25.n1.rds")
r.e.l1.alto.25.n2 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.alto.25.n2.rds")
r.e.l1.alto.25.n3 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.alto.25.n3.rds")
r.e.l1.alto.25.n4 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.alto.25.n4.rds")

r.e.l1.medio.25.n1 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.medio.25.n1.rds")
r.e.l1.medio.25.n2 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.medio.25.n2.rds")
r.e.l1.medio.25.n3 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.medio.25.n3.rds")
r.e.l1.medio.25.n4 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.medio.25.n4.rds")

r.e.l1.bajo.25.n1 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.bajo.25.n1.rds")
r.e.l1.bajo.25.n2 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.bajo.25.n2.rds")
r.e.l1.bajo.25.n3 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.bajo.25.n3.rds")
r.e.l1.bajo.25.n4 <- readRDS("C:/Users/natal/OneDrive - correounivalle.edu.co/Trabajo de Grado/r.e.l1.bajo.25.n4.rds")

#Función para estandarizar la matriz de datos

estandarizar = function(z){
  
  media = mean(z)
  desviacion = sd(z)
  
  (z - media) / desviacion
  
}

#-------------------------------------------------------------------------------
#Construcción de Escalas con el MÉTODO 1 
i.acp.pesos <- function(x) {
  
  I <- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for(i in 1:1000){
    dat <- as.data.frame(x[[i]])
    
    acp <- FactoMineR::PCA(dat, graph = FALSE)
    
    # Extraer vectores propios cuantitativos 
    VectorPropio<-round(acp$svd$V,2)  #Vectores propios
    I[[i]]<- VectorPropio[,1]  # Primer vector propio (Pesos de las variables observables)
  } 
  I
}


i.acp.brutos <- function(x, I) {
  
  iacp <- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for(i in 1:1000){
    z = as.data.frame(x[[i]])
    z = as.matrix(z)
    
    z = estandarizar(z) #Se estandariza la matriz de datos original
    
    pacp <- diag(plyr::ldply(I[[i]])[,1]) # Diagonal de pesos
    
    matriz <- as.matrix(z) # Matriz de datos estandarizada
    
    xp_acp <- matriz%*%pacp # multiplicación matriz datos * diagonal pesos
    
    iacp[[i]] <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2) #Suma por filas xp_acp
    
  } 
  iacp #Escalas MÉTODO 2
}
#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario alto (Escenario 1) y Método 1 

acp1000 = i.acp.pesos(r.3.l1.alto.25.n4) #n=1000
acp500 = i.acp.pesos(r.3.l1.alto.25.n3)  #n=500
acp250 = i.acp.pesos(r.3.l1.alto.25.n2)  #n=250
acp50 = i.acp.pesos(r.3.l1.alto.25.n1)   #n=50


#Escalas para las 1000 matrices del escenario alto (Escenario 1) y Método 1 
E_acp_alto_global_n4 = i.acp.brutos(r.3.l1.alto.25.n4, acp1000) #n=1000
E_acp_alto_global_n3 = i.acp.brutos(r.3.l1.alto.25.n3, acp500)  #n=500
E_acp_alto_global_n2 = i.acp.brutos(r.3.l1.alto.25.n2, acp250)  #n=250
E_acp_alto_global_n1 = i.acp.brutos(r.3.l1.alto.25.n1, acp50)   #n=50

#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario medio (Escenario 2) y Método 1 
acp_medio1000 = i.acp.pesos(r.3.l1.medio.25.n4) #n=1000
acp_medio500 = i.acp.pesos(r.3.l1.medio.25.n3)  #n=500
acp_medio250 = i.acp.pesos(r.3.l1.medio.25.n2)  #n=250
acp_medio50 = i.acp.pesos(r.3.l1.medio.25.n1)   #n=50

#Escalas para las 1000 matrices del escenario medio (Escenario 2) y Método 1  
E_acp_medio_global_n4 = i.acp.brutos(r.3.l1.medio.25.n4, acp_medio1000) #n=1000
E_acp_medio_global_n3 = i.acp.brutos(r.3.l1.medio.25.n3, acp_medio500)  #n=500
E_acp_medio_global_n2 = i.acp.brutos(r.3.l1.medio.25.n2, acp_medio250)  #n=250
E_acp_medio_global_n1 = i.acp.brutos(r.3.l1.medio.25.n1, acp_medio50)   #n=50

#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario bajo (Escenario 3) y Método 1 

acp_bajo1000 = i.acp.pesos(r.3.l1.bajo.25.n4) #n=1000
acp_bajo500 = i.acp.pesos(r.3.l1.bajo.25.n3)  #n=500
acp_bajo250 = i.acp.pesos(r.3.l1.bajo.25.n2)  #n=250
acp_bajo50 = i.acp.pesos(r.3.l1.bajo.25.n1)   #n=50

#Escalas para las 1000 matrices del escenario bajo (Escenario 3) y Método 1 
E_acp_bajo_global_n4 = i.acp.brutos(r.3.l1.bajo.25.n4, acp_bajo1000) #n=1000
E_acp_bajo_global_n3 = i.acp.brutos(r.3.l1.bajo.25.n3, acp_bajo500)  #n=500
E_acp_bajo_global_n2 = i.acp.brutos(r.3.l1.bajo.25.n2, acp_bajo250)  #n=250
E_acp_bajo_global_n1 = i.acp.brutos(r.3.l1.bajo.25.n1, acp_bajo50)   #n=50

#-------------------------------------------------------------------------------

#Construcción de Escalas con el MÉTODO 3
i.afm.pesos <- function(x) {
  
  I <- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for(i in 1:1000){
    
    dat <- as.data.frame(x[[i]])
    
    afm <- FactoMineR::MFA(dat,  group = c(8,7,10), type = c("s", "s", "s") , graph = FALSE) #AFM matriz de datos
    
    # Extraer vectores propios cuantitativos 
    I[[i]] <- round(afm$global.pca$svd$V[,1], 2)  # Primer vector propio (Pesos de las variables observables) 
  } 
  I
}


i.afm.brutos <- function(x, I) {
  
  iafm <- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for(i in 1:1000){
    
    z = as.data.frame(x[[i]])
    z = as.matrix(z)
    
    z = estandarizar(z) #Se estandariza la matriz de datos original
    
    pacp <- diag(plyr::ldply(I[[i]])[,1]) # Diagonal de pesos
    
    
    matriz <- as.matrix(z) # Matriz de datos estandarizada
    
    xp_acp <- matriz%*%pacp # multiplicación matriz datos * diagonal pesos
    
    iafm[[i]] <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2)  #Suma por filas xp_acp
  } 
  iafm c
}

#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario alto (Escenario 1) y Método 3 

afm1000 = i.afm.pesos(r.3.l1.alto.25.n4) 
afm500 = i.afm.pesos(r.3.l1.alto.25.n3) 
afm250 = i.afm.pesos(r.3.l1.alto.25.n2) 
afm50 = i.afm.pesos(r.3.l1.alto.25.n1) 

#Escalas para las 1000 matrices del escenario alto (Escenario 1) y Método 3
E_afm_alto_global_n4 = i.afm.brutos(r.3.l1.alto.25.n4, afm1000) 
E_afm_alto_global_n3 = i.afm.brutos(r.3.l1.alto.25.n3, afm500) 
E_afm_alto_global_n2 = i.afm.brutos(r.3.l1.alto.25.n2, afm250) 
E_afm_alto_global_n1 = i.afm.brutos(r.3.l1.alto.25.n1, afm50) 

#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario medio (Escenario 2) y Método 3 

afm_medio1000 = i.afm.pesos(r.3.l1.medio.25.n4) 
afm_medio500 = i.afm.pesos(r.3.l1.medio.25.n3) 
afm_medio250 = i.afm.pesos(r.3.l1.medio.25.n2) 
afm_medio50 = i.afm.pesos(r.3.l1.medio.25.n1) 

#Escalas para las 1000 matrices del escenario medio (Escenario 2) y Método 3

E_afm_medio_global_n4 = i.afm.brutos(r.3.l1.medio.25.n4, afm_medio1000) 
E_afm_medio_global_n3 = i.afm.brutos(r.3.l1.medio.25.n3, afm_medio500) 
E_afm_medio_global_n2 = i.afm.brutos(r.3.l1.medio.25.n2, afm_medio250) 
E_afm_medio_global_n1 = i.afm.brutos(r.3.l1.medio.25.n1, afm_medio50) 

#-------------------------------------------------------------------------------
#Pesos para las 1000 matrices del escenario alto (Escenario 3) y Método 3 

afm_bajo1000 = i.afm.pesos(r.3.l1.bajo.25.n4) 
afm_bajo500 = i.afm.pesos(r.3.l1.bajo.25.n3) 
afm_bajo250 = i.afm.pesos(r.3.l1.bajo.25.n2) 
afm_bajo50 = i.afm.pesos(r.3.l1.bajo.25.n1) 

#Escalas para las 1000 matrices del escenario medio (Escenario 3) y Método 3

E_afm_bajo_global_n4 = i.afm.brutos(r.3.l1.bajo.25.n4, afm_bajo1000) 
E_afm_bajo_global_n3 = i.afm.brutos(r.3.l1.bajo.25.n3, afm_bajo500) 
E_afm_bajo_global_n2 = i.afm.brutos(r.3.l1.bajo.25.n2, afm_bajo250) 
E_afm_bajo_global_n1 = i.afm.brutos(r.3.l1.bajo.25.n1, afm_bajo50) 


#-------------------------------------------------------------------------------
#Construcción de Escalas con el MÉTODO 2

EscalasACP.dim <-function(datos){
  
  Vectorpropio1dim<- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for (i in 1:1000){
    matriz <- as.data.frame(datos[[i]])
    matriz1 <- as.matrix(matriz[ ,1:8])      #Matriz de datos grupo 1
    matriz2 <- as.matrix(matriz[ ,9:15])     #Matriz de datos grupo 2
    matriz3 <- as.matrix(matriz[ ,16:25])    #Matriz de datos grupo 3
    result1 <- PCA(matriz1, scale.unit = TRUE, graph = FALSE)     #ACP grupo 1
    result2 <- PCA(matriz2, scale.unit = TRUE, graph = FALSE)     #ACP grupo 2
    result3 <- PCA(matriz3, scale.unit = TRUE, graph = FALSE)     #ACP grupo 3
    Vectorpropio1<- c(result1$svd$V[,1], result2$svd$V[,1], result3$svd$V[,1])    #Primer vector propio de cada grupo 
    Vectorpropio1dim[[i]] <- round(Vectorpropio1, 2) #Pesos de las variables observables
  }
  
  
  nueva_data_ACP<- list()
  #Se realiza un ciclo for para que itere sobre las 1000 matrices por cada escenario
  for (i in 1:1000){
    
    matriz <- as.data.frame(datos[[i]])
    matriz1 <- as.matrix(as.data.frame(matriz[ ,1:8]))    #Matriz de datos grupo 1
    matriz2 <- as.matrix(as.data.frame(matriz[ ,9:15]))     #Matriz de datos grupo 2
    matriz3 <- as.matrix(as.data.frame(matriz[ ,16:25]))    #Matriz de datos grupo 3
    
    
    z1 = estandarizar(matriz1) #Estandarización grupo 1
    z2 = estandarizar(matriz2) #Estandarización grupo 2
    z3 = estandarizar(matriz3) #Estandarización grupo 3
    
    
    pesosacp1 <- diag(Vectorpropio1dim[[i]][1:8]) # Diagonal de pesos grupo 1 
    pesosacp2 <- diag(Vectorpropio1dim[[i]][9:15]) #Diagonal de pesos grupo 2
    pesosacp3 <- diag(Vectorpropio1dim[[i]][16:25]) #Diagonal de pesos grupo 3 
    
    matriz1 <- as.matrix(z1)      #Matriz de datos grupo 1 estandarizada
    matriz2 <- as.matrix(z2)     #Matriz de datos grupo 2 estandarizada
    matriz3 <- as.matrix(z3)    #Matriz de datos grupo 3 estandarizada
    
    
    
    xp_acp1 <- matriz1%*%pesosacp1 # multiplicación matriz grupo 1 * diagonal pesos grupo 1
    xp_acp2 <- matriz2%*%pesosacp2 # multiplicación matriz grupo 2 * diagonal pesos grupo 2
    xp_acp3 <- matriz3%*%pesosacp3 # multiplicación matriz grupo 3 * diagonal pesos grupo 3
    
    xp_acp1 <- round(as.numeric(apply(xp_acp1, 1, sum)), digits = 2)  #Suma por filas xp_acp1
    xp_acp2 <- round(as.numeric(apply(xp_acp2, 1, sum)), digits = 2)  #Suma por filas xp_acp2
    xp_acp3 <- round(as.numeric(apply(xp_acp3, 1, sum)), digits = 2)  #Suma por filas xp_acp
    
    matriz_dim1 <- as.matrix(as.data.frame(xp_acp1))  #Escalas grupo 1
    matriz_dim2 <- as.matrix(as.data.frame(xp_acp2))  #Escalas grupo 2
    matriz_dim3 <- as.matrix(as.data.frame(xp_acp3))  #Escalas grupo 3
    
    nueva_data_ACP[[i]] <- cbind(matriz_dim1, matriz_dim2,matriz_dim3) #Matriz con las 3 nuevas variables 
    #conformadas por escalas
  }
  
  
  VectorPropio1_data<- list()
  for (i in 1:1000){
    matriz <- nueva_data_ACP[[i]]
    result <- PCA(matriz, scale.unit = TRUE, graph = FALSE) #ACP nueva base de datos
    VectorPropio1_data[[i]]<- round(result$svd$V[,1],2)  #primer vector propio (Pesos)
  }
  
  eacpdim<- list()
  for (i in 1:1000){
    
    z = as.data.frame(nueva_data_ACP[[i]])
    z = as.matrix(z) #Nueva base de datos
    
    z = estandarizar(z) #Nueva base de datos estandarizada
    
    pacp <- diag(plyr::ldply(VectorPropio1_data[[i]])[,1]) # Diagonal de pesos
    
    
    matriz <- as.matrix(z) # Matriz de datos
    
    xp_acp <- matriz%*%pacp # Producto entre matriz datos * diagonal pesos
    
    eacpdim[[i]] <- round(as.numeric(apply(xp_acp, 1, sum)), digits = 2) #Suma 
  }
  
  eacpdim #Escalas MÉTODO 2
}

#-------------------------------------------------------------------------------
#Escalas para las 1000 matrices del escenario alto (Escenario 1) y Método 2

eacpdim_alto_n4 <- EscalasACP.dim(r.3.l1.alto.25.n4)
eacpdim_alto_n3 <- EscalasACP.dim(r.3.l1.alto.25.n3)
eacpdim_alto_n2 <- EscalasACP.dim(r.3.l1.alto.25.n2)
eacpdim_alto_n1 <- EscalasACP.dim(r.3.l1.alto.25.n1)

#Escalas para las 1000 matrices del escenario medio (Escenario 2) y Método 2

eacpdim_medio_n4 <- EscalasACP.dim(r.3.l1.medio.25.n4)
eacpdim_medio_n3 <- EscalasACP.dim(r.3.l1.medio.25.n3)
eacpdim_medio_n2 <- EscalasACP.dim(r.3.l1.medio.25.n2)
eacpdim_medio_n1 <- EscalasACP.dim(r.3.l1.medio.25.n1)

#Escalas para las 1000 matrices del escenario medio (Escenario 3) y Método 2

eacpdim_bajo_n4 <- EscalasACP.dim(r.3.l1.bajo.25.n4)
eacpdim_bajo_n3 <- EscalasACP.dim(r.3.l1.bajo.25.n3)
eacpdim_bajo_n2 <- EscalasACP.dim(r.3.l1.bajo.25.n2)
eacpdim_bajo_n1 <- EscalasACP.dim(r.3.l1.bajo.25.n1)

#-------------------------------------------------------------------------------
#Construcción de Escalas con el MÉTODO 4

EscalasAFM.dim <- function(datos){
  
  VectorpropioAFM_dim <- list()
  for (i in 1:1000){
    matriz <- as.data.frame(datos[[i]])
    Afm <-MFA(matriz, c(8, 7, 10), c("s", "s", "s"), graph = FALSE)  #Realización AFM
    VectorpropioAFM_dim[[i]]<- c(Afm$separate.analyses$Gr1$svd$V[,1], Afm$separate.analyses$Gr2$svd$V[,1], Afm$separate.analyses$Gr3$svd$V[,1])
    #Primer vector propio (pesos unidos)
    
  }
  eafm_Separado <- list()
  for (i in 1:1000){
    
    z = as.data.frame(datos[[i]])
    z = as.matrix(z)
    
    z = estandarizar(z) #Matriz de datos estandarizada
    
    pesosafm <- diag(plyr::ldply(VectorpropioAFM_dim[[i]])[,1])# Diagonal de pesos
    
    
    matriz <- as.matrix(z) # Matriz de datos estandarizada
    
    xp_afm <- matriz%*%pesosafm  #Producto entre matriz datos * diagonal pesos
    
    eafm_Separado[[i]] <- round(as.numeric(apply(xp_afm, 1, sum)), digits = 2) #Suma xp_afm
  }
  eafm_Separado #Escalas MÉTODO 4
}


#-------------------------------------------------------------------------------
#Escalas para las 1000 matrices del escenario alto (Escenario 1) y Método 4

eafmdim_alto_n4 <- EscalasAFM.dim(r.3.l1.alto.25.n4)
eafmdim_alto_n3 <- EscalasAFM.dim(r.3.l1.alto.25.n3)
eafmdim_alto_n2 <- EscalasAFM.dim(r.3.l1.alto.25.n2)
eafmdim_alto_n1 <- EscalasAFM.dim(r.3.l1.alto.25.n1)

#Escalas para las 1000 matrices del escenario medio (Escenario 2) y Método 4 

eafmdim_medio_n4 <- EscalasAFM.dim(r.3.l1.medio.25.n4)
eafmdim_medio_n3 <- EscalasAFM.dim(r.3.l1.medio.25.n3)
eafmdim_medio_n2 <- EscalasAFM.dim(r.3.l1.medio.25.n2)
eafmdim_medio_n1 <- EscalasAFM.dim(r.3.l1.medio.25.n1)

#Escalas para las 1000 matrices del escenario bajo (Escenario 3) y Método 4

eafmdim_bajo_n4 <- EscalasAFM.dim(r.3.l1.bajo.25.n4)
eafmdim_bajo_n3 <- EscalasAFM.dim(r.3.l1.bajo.25.n3)
eafmdim_bajo_n2 <- EscalasAFM.dim(r.3.l1.bajo.25.n2)
eafmdim_bajo_n1 <- EscalasAFM.dim(r.3.l1.bajo.25.n1)

#-------------------------------------------------------------------------------
#TRANSFORMACIÓN DE ESCALAS
#-------------------------------------------------------------------------------

#Normalizar 0 a 1000

normalizar <- function(datos){
  norm <- list()
  
  # Iterar sobre las 1000 matrices
  for (i in 1:1000) {
    data <- as.data.frame(datos[[i]])
    EscalaMinima <- min(data) #Puntaje mínimo de la escala
    EscalaMaxima <- max(data) #Puntaje máximo de la escala
    Escala <- round(((data - EscalaMinima) / (EscalaMaxima - EscalaMinima)) * 100, 2)
    Escala <- as.vector(Escala)
    norm[[i]] <- Escala
  }
  norm
}

#Escalas normalizadas escenario alto (Escenario 1) MÉTODO 1
E_acp_alto_global_n4_NORM = normalizar(E_acp_alto_global_n4) 
E_acp_alto_global_n3_NORM = normalizar(E_acp_alto_global_n3) 
E_acp_alto_global_n2_NORM = normalizar(E_acp_alto_global_n2) 
E_acp_alto_global_n1_NORM = normalizar(E_acp_alto_global_n1)

#Escalas normalizadas escenario medio (Escenario 2) MÉTODO 1
E_acp_medio_global_n4_NORM = normalizar(E_acp_medio_global_n4) 
E_acp_medio_global_n3_NORM = normalizar(E_acp_medio_global_n3) 
E_acp_medio_global_n2_NORM = normalizar(E_acp_medio_global_n2) 
E_acp_medio_global_n1_NORM = normalizar(E_acp_medio_global_n1) 

#Escalas normalizadas escenario bajo (Escenario 3) MÉTODO 1
E_acp_bajo_global_n4_NORM = normalizar(E_acp_bajo_global_n4)
E_acp_bajo_global_n3_NORM = normalizar(E_acp_bajo_global_n3)
E_acp_bajo_global_n2_NORM = normalizar(E_acp_bajo_global_n2)
E_acp_bajo_global_n1_NORM = normalizar(E_acp_bajo_global_n1)


################################################################################

#Escalas normalizadas escenario alto (Escenario 1) MÉTODO 3
E_afm_alto_global_n4_NORM = normalizar(E_afm_alto_global_n4) 
E_afm_alto_global_n3_NORM = normalizar(E_afm_alto_global_n3) 
E_afm_alto_global_n2_NORM = normalizar(E_afm_alto_global_n2) 
E_afm_alto_global_n1_NORM = normalizar(E_afm_alto_global_n1)

#Escalas normalizadas escenario medio (Escenario 2) MÉTODO 3
E_afm_medio_global_n4_NORM = normalizar(E_afm_medio_global_n4) 
E_afm_medio_global_n3_NORM = normalizar(E_afm_medio_global_n3) 
E_afm_medio_global_n2_NORM = normalizar(E_afm_medio_global_n2) 
E_afm_medio_global_n1_NORM = normalizar(E_afm_medio_global_n1) 

#Escalas normalizadas escenario bajo (Escenario 3) MÉTODO 3
E_afm_bajo_global_n4_NORM = normalizar(E_afm_bajo_global_n4)
E_afm_bajo_global_n3_NORM = normalizar(E_afm_bajo_global_n3)
E_afm_bajo_global_n2_NORM = normalizar(E_afm_bajo_global_n2)
E_afm_bajo_global_n1_NORM = normalizar(E_afm_bajo_global_n1)

################################################################################

#Escalas normalizadas escenario alto (Escenario 1) MÉTODO 2
E_acp_alto_dim_n4_NORM = normalizar(eacpdim_alto_n4) 
E_acp_alto_dim_n3_NORM = normalizar(eacpdim_alto_n3) 
E_acp_alto_dim_n2_NORM = normalizar(eacpdim_alto_n2) 
E_acp_alto_dim_n1_NORM = normalizar(eacpdim_alto_n1) 

#Escalas normalizadas escenario medio (Escenario 2) MÉTODO 2
E_acp_medio_dim_n4_NORM = normalizar(eacpdim_medio_n4) 
E_acp_medio_dim_n3_NORM = normalizar(eacpdim_medio_n3) 
E_acp_medio_dim_n2_NORM = normalizar(eacpdim_medio_n2) 
E_acp_medio_dim_n1_NORM = normalizar(eacpdim_medio_n1) 

#Escalas normalizadas escenario bajo (Escenario 3) MÉTODO 2
E_acp_bajo_dim_n4_NORM = normalizar(eacpdim_bajo_n4) 
E_acp_bajo_dim_n3_NORM = normalizar(eacpdim_bajo_n3) 
E_acp_bajo_dim_n2_NORM = normalizar(eacpdim_bajo_n2) 
E_acp_bajo_dim_n1_NORM = normalizar(eacpdim_bajo_n1) 


################################################################################

#Escalas normalizadas escenario alto (Escenario 1) MÉTODO 4
E_afm_alto_dim_n4_NORM = normalizar(eafmdim_alto_n4) 
E_afm_alto_dim_n3_NORM = normalizar(eafmdim_alto_n3) 
E_afm_alto_dim_n2_NORM = normalizar(eafmdim_alto_n2) 
E_afm_alto_dim_n1_NORM = normalizar(eafmdim_alto_n1) 

#Escalas normalizadas escenario medio (Escenario 2) MÉTODO 4
E_afm_medio_dim_n4_NORM = normalizar(eafmdim_medio_n4) 
E_afm_medio_dim_n3_NORM = normalizar(eafmdim_medio_n3) 
E_afm_medio_dim_n2_NORM = normalizar(eafmdim_medio_n2) 
E_afm_medio_dim_n1_NORM = normalizar(eafmdim_medio_n1) 

#Escalas normalizadas escenario bajo (Escenario 3) MÉTODO 4
E_afm_bajo_dim_n4_NORM = normalizar(eafmdim_bajo_n4) 
E_afm_bajo_dim_n3_NORM = normalizar(eafmdim_bajo_n3) 
E_afm_bajo_dim_n2_NORM = normalizar(eafmdim_bajo_n2) 
E_afm_bajo_dim_n1_NORM = normalizar(eafmdim_bajo_n1) 

#-------------------------------------------------------------------------------
#PSEUDO ERROR CUADRÁTICO MEDIO
#-------------------------------------------------------------------------------

PECM <- function(l1,l2,l3,l4){
  
  datos1<- unlist(l1) #Escalas n=1000
  datos2<- unlist(l2) #Escalas n=500
  datos3<- unlist(l3) #Escalas n=250
  datos4<- unlist(l4) #Escalas n=50
  
  datos <-c(datos1, datos2, datos3, datos4) #Unión de las escalas construidas con
  # el mismo método
  MediaItems <-mean(datos)  #Media de todos los puntajes de las escalas unidas
  
  datosEstimador<- data.frame(l1) #Escalas n=1000
  
  MediaEstimador= as.numeric(as.vector(apply(datosEstimador, 2, mean))) #Media de cada lista de escalas
  
  Sesgo_cuadrado = (mean(MediaEstimador) - MediaItems)^2   #Calculo del sesgo al cuadrado n=1000
  
  ECM<-0
  for (i in 1:1000) {
    ECM<-ECM+((MediaEstimador[i]-MediaItems)^2/1000) + Sesgo_cuadrado 
  }
  ECM #PECM
}

#PECM Método 1------------------------------------------------------------------
PECM_alto_ACP.items <- PECM(l1=E_acp_alto_global_n4_NORM,l2=E_acp_alto_global_n3_NORM,l3=E_acp_alto_global_n2_NORM,l4=E_acp_alto_global_n1_NORM)
PECM_medio_ACP.items <- PECM(l1=E_acp_medio_global_n4_NORM,l2=E_acp_medio_global_n3_NORM,l3=E_acp_medio_global_n2_NORM,l4=E_acp_medio_global_n1_NORM)
PECM_bajo_ACP.items <- PECM(l1=E_acp_bajo_global_n4_NORM,l2=E_acp_bajo_global_n3_NORM,l3=E_acp_bajo_global_n2_NORM,l4=E_acp_bajo_global_n1_NORM)

#PECM Método 3------------------------------------------------------------------
PECM_alto_AFM.items <- PECM(l1=E_afm_alto_global_n4_NORM,l2=E_afm_alto_global_n3_NORM,l3=E_afm_alto_global_n2_NORM,l4=E_afm_alto_global_n1_NORM)
PECM_medio_AFM.items <- PECM(l1=E_afm_medio_global_n4_NORM,l2=E_afm_medio_global_n3_NORM,l3=E_afm_medio_global_n2_NORM,l4=E_afm_medio_global_n1_NORM)
PECM_bajo_AFM.items <- PECM(l1=E_afm_bajo_global_n4_NORM,l2=E_afm_bajo_global_n3_NORM,l3=E_afm_bajo_global_n2_NORM,l4=E_afm_bajo_global_n1_NORM)

#PECM Método 2------------------------------------------------------------------
PECM_alto_ACP.dim <- PECM(l1=E_acp_alto_dim_n4_NORM,l2=E_acp_alto_dim_n3_NORM,l3=E_acp_alto_dim_n2_NORM,l4=E_acp_alto_dim_n1_NORM)
PECM_medio_ACP.dim <- PECM(l1=E_acp_medio_dim_n4_NORM,l2=E_acp_medio_dim_n3_NORM,l3=E_acp_medio_dim_n2_NORM,l4=E_acp_medio_dim_n1_NORM)
PECM_bajo_ACP.dim <- PECM(l1=E_acp_bajo_dim_n4_NORM,l2=E_acp_bajo_dim_n3_NORM,l3=E_acp_bajo_dim_n2_NORM,l4=E_acp_bajo_dim_n1_NORM)

#PECM Método 4------------------------------------------------------------------
PECM_alto_AFM.dim <- PECM(l1=E_afm_alto_dim_n4_NORM,l2=E_afm_alto_dim_n3_NORM,l3=E_afm_alto_dim_n2_NORM,l4=E_afm_alto_dim_n1_NORM)
PECM_medio_AFM.dim <- PECM(l1=E_afm_medio_dim_n4_NORM,l2=E_afm_medio_dim_n3_NORM,l3=E_afm_medio_dim_n2_NORM,l4=E_afm_medio_dim_n1_NORM)
PECM_bajo_AFM.dim <- PECM(l1=E_afm_bajo_dim_n4_NORM,l2=E_afm_bajo_dim_n3_NORM,l3=E_afm_bajo_dim_n2_NORM,l4=E_afm_bajo_dim_n1_NORM)


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
  
  # Iterar sobre las 1000 matrices 
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
  Medias <-cbind(mediasMinimo, mediasMaximo) # Media del mínimo y máximo de cada cluster
}

# Resumen de la media del mínimo y máximo de cada cluster de los puntajes de las escalas
#del MÉTODO 1 y 2
resumenAlto_ACP_items <- PuntosCorte(E_acp_alto_global_n4)
resumenMedio_ACP_items <- PuntosCorte(E_acp_medio_global_n4)
resumenBajo_ACP_items <- PuntosCorte(E_acp_bajo_global_n4)
resumenAlto_ACP_dim <- PuntosCorte(eacpdim_alto_n4)
resumenMedio_ACP_dim <- PuntosCorte(eacpdim_medio_n4)
resumenBajo_ACP_dim <- PuntosCorte(eacpdim_bajo_n4)

# Resumen de la media del mínimo y máximo de cada cluster de los puntajes de las escalas
#del MÉTODO 3 y 4
resumenAlto_AFM_items <- PuntosCorte(E_afm_alto_global_n4)
resumenMedio_AFM_items <- PuntosCorte(E_afm_medio_global_n4)
resumenBajo_AFM_items <- PuntosCorte(E_afm_bajo_global_n4)
resumenAlto_AFM_dim <- PuntosCorte(eafmdim_alto_n4)
resumenMedio_AFM_dim <- PuntosCorte(eafmdim_medio_n4)
resumenBajo_AFM_dim <- PuntosCorte(eafmdim_bajo_n4)

