#Preparing the data
secuencias <- c("TTGCAAAT","ATGCATAT","ATTCAAAT","ATGCAAAA","ATGCAAAT","ATGCAATT","ATGCAATT","AAGCAAAG")
nt <- c("A","C","G","T")

bases <- c("A","C","G","T")
longseq <- nchar(secuencias[1])

#Creating the Position Weight Matrix and initializing it
PWM <- matrix(data = 0,
  nrow = length(bases),
  ncol = longseq
)
rownames(PWM) <- bases
colnames(PWM) <- 1:longseq

#Filling the Position Weight Matrix
for (secuencia in secuencias){
  sec <- strsplit(secuencia, split="")[[1]] 
    for (posicion in 1:longseq){
    base_actual <- sec[posicion]
    PWM[base_actual,posicion] <- PWM[base_actual,posicion] + 1
  }
}  


#Creating the Position Weight Matrix with values over 1 -> PWM2
PWM2 <- matrix( data = PWM,
                nrow = length(bases),
                ncol = longseq)

rownames(PWM2) <- bases
colnames(PWM2) <- 1:longseq

#Counting the number of sequences compared
nseq <- length(secuencias)


#Filling the PWM2 Matrix
for (fila in bases){
  for (columna in 1:longseq){
    valor1 <- PWM2[fila,columna]
    PWM2[fila,columna] <- valor1/nseq
  }
}


#Creating the Position Weight Matrix with pseudocounts (p=1) -> PPM (Position Probability Matrix)
PPM <- matrix(data=PWM2,
               nrow = length(bases),
               ncol = longseq
               )

rownames(PPM) <- bases
colnames(PPM) <- 1:longseq

#Function pseudocount (se aplica a PWM)
pseudocount <- function(valor2,numseq, p = 1, n = 4){
  numerador <-  valor2 + (p/n)
  denominador <- numseq + p
  probabilidad <- numerador/denominador
  return(probabilidad)
}


for (fila in bases){
  for (columna in 1:longseq){
    valor3 <- PWM[fila,columna]
    PPM[fila,columna] <- pseudocount(valor2 = valor3, numseq = nseq)
  }
}

#Creating the Null Matrix
nullM <- matrix(data = 0,
                nrow = length(bases),
                ncol = longseq
                )

rownames(nullM) <- bases
colnames(nullM) <- 1:longseq


#Filling the Null Matrix
for (fila in bases){
  for (columna in 1:longseq){
    #Filling cells as: P(obtaining a base by chance) = 1/number bases (4)
    nullM[fila,columna] <- 1/length(bases)
  }
}


#Creating the Log Likelihood Ratio Matrix
LogM = matrix(data = 0,
                nrow = length(bases),
                ncol = longseq
                )
row.names(LogM) <- bases
colnames(LogM) <- 1:longseq

#Filling the Log Likelihood Ratio Matrix
for (fila in bases){
  for (columna in 1:longseq){
    #Filling cells as: ln(P(PPM)/P(nullM))
    LogM[fila,columna] <- log(PPM[fila,columna]/nullM[fila,columna])
  }
}


#Creating the Information Content Matrix (ICM)
ICM <- matrix(data = 0,
              nrow = length(bases),
              ncol = nseq)
rownames(ICM) <- bases
colnames(ICM) <- 1:longseq

#Calculating ICfinal
ICtotal <- 2
U <- function(P){
  producto <- ifelse ( 
    #log2(0) can't be calculated -> if P = 0, product = 0 
    P == 0,
    0,
    #if P != 0, calculate the product
    P*log2(P)
  )
  resultado <- sum(producto)
  return(-resultado)
}
#Vector with ICfinal-values for each column
vector_ICfinal <- c()

#Transforming each column in vectors
for (i in 1:longseq){
  columna <- PWM2[,i]
  U_i <- U(columna)
  ICfinal <- ICtotal- U_i 
  vector_ICfinal <- c(vector_ICfinal, ICfinal)
}


#Filling the Information Content Matrix (ICM)
for (fila in bases){
  for (columna in 1:longseq){
    #Calculating ICfinal
    ICM[fila,columna] <- PWM2[fila,columna]*vector_ICfinal[columna]
  }
}

print(PWM)
print(PPM)
print(LogM)
print(ICM)

#Checking the results of the Information Content Matrix (ICM) with the corresponding Sequence Logo
library(ggseqlogo)
ggseqlogo(secuencias)