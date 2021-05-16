codgen <- function(DNA = c("A", "T", "G", "C")){
  
  # Mayo de 2021
  # Creado por Melissa Bedoya, Yesica Montoya y David Quiroz
  
  # DESCRIPCION:
  # Esta función tiene como objetivo realizar automáticamente todos estos procesos,
  # partiendo de una secuencia de bases nitrogenadas de DNA, transformándola en RNA
  # y en la cadena de aminoácidos correspondiente.
  
  # El codigo usa condicionales anidados del tipo `if ... else` para los criterios
  # de seleccion del proceso de transcripcion de DNA a RNA.
  
  # El codigo usa el comando `sub` para reemplazar las tripletas de bases nitrogenadas
  # del RNA por las siglas del aminoacido correspondiente.
  
  # ARGUMENTO:
  # DNA: debe ser un vector de longitud igual a 3 o múltiplo de 3, que contenga las
  # letras 'A', 'T', 'C' o 'G', en cualquier orden y sin importar la cantidad de
  # repeticiones.
  
  
  # CODIGO:
  n <- length(DNA)
  RNA <- rep(NA, n)
  
  # Validacion del argumento DNA:
  prueba <- c("A", "T", "G", "C")
  h <- is.element(el = DNA, set = prueba)
  f <- sum( h ) 
  if(f < n) {
    i <- which(!h)
    mensaje <- paste0(
      "Las posiciones ", paste0(i, collapse = ", "),
      " del vector de entrada\ncontienen elementos diferentes de: ",
      paste0(prueba, collapse = ", ")
    )
    stop(mensaje) 
  }
  
  
  mensaje2 <- "El argumento 'DNA' debe ser un vector\nde longitud igual a 3 o multiplos de 3"
  check <- function(x) x == round(x)
  a <- length(DNA)/3
  if(check(a) == F) stop(mensaje2)
  
  
  # Ciclo para transcribir el DNA al RNA:
  for(i in 1:n){
    RNA[i] <- if(DNA[i] == "A"){
      "U"
    } else{
      if(DNA[i] == "T"){
        "A"
      } else{
        if(DNA[i] == "G"){
          "C"
        } else{
          if(DNA[i] == "C"){
            "G"
          }
        }
      }
    }
  }
  
  # Organizacion del RNA en tripletas:
  RNAm <- matrix(data = RNA, ncol = 3, byrow = TRUE)
  RNA_final <- matrix(data = NA, ncol = 1, nrow = nrow(RNAm))
  k <- nrow(RNAm)
  for(i in 1:k) RNA_final[i,] <- paste0(RNAm[i,], collapse = "")
  
  # Traduccion de RNA a aminoacidos:
  aa <- c(RNA_final[ ,1])
  
  aa <- sub(pattern = "UUU", replacement = "Phe", x = aa)
  aa <- sub(pattern = "UUC", replacement = "Phe", x = aa)
  aa <- sub(pattern = "UUA", replacement = "Leu", x = aa)
  aa <- sub(pattern = "UUG", replacement = "Leu", x = aa)
  aa <- sub(pattern = "CUU", replacement = "Leu", x = aa)
  aa <- sub(pattern = "CUC", replacement = "Leu", x = aa)
  aa <- sub(pattern = "CUA", replacement = "Leu", x = aa)
  aa <- sub(pattern = "CUG", replacement = "Leu", x = aa)
  aa <- sub(pattern = "AUU", replacement = "Ile", x = aa)
  aa <- sub(pattern = "AUC", replacement = "Ile", x = aa)
  aa <- sub(pattern = "AUA", replacement = "Ile", x = aa)
  aa <- sub(pattern = "AUG", replacement = "Met", x = aa)
  aa <- sub(pattern = "GUU", replacement = "Val", x = aa)
  aa <- sub(pattern = "GUC", replacement = "Val", x = aa)
  aa <- sub(pattern = "GUA", replacement = "Val", x = aa)
  aa <- sub(pattern = "GUG", replacement = "Val", x = aa)
  aa <- sub(pattern = "UCU", replacement = "Ser", x = aa)
  aa <- sub(pattern = "UCC", replacement = "Ser", x = aa)
  aa <- sub(pattern = "UCA", replacement = "Ser", x = aa)
  aa <- sub(pattern = "UCG", replacement = "Ser", x = aa)
  aa <- sub(pattern = "CCU", replacement = "Pro", x = aa)
  aa <- sub(pattern = "CCC", replacement = "Pro", x = aa)
  aa <- sub(pattern = "CCA", replacement = "Pro", x = aa)
  aa <- sub(pattern = "CCG", replacement = "Pro", x = aa)
  aa <- sub(pattern = "ACU", replacement = "Thr", x = aa)
  aa <- sub(pattern = "ACC", replacement = "Thr", x = aa)
  aa <- sub(pattern = "ACA", replacement = "Thr", x = aa)
  aa <- sub(pattern = "ACG", replacement = "Thr", x = aa)
  aa <- sub(pattern = "GCU", replacement = "Ala", x = aa)
  aa <- sub(pattern = "GCC", replacement = "Ala", x = aa)
  aa <- sub(pattern = "GCA", replacement = "Ala", x = aa)
  aa <- sub(pattern = "GCG", replacement = "Ala", x = aa)
  aa <- sub(pattern = "UAU", replacement = "Tyr", x = aa)
  aa <- sub(pattern = "UAC", replacement = "Tyr", x = aa)
  aa <- sub(pattern = "UAA", replacement = "STOP", x = aa)
  aa <- sub(pattern = "UAG", replacement = "STOP", x = aa)
  aa <- sub(pattern = "CAU", replacement = "His", x = aa)
  aa <- sub(pattern = "CAC", replacement = "His", x = aa)
  aa <- sub(pattern = "CAA", replacement = "Gln", x = aa)
  aa <- sub(pattern = "CAG", replacement = "Gln", x = aa)
  aa <- sub(pattern = "AAU", replacement = "Asn", x = aa)
  aa <- sub(pattern = "AAC", replacement = "Asn", x = aa)
  aa <- sub(pattern = "AAA", replacement = "Lys", x = aa)
  aa <- sub(pattern = "AAG", replacement = "Lys", x = aa)
  aa <- sub(pattern = "GAU", replacement = "Asp", x = aa)
  aa <- sub(pattern = "GAC", replacement = "Asp", x = aa)
  aa <- sub(pattern = "GAA", replacement = "Glu", x = aa)
  aa <- sub(pattern = "GAG", replacement = "Glu", x = aa)
  aa <- sub(pattern = "UGU", replacement = "Cys", x = aa)
  aa <- sub(pattern = "UGC", replacement = "Cys", x = aa)
  aa <- sub(pattern = "UGA", replacement = "STOP", x = aa)
  aa <- sub(pattern = "UGG", replacement = "Trp", x = aa)
  aa <- sub(pattern = "CGU", replacement = "Arg", x = aa)
  aa <- sub(pattern = "CGC", replacement = "Arg", x = aa)
  aa <- sub(pattern = "CGA", replacement = "Arg", x = aa)
  aa <- sub(pattern = "CGG", replacement = "Arg", x = aa)
  aa <- sub(pattern = "AGU", replacement = "Ser", x = aa)
  aa <- sub(pattern = "AGC", replacement = "Ser", x = aa)
  aa <- sub(pattern = "AGA", replacement = "Arg", x = aa)
  aa <- sub(pattern = "AGG", replacement = "Arg", x = aa)
  aa <- sub(pattern = "GGU", replacement = "Gly", x = aa)
  aa <- sub(pattern = "GGC", replacement = "Gly", x = aa)
  aa <- sub(pattern = "GGA", replacement = "Gly", x = aa)
  aa <- sub(pattern = "GGG", replacement = "Gly", x = aa)
  
  # Organizacion de los datos en tripletas para mostrar equivalencias:
  DNA1 <- matrix(data = DNA, ncol = 3, byrow = TRUE)
  DNA2 <- matrix(data = NA, ncol = 1, nrow = nrow(DNA1))
  for(i in 1:nrow(DNA1)) DNA2[i,] <- paste0(DNA1[i,], collapse = "")
  
  # Impresion de resultados:
  data.frame(DNA = DNA2, RNA = RNA_final, aa)
}
