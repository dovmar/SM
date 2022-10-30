# Dovydas Martinkus
# Duomenų Mokslas 4k. 1gr.
# 3 uzduotis

###


# Perkelties metodas


## Funkciju aprasymas


vyraujanti <- function(A) {
  result <- TRUE
  A <- abs(A)
  
  if (A[1,1] <= A[1,2]) {
    result <- FALSE
  }
  
  for ( i in 2:(nrow(A)-1) ) {
    if(A[i,i] < A[i,i-1] + A[i,i+1]) {
      result <- FALSE
    }
  }
  n <- nrow(A)
  if (A[n,n] < A[n,n-1]) {
    result <- FALSE
  }
  
  return(result)
}




perkelties <- function(A,B) {
  
  if (!vyraujanti(A)) {
    print("TLS isstrizaine nera vyraujanti")
  }
  
  p <- -1 * A[1,2] / A[1,1]
  q <- B[1] / A[1,1]
  print(A)
  print(B)
  
  for ( i in 2:(nrow(A)-1) ) {
    p_i <- -1* A[i,1+i] / (A[i,i] + A[i,i-1]*p[i-1])
    q_i <- (B[i]-A[i,i-1]*q[i-1]) / (A[i,i] + A[i,i-1]*p[i-1])
    
    p <- c(p,p_i)
    q <- c(q,q_i)
  }
  n <- nrow(A)
  q_n <- (B[n] - A[n,n-1]*q[n-1]) / (A[n,n] + A[n,n-1]*p[n-1])
  
  
  x <- numeric(n)
  x[n] <- q_n
  
  for (i in seq(n-1,1)) {
    x[i] <- p[i]*x[i+1] + q[i]  
  }
  
  return(x) 
}





A <- matrix(c(3, 1, 0, 0,
              -1, 4, 3, 0,
              0, 2, 4, -1,
              0, 0, 2, -3),
            ncol=4,nrow=4,byrow=TRUE)

B <- matrix(c(2,-2,1,-1),ncol=1)




x <- perkelties(A,B)


## gauto sprendinio patikrinimas
palyginimas <- cbind(A %*% matrix(x,ncol=1),B)
colnames(palyginimas) <- c("Gautas B","Norimas B")
t(palyginimas)












# Kubinis splainas




funkcija <- function(x) {
  exp(-x)*(x^3+2)
}


interpoliavimo_taskai <- function(func,n,a,b) {
  step <-(b-a)/n
  x <- a + step*(0:10)
  y <- funkcija(x)
  return(data.frame(x=x,y=y))
}




kubinis_splainas <- function(x,y) {
  n <- length(x)-1
  h <- diff(x)
  y_diff <- diff(y)
  B <- numeric(n-1) 
  A <- matrix(nrow=n-1,ncol=n-1)
  
  print(h)
  print(y_diff)
  for ( i in 1:(n-1) ) {
    
    if (i == 1) {
      row <- c(2*(h[i]+h[i+1]),h[i+1],rep(0,n-1-i-1))
    }
    
    else if (i == n-1) {
      row <- c(rep(0,n-1-2),h[i],2*(h[i]+h[i+1]))
    }
    
    else {
      row <- c(c(rep(0,i-2),h[i],2*(h[i]+h[i+1]),h[i+1],rep(0,n-1-i-1)))
    }
    
    b_row <- 6*((y[i+2]-y[i+1])/h[i+1] - (y[i+1]-y[i])/h[i])
    
    print(i)
    A[i,] <- row
    B[i] <- b_row
  } 

  
  g <- c(0,perkelties(A,B),0)
  
  G <- g[1:length(g)-1] / 2
  
  e <- y_diff / h - 1/6*g[2:length(g)]*h - 1/3*g[1:length(g)-1]*h
  
  H <- diff(g) / (6 * h)
  

  func <- function(z) {
    results <- c()
    for (zz in z) {
      for ( i in 0:(n-1) ) {
        if (x[i+1] <= zz & zz <= x[i+2]) {
          results <- c(results,y[i+1] + e[i+1]*(zz - x[i+1]) + G[i+1]*(zz - x[i+1])^2 + H[i+1]*(zz - x[i+1])^3)
          break
        }
        if (i == n-1) {
          print("Funkcijos argumentas ne is tinkamo intervalo")
          results <- c(results,NA)
        }
      }
    }
    return(results)
  }
  
}




a <- -1.
b <- 3
n <- 10

lentele <- interpoliavimo_taskai(funkcija,n,a,b)




gautas_splainas <- kubinis_splainas(lentele$x,lentele$y)

r_splainas <- splinefun(lentele$x,lentele$y)
library(tidyverse)
library(latex2exp)


rezultatai <- tibble(x = seq(-1, 3, 0.1),
           `Funkcija` = funkcija(x),
           `Apskaičiuotas splainas` = gautas_splainas(x),
           `splinefun splainas` = r_splainas(x))


rezultatai2 <- tibble(x = lentele$x,
                      `Funkcija` = funkcija(x),
                      `Apskaičiuotas splainas` = gautas_splainas(x),
                      `splinefun splainas` = r_splainas(x))

rezultatai <- rezultatai %>% pivot_longer(2:4,names_to = " ",values_to="y")
rezultatai$` `<- factor(rezultatai$` `,levels=c("Funkcija","Apskaičiuotas splainas","splinefun splainas"))


rezultatai2 <- rezultatai2 %>% pivot_longer(2:4,names_to = " ",values_to="y")
rezultatai2$` `<- factor(rezultatai2$` `,levels=c("Funkcija","Apskaičiuotas splainas","splinefun splainas"))

ggplot(rezultatai, aes(x,y,color=` `)) +
  geom_line() + geom_point(data = rezultatai2, aes(x,y,color=` `)) +
  labs(title=TeX("e^{-x}(x^3+2)"),
       subtitle = "Funkcijos ir jos interpoliavimo naudojant kubinius splainaus grafikai") +
  theme_minimal(base_size = 16) 
  



ggplot(rezultatai, aes(x,y,color=` `)) +
  geom_line() +  geom_point(data = rezultatai2, aes(x,y,color=` `)) +
  labs(title=TeX("e^{-x}(x^3+2)"),
       subtitle = "Funkcijos ir jos interpoliavimo naudojant kubinius splainaus grafikai") +
  theme_minimal(base_size = 16) + facet_wrap(vars(` `)) 

efun(lentele$x,lentele$y)







