# Dovydas Martinkus
# Duomenų Mokslas 4k. 1gr.
# 2 uzduotis

###


# Funkciju aprasymas

zeidelio <- function(x0,D,L,U,B,eps) {
  x <- matrix(x0)
  n <- 0

  repeat {
    x_n <- solve(D-L) %*% (U %*% matrix(x[,n+1])+B)
    x <- cbind(x, x_n)

    
    if (norm(matrix(matrix(x[,n+2]) - matrix(x[,n+1]),
                    byrow=TRUE),
             type = "M") > eps) {
       n <- n + 1
    }
    else {
      return(t(x))
    }
  }
}




jungtiniu_gradientu <- function(x0,A,B,eps) {
  x <- matrix(x0)
  z <- A %*% x - B
  p <- z
  r <- A %*% p
  tau <- c((t(z) %*% z) / (t(r) %*% p))
  n <- 0
  
  repeat {
    
    
    if( n>0 ){
      r_n <- A %*% matrix(p[,n+1])
      r <- cbind(r, r_n)
      tau_n <- (t(matrix(z[,n+1])) %*% matrix(z[,n+1])) / (t(matrix(r[,n+1])) %*% matrix(p[,n+1]))
      tau <- cbind(tau, tau_n)
     }

    
    x_n <- matrix(x[,n+1]) - tau[n+1]*matrix(p[,n+1])
    z_n <- matrix(z[,n+1]) - tau[n+1]*matrix(r[,n+1])
    
  
    x <- cbind(x, x_n)
    z <- cbind(z, z_n)
    
    
    if ( t(z_n) %*% z_n >= eps^2 ) {
      beta_n <- c((t(z_n) %*% z_n) / ( t(matrix(z[,n+1])) %*%  matrix(z[,n+1])))
      p_n <- z_n + beta_n*matrix(p[,n+1])
      p <- cbind(p, p_n)
      n <- n + 1
    }
    else {
      return(list(a=t(x),b=t(p)))
    }
  }
}





## lentele, parodanti artiniu paklaidas pagal iteracijas
lentele <- function(xn,A,B) {
  
  norms <- apply(xn,1,function(x) norm(matrix(A %*% x - B,byrow=TRUE),type="M"))
  xn <- cbind(0:(nrow(xn)-1),xn,norms)
  return(xn)
}





D <- matrix(c(0.18, 0.57, 0.38, 0.42,
              0.57, 0.95, 0.70, 0.44,
              0.38, 0.70, 0.37, 0.18,
              0.42, 0.44, 0.18, 0.40),
            ncol=4,nrow=4)


E <- diag(1, 4, 4)


B <- matrix(c(1.2,2,3,1.5))


A <- D + 0.1*(16+3)*E


U <- A
U[lower.tri(U,TRUE)] <- 0
U <- U * -1


L <- A
L[upper.tri(L,TRUE)] <- 0
L <- L * -1


D <- A
D[!(lower.tri(D,TRUE) & upper.tri(D,TRUE))] <- 0
D <- D



# konvergavimo sąlygų patikrinimas

max(abs(eigen(solve(D-L)%*%U)$values))

norm(solve(D-L)%*%U)





eps <- 0.0001

x0 <- rep(0,length(B)) # bet koks pradinis artinys






# Zeidelio metodas

xn_1 <- zeidelio(x0,D,L,U,B,eps)

round(lentele(xn_1,A,B),5)



palyginimas <- cbind(A %*% matrix(xn_1[nrow(xn_1),]),B)
colnames(palyginimas) <- c("Gautas B","Norimas B")
round(t(palyginimas),6)







# Jungtiniu gradientu metodas


result <- jungtiniu_gradientu(x0,A,B,eps)

xn_2 <- result[[1]]
p_n <- result[[2]]


round(lentele(xn_2,A,B),5)



palyginimas <- cbind(A %*% matrix(xn_2[nrow(xn_2),]),B)
colnames(palyginimas) <- c("Gautas B","Norimas B")
t(palyginimas)




## patikriname ar gauti p_n tikrai sudaro matricos A atzvilgiu jungtiniu vektoriu sistema

for(i in 1:nrow(p_n)) {
 for(j in 1:nrow(p_n)) {
   if(i != j) {
     print(all.equal(matrix(0)
                            ,t(A%*%p_n[i,]) %*% p_n[j,]))
   }
 } 
}







library(Rlinsolve)

# Artiniu palyginimas

palyginimas <- cbind(matrix(xn_1[nrow(xn_1),]),matrix(xn_2[nrow(xn_2),]),lsolve.jacobi(A,B)$x, lsolve.cg(A,B)$x)
colnames(palyginimas) <- c("Zeidelio","Jungtinių gradientų","R iteracinis metodas lsolve.jacobi","R variacinis metodas lsove.cg")
t(palyginimas)

