# Dovydas Martinkus
# Duomenų Mokslas 4k. 1gr.
# 2 uzduotis

###


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


zeidelio_lentele <- function(xn,A,B) {
  
  norms <- apply(xn,1,function(x) norm(matrix(A %*% x - B,byrow=TRUE)))
  xn <- cbind(xn,norms)
  df <- data.frame(xn)
  names_df <- names(df)
  names_df[length(names_df)] <- "max norm"
  names(df) <- names_df
  return(df)
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



eps <- 0.01

x0 <- rep(0,length(B))




xn <- zeidelio(x0,D,L,U,B,eps)

zeidelio_lentele(xn,A,B)




# gauto artinio patikrinimas
palyginimas <- cbind(A %*% matrix(xn[nrow(xn),]),B)
colnames(palyginimas) <- c("Gautas B","Norimas B")
t(palyginimas)





library(Rlinsolve)

# artiniu palyginimas
palyginimas <- cbind(matrix(xn[nrow(xn),]),lsolve.jacobi(A,B)$x)
colnames(palyginimas) <- c("Zeidelio","R funkcija gautas")
t(palyginimas)

