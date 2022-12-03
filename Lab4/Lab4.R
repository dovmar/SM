# Dovydas Martinkus
# Duomenų Mokslas 4k. 1gr.
# 4 uzduotis

###



## Funkciju aprasymas

funkcija <- function(x,u) {
  x*sin(2*u) + x^2
}

tinklas <- function(start,end,step) {
  N <- (end-start) / step
  if (N%%1!=0) {
    print('Netinkamas zingsnio dydis')
  }
  step * 0:N
}



runges_kuto <- function(func,m,u0,zingsnis,a,b,sigma,start=0,
                        end=1) {
 
  tinklas <- tinklas(start,end,zingsnis)
  n <- length(tinklas)-1
  
  y <- u0
  for (i in 1:n) {
   y_i <- y[i] + zingsnis * k_m(func,m,y[i],tinklas[i],
                         zingsnis,a,b,sigma) 
  y <- c(y,y_i)
  }
  data.frame(x=tinklas,y=y,zingsnis=zingsnis,m=m)
}



k_m <- function(func,m,y_n,t_n,zingsnis,
                a_array,b_matrix,sigma_array) {
  k <- numeric()
  for (i in 1:m) {
    if (i == 1) {
      a = 1
      b_sum = 1
    }
    else {
      a = a_array[i]
      b_sum <- b_matrix[i,1:(i-1)] %*% k 
    }
    k_i <- func(t_n + zingsnis*a, y_n + zingsnis*b_sum)
    k <- c(k,k_i)
  }
  k %*% sigma_array[1:m]
}



runges_paklaida <- function(func,m,u0,zingsnis,a,b,sigma,p,start=0,
                            end=1) {
  result_a<-runges_kuto(func,m,u0,zingsnis,a,b,sigma,start=0,
              end=1)
  result_b<-runges_kuto(func,m,u0,zingsnis*2,a,b,sigma,start=0,
                 end=1)
  
  y_tau <- result_a[lengths(result_a)[1],2]
  y_2tau <- result_b[lengths(result_b)[1],2]
  data.frame(m=m,zingsnis=zingsnis,paklaida=abs(y_2tau-y_tau) / (2^p-1))
}





# 3-pakopis


a_3 <- c(0,1/2,1)
sigma_3 <- c(1/6,4/6,1/6)
b_3 <- matrix(c(0,0,
            1/2,0,
            -1,2),nrow=3,byrow = TRUE)


runges_kuto(funkcija,3,0,0.025,a_3,b_3,sigma_3)

runges_paklaida(funkcija,3,0,0.025,a_3,b_3,sigma_3,3)




# 4-pakopis 


a_4 <- c(0,1/2,1/2,1)
sigma_4 <- c(1/6,2/6,2/6,1/6)
b_4 <- matrix(c(0,0,0,
            1/2,0,0,
            0,1/2,0,
            0,0,1),nrow = 4,byrow = TRUE)



runges_kuto(funkcija,4,0,0.025,a_4,b_4,sigma_4)


runges_paklaida(funkcija,4,0,0.025,a_4,b_4,sigma_4,4)




# Rezultatu palyginimas

library(deSolve)


funkcija2 <- function(x,u,parms=NULL) {
  list(funkcija(x,u))
}

r_runges_kuto <- function(func,rk,u0,zingsnis,start=0,end=1) {
  tinklas <- tinklas(start,end,zingsnis)
  y<- ode(times = tinklas, y = u0, func = func,
       parms = NULL, method = rkMethod('rk4'))
  data.frame(x=tinklas,y=y[,2],zingsnis=zingsnis,m=rk)
}





rezultatai <- NULL
paklaida <- NULL

for (i in c(0.1,0.05,0.025,0.0125)) {
  rezultatai <-rbind(rezultatai,
                     runges_kuto(funkcija,3,0,i,a_3,b_3,sigma_3))
  rezultatai <-rbind(rezultatai,
                     runges_kuto(funkcija,4,0,i,a_4,b_4,sigma_4))
  rezultatai <-rbind(rezultatai,
                     r_runges_kuto(funkcija2,'rk4',0,i))
  
  paklaida <- rbind(paklaida,
                    runges_paklaida(funkcija,3,0,i,a_3,b_3,sigma_3,3),
                    runges_paklaida(funkcija,4,0,i,a_4,b_4,sigma_4,4))
  
}

rezultatai$zingsnis <- factor(rezultatai$zingsnis)
rezultatai$m <- factor(rezultatai$m)



library(ggplot2)

ggplot(subset(rezultatai,1==1),aes(x,y,color=m)) +
  geom_line()  + facet_wrap(vars(zingsnis),labeller = 'label_both') +
  theme_minimal(base_size=20) + labs(title='')













