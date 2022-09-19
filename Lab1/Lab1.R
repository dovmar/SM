# Dovydas Martinkus
# Duomenų Mokslas 4k. 1gr.


func <- function(x) {
 x^5 + x - 3 
}


derivative <- function(x) {
  5*x^4 + 1
}



####


intervalas <- function(an, bn, cn, func) {
  if (func(an) * func(cn) < 0) {
    return(c(an, cn))
  } else {
    return(c(cn, bn))
  }
}


pusiaukirtos <- function(a0, b0, func, eps) {
  n <- 0
  c0 <- mean(c(a0, b0))
  a <- a0
  b <- b0
  c <- c0

  repeat {
    
    if ((abs(a[n+1] - b[n+1]) / 2 > eps)) {
      
      if (func(c[n+1]) == 0) {
        return(data.frame(a,b,c))
      }

      naujas_intervalas <- intervalas(a[n+1], b[n+1], c[n+1], func)
      a <- c(a,naujas_intervalas[1])
      b <- c(b,naujas_intervalas[2])
      c_n <- mean(c(a[n+2],b[n+2]))
      c <- c(c, c_n)

      n <- n + 1
      
    } else {
      return(data.frame(a,b,c))
    }
  }
}



pusiaukirtos_lentele <- function(x) {
  cbind(n=seq(0,lengths(x)[1]-1),
        x,
        y=func(x$c),
        abs(x$a-x$b)/2)
}






###


niutono <- function(x0, func, deriv, eps) {
  x <- x0
  n <- 0

  repeat {
    x_n <- x[n+1] - func(x[n+1]) / deriv(x[n+1])
    x <- c(x, x_n)

    if ((abs(x[n+2] - x[n+1]) > eps)) {
      n <- n + 1
    } else {
      return(x)
    }
  }
}



niutono_lentele <- function(x) {
  data.frame(n=seq(0,length(x)-1),
             x_n =x,
             y = func(x),
             `abs(x_n-x_n+1)`=c(0,diff(x)),
             check.names = FALSE)
}





library(ggplot2)
library(ggrepel)
library(latex2exp)



#####

eps <-  0.0001


# pradinis grafikas
ggplot(data.frame(x = seq(-3, 3, 0.1)), aes(x)) +
  geom_function(fun = func, colour = "black") + 
  geom_hline(yintercept = 0, color = "red") + 
  theme_minimal(base_size = 16) +
  labs(title=TeX("x^5+x-3"),subtitle = "Pradinis funkcijos grafikas")


ggplot(data.frame(x = seq(-3, 3, 0.1)), aes(x)) +
  geom_function(fun = derivative, colour = "black") + 
  geom_hline(yintercept = 0, color = "red") + 
  theme_minimal(base_size = 16) +
  labs(title=TeX("5x^4+1"),subtitle = "Išvestinės grafikas")


####


res <- pusiaukirtos(-2,0,func)
func(res$c[length(res$c)])


ggplot(data.frame(x = seq(-2, 2, 0.1)), aes(x)) +
  geom_function(fun = func, colour = "black") +
  geom_hline(yintercept = 0, color = "red") + 
  theme_minimal(base_size = 16) +
  geom_point(data=pusiaukirtos_lentele(res)[1:8,],
             aes(x=c,y=y)) + 
    geom_text_repel(data=pusiaukirtos_lentele(res)[1:8,],aes(x=c,y=y,label=paste0("x",n))) +
  labs(title=TeX("x^5+x-3"),subtitle = "Pusiaukirtos metodas")



pusiaukirtos_lentele(res)

####


xn <- niutono(0,func,derivative,eps)
func(xn[length(xn)])


ggplot(data.frame(x = seq(-3, 3, 0.1)), aes(x)) +
  geom_function(fun = func, colour = "black") +
  geom_hline(yintercept = 0, color = "red") + 
  theme_minimal(base_size = 16) +
  geom_point(data=niutono_lentele(xn)[1:8,],
             aes(x=x_n,y=y)) + 
  geom_text_repel(data=niutono_lentele(xn)[1:8,],aes(x=x_n,y=y,label=paste0("x",n))) +
  labs(title=TeX("x^5+x-3"),subtitle = "Niutono metodas")



niutono_lentele(xn)



####



palyginimas <- c(tail(res$c,1),tail(xn,1),uniroot(func,c(-3,3))$root)
names(palyginimas) <- c("Pusiaukirtos","Niutono","Uniroot()")
palyginimas



library(readr)


write_csv(round(niutono_lentele(xn),7),"n.csv")



write_csv(round(pusiaukirtos_lentele(res),7),"p.csv")


palyginimas

