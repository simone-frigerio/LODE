#Esempio di RW con 2 AO e 2 LS
n=100
jmp<-numeric(100)
ao<-numeric(100)
sig_eps <- runif(1,1,100)
sig_eta <- runif(1,1,4)*sig_eps/2
jmp_size <- runif(1,2.5,7.5)
ao_size <-runif(1,2.5,7.5) 
x_jmp<-sample(5:100,1)
xx_jmp<-sample((5:100)[-c(x_jmp-1,x_jmp,x_jmp+1)],1)
jmp[c(x_jmp,xx_jmp)]<-sample(c(-1,1),2,replace=T)*jmp_size*sig_eta
x_ao<-sample((5:100)[-c(x_jmp,xx_jmp)],2)
ao[x_ao] <- sample(c(-1,1),2,replace = T) * ao_size * sig_eps
mu<-cumsum(rnorm(n,sd=sig_eta)+jmp)+runif(1,-100*sig_eps,100*sig_eps)
y <- mu + rnorm(n, sd = sig_eps) + ao
plot(y,type="l",main="Esempio di RW con 2 AO e 2 LS")
lines(mu,col="red")