n<-100
sig_eps <- runif(1,1,100)
sig_eta <- runif(1,1,4)*sig_eps/2
c_eps <- runif(1,2,6)
c_eta <- runif(1,2,6)
disturbi_eps<-cbind(rnorm(n,0,sd=sig_eps),rnorm(n,0,sd=c_eps*sig_eps))
eps<-apply(disturbi_eps, 1, sample,prob=c(0.9,0.1),size=1)
distrurbi_eta<-cbind(rnorm(n,0,sd=sig_eta),rnorm(n,0,sd=c_eta*sig_eta))
eta<-apply(distrurbi_eta, 1, sample,prob=c(0.9,0.1),size=1)

sig_beta<-runif(1,0.2,0.8)
b<-cumsum(rnorm(n, sd = sig_beta))+runif(1,-5,5) 
b<-b*sig_eps/3

periodicita<-4
stg<-numeric()
sig_stg<-runif(1,0.2,0.8)
stg[1:(periodicita-1)]<-runif(periodicita-1,-5,5)
for(i in periodicita:n){
  stg[i]<--sum(stg[(i-periodicita+1):(i-1)])+rnorm(1,0,sig_stg)
}
stg<-stg*sig_eps

mu <- cumsum(eta+b)+runif(1,-100*sig_eps,100*sig_eps)
y <- mu + eps + stg

plot(y,type="l",main="Esempio di LLT più stagionalità")
lines(mu,col="red")