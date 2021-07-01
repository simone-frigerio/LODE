library(KFAS)
library(dplyr)
library(DiceKriging)
library(doParallel)
library(lhs)

#Iperettangolo dei lambda
lower <- c(0.1,0.1)
upper <- c(2,2)
#Numero iterazioni della BO
n.trials<-50
#Numero di valori iniziali da generare con LHS
vallhs<-17

#Plot
timeUKd<-time(UKDriverDeaths)[79:192]
UKd<-(UKDriverDeaths[79:192])
plot(timeUKd,UKd,type="l",xlab="Anno",ylab="Numero di feriti gravi o morti",main="Incidenti stradali UK")

y<-log(UKd)
n<-length(y)
per<-12 #periodicita

#Funzione obiettivo da minimizzare
obj <- function(pars, mod, l1, l2){
  n <- (length(pars) - 4) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  mod$Q[2, 2, ]<- array(pars[(2*n+3)]^2,c(1, 1, n))
  mod$Q[3, 3, ]<- array(pars[(2*n+4)]^2,c(1, 1, n))
  -KFS(mod, filtering = "none", smoothing = "none")$logLik +
    l1 * sum(abs(pars[2:(n + 1)])) +
    l2 * sum(abs(pars[(n + 3):(2*n + 2)]))
}
#Funzione che calcola lo smoother con le varianze aggiuntive
make_mod <- function(pars, mod, kfs = TRUE) {
  n <- (length(pars) - 4) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  mod$Q[2, 2, ]<- array(pars[(2*n+3)]^2,c(1, 1, n))
  mod$Q[3, 3, ]<- array(pars[(2*n+4)]^2,c(1, 1, n))
  if (kfs) return(list(mod = mod, kfs = KFS(mod, smoothing = "state")))
  list(mod = mod, kfs = NULL)
}
#Gradiente della funzione obiettivo
grad <- function(pars, mod, l1, l2) {  
  n <- (length(pars) - 4) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  mod$Q[2, 2, ]<- array(pars[(2*n+3)]^2,c(1, 1, n))
  mod$Q[3, 3, ]<- array(pars[(2*n+4)]^2,c(1, 1, n))
  kfs <- KFS(mod, smoothing = "disturbance", simplify = FALSE)
  h <- mod$H[1, 1, ]
  e2D <- (as.numeric(kfs$epshat)^2 - (h - as.numeric(kfs$V_eps))) / (h^2)
  r2N <- (kfs$r)^2 - apply(kfs$N,3,diag)
  d_eps <- -sum(e2D) * pars[1]
  d_eta <- -sum(r2N[1,]) * pars[n + 2]
  v_eps <- -e2D * pars[2:(n + 1)] +
    l1 * sign(pars[2:(n + 1)])
  v_eta <- -r2N[1,2:(n + 1)] * pars[(n + 3):(2*n + 2)] +
    l2 * sign(pars[(n + 3):(2*n + 2)])
  beta <- -sum(r2N[2,]) * pars[(2*n+3)]
  stag <- -sum(r2N[3,]) * pars[(2*n+4)]
  c(d_eps, v_eps, d_eta, v_eta,beta,stag)
}


#Log(verosimiglianza)
loglik <- function(pars, mod) {
  mod$H[1, 1, ] <- pars[1]^2
  mod$Q[1, 1, ] <- pars[2]^2
  mod$Q[2, 2, ] <- pars[3]^2
  mod$Q[3, 3, ] <- pars[4]^2
  -KFS(mod, filtering = "none", smoothing = "none")$logLik
}
#Riscaliamo la serie
mod <- SSModel(y ~ SSMtrend(2, list(array(0, c(1, 1, (n))),0), a1 = c(y[1],0), P1 = c(var(y),0,0,var(y)), P1inf = c(0,0,0,0))+SSMseasonal(per,1,"dummy"),
               H = array(0, c(1, 1, n)))
mle <- optim(c(sd(diff(y,differeces=per))/10, sd(diff(y,differeces=per))/10, sd(diff(diff(y,differeces=per)))/10, sd(diff(y,differeces=per))/10), loglik, mod = mod, method = "L-BFGS-B")
mle$convergence
y1<-y[1]
risc<-abs(mle$par[1])
y<-(y-y[1])*5/abs(mle$par[1])
mod <- SSModel(y ~ SSMtrend(2, list(array(0, c(1, 1, (n))),0), a1 = c(y[1],0), P1 = c(var(y),0,0,var(y)), P1inf = c(0,0,0,0))+SSMseasonal(per,1,"dummy"),
               H = array(0, c(1, 1, n)))
mle <- optim(c(c(sd(diff(y,differeces=per))/10, sd(diff(y,differeces=per))/10, sd(diff(diff(y,differeces=per)))/10, sd(diff(y,differeces=per))/10)), loglik, mod = mod, method = "L-BFGS-B")
print(abs(mle$par[1]))

#BO
#Funzione che definisce valori iniziali e dominio
fun_iniz_dom<-function(y){
  mle <- optim(c(sd(diff(y))/2, sd(diff(y))/2, 0.5, 0.3)*2, loglik, mod = mod, method = "L-BFGS-B",control = list(maxit = 500))
  mod_mle <- mod
  mod_mle$H[1, 1, ] <- mle$par[1]^2
  mod_mle$Q[1, 1, ] <- mle$par[2]^2
  mod_mle$Q[2, 2, ] <- mle$par[3]^2
  mod_mle$Q[3, 3, ] <- mle$par[4]^2
  kfs_dist <- KFS(mod_mle,smoothing = "disturbance")
  ineta<-abs(kfs_dist$etahat[,1])
  ineps<-abs(kfs_dist$epshat)
  init<-c(sqrt(mod_mle$H[1,1,1])*0.666,ineps,sqrt(mod_mle$Q[1,1,1])*0.666,ineta,sqrt(mod_mle$Q[2,2,1]),sqrt(mod_mle$Q[3,3,1]))
  low<--8*init
  up<-8*init
  dom<-matrix(c(low,up),ncol=2)
  newList <- list("iniziali" = init, "up" = up,"low"=low)
  return(newList)
}

iniz_dom<-fun_iniz_dom(y)#rida i valori iniziale per mettere in input nella prossima funzione

#Funzione che calcola il BIC
bic_kalman<-function(iniz_dom,l1,l2){
  k1<-100#servono per definire il bic
  k2<-1
  init<-iniz_dom$iniziali#valori iniziali
  up<-iniz_dom$up#dominio superiore
  low<-iniz_dom$low#dominio inferiore
  
  foo1 <- try(optim(init, obj, grad,mod = mod, l1 = l1, l2 = l2,method = c("BFGS"),control = list(maxit = 50*2*n,factr=1e-10)),silent=TRUE)
  if (class(foo1) == "try-error"){
    foo1 <- try(optim(init, obj, grad,mod = mod, l1 = l1, l2 = l2,lower=low,upper=up,method = c("L-BFGS-B"),control = list(maxit = 50*2*n,factr=1e-10)))
    print("L-BFGS-B")   
  }
  foo1$loglik <- make_mod(foo1$par, mod, TRUE)$kfs$logLik
  x<-foo1
  bic<--2*x$loglik + log(n)*k2*sum(abs(x$par[1:n+1])>=abs(x$par[1]/k1))+log(n)*k2*sum(abs(x$par[(n+2):(2*n+2)])>=abs(x$par[(n+2)]/k1))
  foo1$bic<-bic
  return(foo1)
}
#Funzione di acquisizione
CB <- function( x, gp, beta_=1, UCB=T ) {
  if( !is.matrix(x) )
    x <- data.frame( t(x) )
  pred <- predict( gp, x, "UK", checkNames=F )
  if( UCB )
    res <- pred$mean + sqrt(beta_)*pred$sd
  else
    res <- - (pred$mean - sqrt(beta_)*pred$sd) # because all the acquisitions are maximized
  
  return(res)
}


if(tolower(.Platform$OS.type) != "windows"){
  cl <- makeCluster(spec=detectCores(), type="FORK", outfile="")  
} else {
  cl <- makeCluster(spec=detectCores(), outfile="")
}

setDefaultCluster(cl=cl)
registerDoParallel(cl)

kernel <- "gauss"
acquisition <- "CB" 


cat("\n> ***** Starting BO *****\n")
n0 <- round(vallhs) 

#Generazione dei punti iniziali con LHS
set.seed(1)
design <- lhs::maximinLHS(n=n0,k=length(lower))
for( i in 1:ncol(design) ){
  design[,i] <- design[,i] * ( upper[i] - lower[i] ) + lower[i]
}

bo.elapsed <- Sys.time()
tot_foo <- foreach( i = 1:n0, .packages="KFAS", .combine="c" ) %dopar% {
  bic_kalman( iniz_dom=iniz_dom, l1=design[i,1], l2=design[i,2] )
}
y.bic<-as.numeric(tot_foo[which(names(tot_foo)=="bic")])
bo.elapsed <- difftime(Sys.time(),bo.elapsed,units="secs")
print(bo.elapsed)

acq.fun <- CB

cat("> Starting Bayesian Optimization (BO)...\n")
cat("> ",n.trials-n0,"remaining trials...\n")
bo.iter.elapsed <- NULL
while( nrow(design)<n.trials ) {
  
  tmp <- Sys.time()
  
  gp.fit <- km( design=data.frame(design), response=y.bic,
                covtype=kernel, nugget.estim=T,
                control=list(trace=0))
  

  fn.scale <- -1 
  x0 <- NULL
  for( i in 1:ncol(design) )
    x0 <- c( x0, runif(1,lower[i],upper[i]) )
  #Ottimizzazione della funzione di acquisizione
  opt.res <- optim( par=x0, fn=acq.fun, gr=NULL,
                    gp=gp.fit, UCB=F, #BIC è minimizzato
                    lower=lower, upper=upper, control=list(trace=0,fnscale=fn.scale) )
  
  design <- rbind( design, opt.res$par )
  #Calcolo il BIC nel nuovo punto
  output<-bic_kalman( iniz_dom=iniz_dom, l1=opt.res$par[1], l2=opt.res$par[2] )
  tot_foo<-c(tot_foo,output)
  y.bic <- c( y.bic,output$bic) 
  
  tmp <- difftime( Sys.time(), tmp, units="secs" )
  bo.iter.elapsed <- c( bo.iter.elapsed, tmp )
  print(nrow(design))# print del numero dell'iterazione
}
cat("> BO time:",bo.elapsed+sum(bo.iter.elapsed),"seconds\n")

setDefaultCluster(cl=NULL)
stopCluster(cl)


tot_par<-tot_foo[which(names(tot_foo)=="par")]
min<-which.min(y.bic)
opt_par<-as.numeric(unlist(tot_par[min]))#parametri con cui si ottiene il BIC minore

#eseguo il kalman smoother e riscalo i risultati
stima_lode <- make_mod(opt_par, mod)$kfs$alphahat[, 1]
stima_lode <- exp(stima_lode*risc/5+y1)

#guardo gli outlier individuati
plot(opt_par,och=20,main="parametri ottimi")
w_ao<-which(abs(opt_par[2:(n+1)])>abs(opt_par[1]/100))
w_level<-which(abs(opt_par[(n+3):(2*n+2)])>abs(opt_par[n+2]/100))
cat("Istante in cui è avvenuto il cambio di livello: ",timeUKd[w_level+1]) 
cat("Istante in cui è presente un outlier additivo: ",timeUKd[w_ao]) 

#plot dei risultati
plot(timeUKd,UKd,type="l",xlab="Anno",ylab="Numero di feriti gravi o morti",main="Incidenti stradali UK")
ylode<-ts(stima_lode,start = timeUKd[1],frequency = 12)
lines(ylode,col="blue",lwd=2)




