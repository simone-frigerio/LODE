library(KFAS)
library(dplyr)
library(DiceKriging)
library(doParallel)
lower <- c(0.1,0.1)
upper <- c(2,2)
n.trials<-50
vallhs<-17
y1<-numeric()
risc<-numeric()

#plot
plot(Nile,main="Nilo",lwd=2,xlab="Anno",ylab="Portata fiume",cex.lab=1.3,cex.main=1.5)

y<-as.numeric(Nile)
n<-length(y)

#funzioni varie
obj <- function(pars, mod, l1, l2) {
  n <- (length(pars) - 2) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  -KFS(mod, filtering = "none", smoothing = "none")$logLik +
    l1 * sum(abs(pars[2:(n + 1)])) +
    l2 * sum(abs(pars[(n + 3):(2*n + 2)]))
}
make_mod <- function(pars, mod, kfs = TRUE) {
  n <- (length(pars) - 2) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  if (kfs) return(list(mod = mod, kfs = KFS(mod, smoothing = "state")))
  list(mod = mod, kfs = NULL)
}

grad <- function(pars, mod, l1, l2) {  
  n <- (length(pars) - 2) / 2
  mod$H[1, 1, ] <- array(pars[1]^2 + pars[2:(n + 1)]^2, c(1, 1, n))
  mod$Q[1, 1, ] <- array(pars[n + 2]^2 + pars[(n + 3):(2*n + 2)]^2, c(1, 1, n))
  kfs <- KFS(mod, smoothing = "disturbance", simplify = FALSE)
  h <- mod$H[1, 1, ]
  e2D <- (as.numeric(kfs$epshat)^2 - (h - as.numeric(kfs$V_eps))) / (h^2)
  r2N <- as.numeric(kfs$r)^2 - as.numeric(kfs$N)
  d_eps <- -sum(e2D) * pars[1]
  d_eta <- -sum(r2N) * pars[n + 2]
  v_eps <- -e2D * pars[2:(n + 1)] +
    l1 * sign(pars[2:(n + 1)])
  v_eta <- -r2N[2:(n + 1)] * pars[(n + 3):(2*n + 2)] +
    l2 * sign(pars[(n + 3):(2*n + 2)])
  
  c(d_eps, v_eps, d_eta, v_eta)
}

loglik <- function(pars, mod) {
  mod$H[1, 1, ] <- pars[1]^2
  mod$Q[1, 1, ] <- pars[2]^2
  -KFS(mod, filtering = "none", smoothing = "none")$logLik
}

#riscaliamo la serie:
mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0),
               H = array(0, c(1, 1, n)))

mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod, method = "L-BFGS-B")
y1<-y[1]
risc<-abs(mle$par[1])
y<-(y-y[1])*5/abs(mle$par[1])
mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
               H = array(0, c(1, 1, n)))
#BO:
#funzione che definisce valori iniziali e dominio
fun_iniz_dom<-function(y){
  mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod, method = "L-BFGS-B")
  mod_mle <- mod
  mod_mle$H[1, 1, ] <- mle$par[1]^2
  mod_mle$Q[1, 1, ] <- mle$par[2]^2
  kfs_dist <- KFS(mod_mle, smoothing = "disturbance")
  ineta<-abs(kfs_dist$etahat[,1])
  ineps<-abs(kfs_dist$epshat)
  init<-c(sqrt(mod_mle$H[1,1,1])*0.666,ineps,sqrt(mod_mle$Q[1,1,1])*0.666,ineta)*1  
  low<--8*init
  up<-8*init
  dom<-matrix(c(low,up),ncol=2)
  newList <- list("iniziali" = init, "up" = up,"low"=low)
  return(newList)
}

iniz_dom<-fun_iniz_dom(y)#rida i valori iniziale per mettere in input nella prossima funzione

#funzione che calcola il BIC
bic_kalman<-function(iniz_dom,l1,l2){
  k1<-100#servono per definire il bic
  k2<-1
  init<-iniz_dom$iniziali#valori iniziali
  up<-iniz_dom$up#dominio superiore
  low<-iniz_dom$low#dominio inferiore
  
  foo1 <- try(optim(init, obj, grad,mod = mod, l1 = l1, l2 = l2,method = c("BFGS"),control = list(maxit = 50*2*n,factr=1e-10)),silent=TRUE)
  if (class(foo1) == "try-error"){
    foo1 <- try(optim(init, obj, grad,mod = mod, l1 = l1, l2 = l2,lower=low,upper=up,method = c("L-BFGS-B"),control = list(maxit = 50*2*n,factr=1e-10)))
  }
  foo1$loglik <- make_mod(foo1$par, mod, TRUE)$kfs$logLik
  x<-foo1
  bic<--2*x$loglik + log(n)*k2*sum(abs(x$par[1:n+1])>=abs(x$par[1]/k1))+log(n)*k2*sum(abs(x$par[(n+2):(2*n+2)])>=abs(x$par[(n+2)]/k1))
  foo1$bic<-bic
  return(foo1)
}

#funzione di acquisizione
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
  
  # optimize the acquisition function
  fn.scale <- -1 # maximization of the acquisition function
  x0 <- NULL
  for( i in 1:ncol(design) )
    x0 <- c( x0, runif(1,lower[i],upper[i]) )
  
  opt.res <- optim( par=x0, fn=acq.fun, gr=NULL,
                    gp=gp.fit, UCB=F, # <-- BIC is minimized
                    lower=lower, upper=upper, control=list(trace=0,fnscale=fn.scale) )
  
  design <- rbind( design, opt.res$par )
  output<-bic_kalman( iniz_dom=iniz_dom, l1=opt.res$par[1], l2=opt.res$par[2] )
  tot_foo<-c(tot_foo,output)
  y.bic <- c( y.bic,output$bic) 
  
  tmp <- difftime( Sys.time(), tmp, units="secs" )
  bo.iter.elapsed <- c( bo.iter.elapsed, tmp )
  print(nrow(design))
}
cat("> BO time:",bo.elapsed+sum(bo.iter.elapsed),"seconds\n")

setDefaultCluster(cl=NULL)
stopCluster(cl)

tot_par<-tot_foo[which(names(tot_foo)=="par")]
min<-which.min(y.bic)
opt_par<-as.numeric(unlist(tot_par[6]))

#eseguo il kalman smoother e riscalo i risultati
stima_lode <- make_mod(opt_par, mod)$kfs$alphahat[, 1]
stima_lode <- stima_lode*risc/5+y1

#guardo gli outlier individuati
w_ao<-which(abs(opt_par[2:(n+1)])>abs(opt_par[1]/100))
w_level<-which(abs(opt_par[(n+3):(2*n+2)])>abs(opt_par[n+2]/100))
cat("Anno in cui è avvenuto il cambio di livello: ",time(Nile)[w_level+1]) 

#plot dei risultati
ylode<-ts(stima_lode,start = time(Nile)[1],frequency = 1)
plot(Nile,cex.lab=1.3,cex.main=1.5,main="Nilo",xlab="Anno",ylab="Portata fiume")
lines(ylode,col="blue",lwd=2)

