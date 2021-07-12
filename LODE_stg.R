# Funzione che implementa il metodo LODE, modellando la serie con un 
# random walk e ricercando AO e LS.

# Iperettangolo dei lambda
# lower <- 0.03
# upper <- 1.3
# Numero iterazioni della BO
# n.trials <- 15
# Numero di valori iniziali da generare con LHS
# vallhs <- 7
# Decidere se effettuare o meno il riscalamento automatico della serie
# Riscalamento = T
# tsy è la serie storica e deve essere una variabile di classe ts
LODE_stg<-function(tsy,lower=0.03,upper=1.3,n.trials=15,vallhs=7,riscalamento=T){
  
  if(vallhs>n.trials){
    stop("il numero di valori iniziali da generare con LHS deve essere minore del numero totale di iterazioni")
  }
  if(class(tsy)!="ts"){
    stop("la classe della serie storica non è ts")
  }
  
  library(KFAS)
  library(dplyr)
  library(DiceKriging)
  library(doParallel)
  library(lhs)
  y<-as.numeric(tsy)
  timey<-time(tsy)
  n<-length(y)
  per<-frequency(tsy) #periodicita

  #Funzione obiettivo da minimizzare
  obj <- function(pars, mod, l4){
    param<-rep(0,n)
    for(i in (per+1):(n-1*per)){
      param[i:(i+per-1)]<-param[i:(i+per-1)]+(rep(pars[i+3],per)^2)
    }
    mod$H[1, 1, ] <- array(pars[1]^2, c(1, 1, n))
    mod$Q[1, 1, ] <- array(pars[2]^2,c(1, 1, n))
    mod$Q[2, 2, ] <- array(pars[3]^2+param,c(1,1,n))
    -KFS(mod, filtering = "none", smoothing = "none")$logLik +
      l4 * sum(abs(pars[4:(n+3)]))#param
    
  }
  
  #Funzione che calcola lo smoother con le varianze aggiuntive
  make_mod <- function(pars, mod, kfs = TRUE) {
    param<-rep(0,n)
    for(i in (per+1):(n-1*per)){
      param[i:(i+per-1)]<-param[i:(i+per-1)]+(rep(pars[i+3],per)^2)
    }
    mod$H[1, 1, ] <- array(pars[1]^2, c(1, 1, n))
    mod$Q[1, 1, ] <- array(pars[2]^2,c(1, 1, n))
    mod$Q[2, 2, ] <- array(pars[3]^2+param,c(1,1,n))
    if (kfs) return(list(mod = mod, kfs = KFS(mod, smoothing = "state")))
    list(mod = mod, kfs = NULL)
  }
  
  #Gradiente della funzione obiettivo
  grad <- function(pars, mod, l4) {  
    param<-rep(0,n)
    for(i in (per+1):(n-1*per)){
      param[i:(i+per-1)]<-param[i:(i+per-1)]+(rep(pars[i+3],per)^2)
    }
    mod$H[1, 1, ] <- array(pars[1]^2, c(1, 1, n))
    mod$Q[1, 1, ] <- array(pars[2]^2,c(1, 1, n))
    mod$Q[2, 2, ] <- array(pars[3]^2+param,c(1,1,n))
    kfs <- KFS(mod, smoothing = "disturbance", simplify = FALSE)
    h <- mod$H[1, 1, ]
    e2D <- (as.numeric(kfs$epshat)^2 - (h - as.numeric(kfs$V_eps))) / (h^2)
    r2N <- (kfs$r)^2 - apply(kfs$N,3,diag)
    d_eps <- -sum(e2D) * pars[1]
    d_eta <- -sum(r2N[1,]) * pars[2]
    d_stg <- -sum(r2N[2,]) * pars[3]
    v_stg <- c(rep(0,per),-(r2N[2,(per+2):(n + 1-per)]+r2N[2,(per+3):(n + 1-per+1)]+r2N[2,(per+4):(n + 1-per+2)]+r2N[2,(per+5):(n + 1-per+3)])
               *pars[(per+4):(n + 3-per)],rep(0,per))+ l4 * sign(pars[(4):(n + 3)])
    c(d_eps, d_eta,d_stg, v_stg)
  }
  
  #Log(verosimiglianza)
  loglik <- function(pars, mod) {
    mod$H[1, 1, ] <- pars[1]^2
    mod$Q[1, 1, ] <- pars[2]^2
    mod$Q[2, 2, ] <- pars[3]^2
    -KFS(mod, filtering = "none", smoothing = "none")$logLik
  }
  
  if(riscalamento==T){
    #Riscaliamo la serie
    mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0)+SSMseasonal(per,100,"dummy"),
                   H = array(0, c(1, 1, n)))
    
    
    mle <- optim(c(sd(diff(y,lag=per))/100, sd(diff(y,lag=per))/100, sd(diff(y,lag=per)))/100, loglik, mod = mod, method = "L-BFGS-B")
    y1<-y[1]
    risc<-abs(mle$par[1])
    y<-(y-y[1])*5/abs(mle$par[1])
  }
    mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0)+SSMseasonal(per,100,"dummy"),
                   H = array(0, c(1, 1, n)))
  
  #BO
  #Funzione che definisce valori iniziali e dominio
    fun_iniz_dom<-function(y){
      mle <- optim(c(sd(diff(y,lag=per))/100, sd(diff(y,lag=per))/100, sd(diff(y,lag=per)))/100, loglik, mod = mod, method = "L-BFGS-B")
      mod_mle <- mod
      mod_mle$H[1, 1, ] <- mle$par[1]^2
      mod_mle$Q[1, 1, ] <- mle$par[2]^2
      mod_mle$Q[2, 2, ] <- mle$par[3]^2
      #iniz con smoother
      kfs_dist <- KFS(mod_mle, smoothing = "disturbance")
      ineta<-abs(kfs_dist$etahat[,2])
      init<-c(sqrt(mod_mle$H[1,1,1]),sqrt(mod_mle$Q[1,1,1]),sqrt(mod_mle$Q[2,2,1])*0.666,ineta)*1  
      # ---- cv with respect to F
      low<--8*init
      up<-8*init
      dom<-matrix(c(low,up),ncol=2)
      newList <- list("iniziali" = init, "up" = up,"low"=low)
      return(newList)
    }
    
    iniz_dom<-fun_iniz_dom(y)#rida i valori iniziale per mettere in input nella prossima funzione
    
    #Funzione che calcola il BIC
    bic_kalman<-function(iniz_dom,l4){
      k1<-100#servono per definire il bic
      k2<-1
      init<-iniz_dom$iniziali#valori iniziali
      up<-iniz_dom$up#dominio superiore
      low<-iniz_dom$low#dominio inferiore
      
      foo1 <- try(optim(init, obj, grad,mod = mod, l4= l4, method = c("BFGS"),control = list(maxit = 50*2*n)),silent=TRUE)
      if (class(foo1) == "try-error"){
        foo1 <- try(optim(init, obj, grad,mod = mod, l4 = l4,lower=low,upper=up,method = c("L-BFGS-B"),control = list(maxit = 50*2*n)))
      }
      foo1$loglik <- make_mod(foo1$par, mod, TRUE)$kfs$logLik
      x<-foo1
      bic<--2*x$loglik +log(n)*k2*sum(abs(x$par[(3):(n+3)])>=abs(x$par[(3)]/k1))
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
  design[] <- design[] * ( upper[] - lower[] ) + lower[]
  
  bo.elapsed <- Sys.time()
  tot_foo <- foreach( i = 1:n0, .packages="KFAS", .combine="c" ) %dopar% {
    bic_kalman( iniz_dom=iniz_dom, l4=design[i])
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
                      lower=lower,method="L-BFGS-B", upper=upper, control=list(trace=0,fnscale=fn.scale) )
    
    design <- rbind( design, opt.res$par )
    #Calcolo il BIC nel nuovo punto
    output<-bic_kalman( iniz_dom=iniz_dom, l4=opt.res$par[1] )
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
  if(riscalamento==T){
  y<-y*risc/5+y1
  mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0)+SSMseasonal(per,100,"dummy"),
                   H = array(0, c(1, 1, n)))
  stima_lode <- make_mod(opt_par, mod)$kfs$alphahat
  }else{stima_lode <- make_mod(opt_par, mod)$kfs$alphahat
  }
  stima_lode_tr <- stima_lode[, 1]
  
  
  #guardo gli outlier individuati
  w_stg<-which(abs(opt_par[(4):(n+3)])>=abs(opt_par[(3)]/100))
  cat("Istanti in cui sono avvenuti cambi di stagionalità: ",timey[w_stg+1],"\n") 

  #plot dei risultati
  ylode<-ts(stima_lode_tr,start = timey[1],frequency = per)
  plot(tsy,cex.lab=1.3,cex.main=1.5,main="Serie analizzata e stima trend LODE",xlab="Tempo",ylab="Valore serie storica")
  lines(ylode,col="blue",lwd=2)
  return(list(trend=ylode,sc=w_stg,componenti=stima_lode))
}



#Esempio veloce di una sua applicazione
set.seed(1)
n<-50
sig_eps <- sig_eta <- 1
sig_stg<-0.6
per<-4
s<-numeric()
s[1:(per-1)]<-c(runif(1,-5,5),runif(1,-5,5),runif(1,-5,5))
for(i in per:n){
  s[i]<--sum(s[(i-per+1):(i-1)])+ifelse(i==25|i==26|i==27|i==28,runif(1,5,10)*sig_stg*sample(c(1,-1),1),rnorm(1,0,sig_stg))
}
s<-s*sig_eps
mu <- cumsum(rnorm(n, sd = sig_eta))+runif(1,-100*sig_eps,100*sig_eps)
y <- mu + rnorm(n, sd = sig_eps)+ s
tsy<-ts(y,start = 1970,frequency = 4)
plot(tsy,type="l")
x<-LODE_stg(tsy)
x$trend #stima del trend con metodo LODE
x$sc #posizione delle varianze aggiuntive della stagionalità diverse da zero
x$componenti #stima della componenti del modello




  