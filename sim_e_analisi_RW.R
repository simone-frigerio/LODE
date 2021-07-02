#Simulazione di serie modellate con un RW ed analisi con LODE, auxres e KFS
library(KFAS)
library(dplyr)
library(DiceKriging)
library(doParallel)

#simulazione ed analisi LODE ----
opts1 <- vector(mode = "list")
val1<-vector(mode = "list")
opts2 <- vector(mode = "list")
val2<-vector(mode = "list")

tsy<-matrix(rep(0,100*100),nrow = 100)
tsmu<-matrix(rep(0,100*100),nrow = 100)
jump<-numeric()
eps<-numeric()
sigeps<-numeric()
sigeta<-numeric()
serie1<-numeric()
y1<-numeric()
risc<-numeric()

#parametri della BO
lower <- c(0.1,0.1)
upper <- c(2,2)
n.trials <- 50
vallhs<-17


for(j in 1:1){
  set.seed(17052021+j)
  if(j==34|j==102|j==236){#cambio alcune serie perchè danno problemi con KFS
    set.seed(j)
  }
  n<-100
  sig_eps <- runif(1,1,100)
  sig_eta <- runif(1,1,4)*sig_eps/2
  eps_sd_fct <- runif(1,2,6)
  jmp_size <- runif(1,2,6)
  disturbiao<-matrix(c(rnorm(n,0,sd=sig_eps),rnorm(n,0,sd=eps_sd_fct*sig_eps)),ncol=2)
  eps<-apply(disturbiao, 1, sample,prob=c(0.9,0.1),size=1)
  disturbij<-matrix(c(rnorm(n,0,sd=sig_eta),rnorm(n,0,sd=jmp_size*sig_eta)),ncol=2)
  eta<-apply(disturbij, 1, sample,prob=c(0.9,0.1),size=1)
  
  mu <- cumsum(eta)+runif(1,-100*sig_eps,100*sig_eps)
  y <- mu + eps
  
  tsy[j,]<-y
  tsmu[j,]<-mu
  
  plot(y, type = "l",main=paste(j))
  lines(mu, col = "red")
 

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
  
  print(j)
  mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0),
                 H = array(0, c(1, 1, n)))
  
  mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod, method = "L-BFGS-B")
  y1[j]<-y[1]
  risc[j]<-abs(mle$par[1])
  y<-(y-y[1])*5/abs(mle$par[1])
  mod <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                 H = array(0, c(1, 1, n)))
  mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod, method = "L-BFGS-B")

  #BO
  #funzione che definisce valori iniziali e dominio
  fun_iniz_dom<-function(y){
    mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod, method = "L-BFGS-B")
    mod_mle <- mod
    mod_mle$H[1, 1, ] <- mle$par[1]^2
    mod_mle$Q[1, 1, ] <- mle$par[2]^2
    #iniz con smoother
    kfs_dist <- KFS(mod_mle, smoothing = "disturbance")
    ineta<-abs(kfs_dist$etahat[,1])
    ineps<-abs(kfs_dist$epshat)
    init<-c(sqrt(mod_mle$H[1,1,1])*0.666,ineps,sqrt(mod_mle$Q[1,1,1])*0.666,ineta)*1  
    # ---- cv with respect to F
    low<--8*init
    up<-8*init
    dom<-matrix(c(low,up),ncol=2)
    newList <- list("iniziali" = init, "up" = up,"low"=low)
    return(newList)
  }
  
  iniz_dom<-fun_iniz_dom(y)#rida i valori iniziale per mettere in input nella prossima funzione
  
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
    
    return(bic)
  }
  
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
  
  seed <- 1
  
  kernel <- "gauss"
  acquisition <- "CB" 
  

  cat("\n> ***** Starting BO *****\n")
  n0 <- round(vallhs) 
  
  set.seed(seed)
  design <- lhs::maximinLHS(n=n0,k=length(lower))
  for( i in 1:ncol(design) ){
    design[,i] <- design[,i] * ( upper[i] - lower[i] ) + lower[i]
  }
  
  bo.elapsed <- Sys.time()
  y.bic <- foreach( i = 1:n0, .packages="KFAS", .combine="c" ) %dopar% {
    bic_kalman( iniz_dom=iniz_dom, l1=design[i,1], l2=design[i,2] )
  }
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
    fn.scale <- -1 # maximization of the acquisition function!
    x0 <- NULL
    for( i in 1:ncol(design) )
      x0 <- c( x0, runif(1,lower[i],upper[i]) )
    
    opt.res <- optim( par=x0, fn=acq.fun, gr=NULL,
                      gp=gp.fit, UCB=F, #BIC is minimized
                      lower=lower, upper=upper, control=list(trace=0,fnscale=fn.scale) )
    
    design <- rbind( design, opt.res$par )
    y.bic <- c( y.bic, bic_kalman( iniz_dom=iniz_dom, l1=opt.res$par[1], l2=opt.res$par[2] )) 
    
    tmp <- difftime( Sys.time(), tmp, units="secs" )
    bo.iter.elapsed <- c( bo.iter.elapsed, tmp )
    print(nrow(design))
  }
  cat("> BO time:",bo.elapsed+sum(bo.iter.elapsed),"seconds\n")
  
  setDefaultCluster(cl=NULL)
  stopCluster(cl)
  
  bo.results <- cbind(design,y.bic)
  colnames(bo.results) <- c("l1","l2","bic")
  
  #50 iterazioni
  ar<-arrange(as.data.frame(bo.results),by=bic)
  l1opt<-ar[1,1]#guardo i parametri migliori dopo 50 iterazioni
  l2opt<-ar[1,2]
  foo1 <- try(optim(iniz_dom$iniziali, obj, grad,mod = mod, l1 = l1opt, l2 = l2opt,method = c("BFGS"),control = list(maxit = 50*2*n,factr=1e-10)),silent=TRUE)
  if (class(foo1) == "try-error"){
    foo1 <- try(optim(iniz_dom$iniziali, obj, grad,mod = mod, l1 = l1opt, l2 = l2opt,lower=iniz_dom$low,upper=iniz_dom$up,method = c("L-BFGS-B"),control = list(maxit = 50*2*n,factr=1e-10)))
  }
  foo1$loglik <- make_mod(foo1$par, mod, TRUE)$kfs$logLik
  opts1[[paste0("j=",j)]] <-  foo1
  val1[[paste0("j=",j)]] <-  foo1$value
  serie1[j]<-length(val1)
  
  #30 iterazioni
  ar<-arrange(as.data.frame(bo.results[1:30,]),by=bic)
  l1opt<-ar[1,1]#guardo i parametri migliori dopo 30 iterazioni
  l2opt<-ar[1,2]
  foo1 <- try(optim(iniz_dom$iniziali, obj, grad,mod = mod, l1 = l1opt, l2 = l2opt,method = c("BFGS"),control = list(maxit = 50*2*n,factr=1e-10)),silent=TRUE)
  if (class(foo1) == "try-error"){
    foo1 <- try(optim(iniz_dom$iniziali, obj, grad,mod = mod, l1 = l1opt, l2 = l2opt,lower=iniz_dom$low,upper=iniz_dom$up,method = c("L-BFGS-B"),control = list(maxit = 50*2*n,factr=1e-10)))
  }
  foo1$loglik <- make_mod(foo1$par, mod, TRUE)$kfs$logLik
  opts2[[paste0("j=",j)]] <-  foo1
  val2[[paste0("j=",j)]] <-  foo1$value

}



#analisi auxres ed estrazone risultati------
m<-1

#inizializzazione delle varie variabili
bst_lasso_mu1<-matrix(rep(0,100*m),nrow=m)
bst_lasso_mu2<-matrix(rep(0,100*m),nrow=m)
kfs_mle<-matrix(rep(0,100*m),nrow=m)
muauxres<-matrix(rep(0,100*m),nrow=m)
erbest1<-numeric()
erbest2<-numeric()
erauxres<-numeric()
err<-numeric()
bic1<-numeric()
pv<-0.995
sw_lev<-numeric()
sw_y<-numeric()
c<-serie1
ii<-list()
ii[[1]]<-1:c[1]
for(i in 2:m){
  ii[[i]]<-(c[i-1]+1):c[i]
}
#parametri che entrano in gioco nel calcolo del BIC
k1<-100
k2<-1

#ciclo sulle varie serie simulate
for(i in 1:m){
  ind<-ii[[i]]
  print(i)
  y<-tsy[i,]
  mu<-tsmu[i,]
  
  #KFS
  mod1 <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0),
                  H = array(0, c(1, 1, n)))
  
  #faccio prima optim e poi fitSSM per comodità nel ricopiare il codice
  mle <- optim(c(sd(y)/4, sd(y)/4), loglik, mod = mod1, method = "L-BFGS-B")
  mod_mle <- mod1
  mod_mle$H[1, 1, ] <- mle$par[1]^2
  mod_mle$Q[1, 1, ] <- mle$par[2]^2
  updtfn <- function(pars, model) {
    model$H[1, 1, ] <- exp(pars[1])
    model$Q[1, 1, ] <- exp(pars[2])
    model
  }
  init <- c(log(mod_mle$H[1, 1, 1]),log(mod_mle$Q[1, 1, 1])) # initial values
  fit1 <- fitSSM(mod1, init, updtfn)
  ssmo1 <- KFS(fit1$model, smoothing="state")
  #stime del trend ed errore KFS
  kfs_mle[i,] <- ssmo1$alphahat[,1]
  err[i]<-mean((tsmu[i,] - kfs_mle[i,])^2)
  
  #LODE50
  bst_lasso_mu1[i,] <-with(opts1[[i]], make_mod(par, mod1)$kfs$alphahat[, 1])
  erbest1[i]<-mean((tsmu[i,] - bst_lasso_mu1[i,])^2)
  #LODE30
  bst_lasso_mu2[i,] <-with(opts2[[i]], make_mod(par, mod1)$kfs$alphahat[, 1])
  erbest2[i]<-mean((tsmu[i,] - bst_lasso_mu2[i,])^2)
  
  #auxres
  dsmo1 <- KFS(fit1$model, smoothing = "disturbance")
  auxres_level <- rstandard(dsmo1, "state")
  smo1 <- KFS(fit1$model)
  auxres_y <- rstandard(smo1, "pearson")
  
  #test SW
  sw_lev[i]<-shapiro.test(auxres_level)$p.value
  sw_y[i]<-shapiro.test(auxres_y)$p.value
  
  #AO
  ndx_ao<- which(abs(auxres_y) > qnorm(pv)) 
  if (length(ndx_ao)==0){
    AO=AO<-matrix(0,nrow=length(y),ncol=1)
  } else {
    AO<-matrix(0,nrow=length(y),ncol=length(ndx_ao))
    for(l in 1:length(ndx_ao)){
      AO[ndx_ao[l],l]<-1
    }
  }
  
  #level-shift
  ndx <- which(abs(auxres_level[1:(n-1)]) > qnorm(pv))
  
  if (length(ndx)==0){
    shift=shift<-matrix(0,nrow=length(y),ncol=1)
    max<-numeric()
  } else {
    #selezioniamo solo uno per gruppo di salti
    spi<-split(ndx, cumsum(c(1, diff(ndx) != 1)))
    max<-numeric()
    mm<-0
    for(g in spi){
      mm<-mm+1
      max[mm]<-which.max(abs(auxres_level[g]))+g[1]-1
    }
    shift<-matrix(0,nrow=length(y),ncol=length(max))
    for(l in 1:length(max)){
      shift[(max[l]+1):length(y),l]<-1
    }
  }
  
  
  #Stimo modello con le dummy
  if(length(ndx_ao)==0){   
    mod2 <- SSModel(y ~ shift+SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                    H = array(0, c(1, 1, n)))
  }else if(length(ndx)==0){
    mod2 <- SSModel(y ~ AO+SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                    H = array(0, c(1, 1, n)))
  }else{ 
    mod2 <- SSModel(y ~ AO+shift+SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                    H = array(0, c(1, 1, n)))
  }
  
  fit2 <- fitSSM(model=mod2, inits=fit1$optim.out$par, updatefn=updtfn)
  fit2$optim.out$convergence
  
  options(warn=2)
  ssmo2 <- try(KFS(fit2$model, smoothing="state"))
  if(class(ssmo2)=="try-error"){
    if(length(ndx_ao)==0){   
      mod2 <- SSModel(y ~ shift+SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                      H = array(0, c(1, 1, n)))
      diag(mod2$P1inf)[1:length(max)]<-0
      diag(mod2$P1)[1:length(max)]<-var(y)*100
    }else if(length(ndx)==0){
      mod2 <- SSModel(y ~ AO+SSMtrend(1, array(0, c(1, 1, n)),  a1 = y[1], P1 = var(y), P1inf = 0),
                      H = array(0, c(1, 1, n)))
      diag(mod2$P1inf)[1:length(ndx_ao)]<-0
      diag(mod2$P1)[1:length(ndx_ao)]<-var(y)*100
    }else{ 
      mod2 <- SSModel(y ~ AO+shift+SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0),
                      H = array(0, c(1, 1, n)))
      diag(mod2$P1inf)[(1):(length(ndx_ao)+length(max))]<-0
      diag(mod2$P1)[(1):(length(ndx_ao)+length(max))]<-var(y)*100
    }
    fit2 <- fitSSM(model=mod2, inits=fit1$optim.out$par, updatefn=updtfn)
    options(warn=1)
    ssmo2 <- try(KFS(fit2$model, smoothing="state"))
  }
  
  options(warn=1)
  #Estrazione delle componenti
  level <- ssmo2$alphahat[,"level"]
  if(length(max)==1){
    shift <- as.numeric(ssmo2$alphahat[,(length(ndx_ao)+1):(length(ndx_ao)+length(max))]) * t(as.matrix(fit2$model$Z[1,(length(ndx_ao)+1):(length(ndx_ao)+length(max)),]))
    shiftot<-as.numeric(shift) 
  }else if(length(ndx)==0){
    shiftot<-rep(0,n)
  }else{
    shift <- as.matrix(ssmo2$alphahat[,(length(ndx_ao)+1):(length(ndx_ao)+length(max))]) * t(as.matrix(fit2$model$Z[1,(length(ndx_ao)+1):(length(ndx_ao)+length(max)),]))
    shiftot<-apply(shift, 1, sum) 
  }
  
  #Grafico della serie originale e delle componenti stimate
  erauxres[i]<-mean((tsmu[i,]-(level + shiftot))^2)
  muauxres[i,]<-level+shiftot
  
  
  plot(tsy[i,],col="green",type="l",main=paste0("Serie ",i),xlab="Tempo",ylab="Valore serie e stime")
  lines(tsmu[i,],type="l",col="black",lwd=1)
  lines(bst_lasso_mu1[i,],col="red",lwd=2)
  lines(muauxres[i,] , col="blue",lwd=2)
  lines(kfs_mle[i,],col="orange",lwd=1.5)

}

#MSE-----
#tutte
eauxres<-erauxres[]/err[] #auxres
elode50<-erbest1[]/err[] #LODE50
elode30<-erbest2[]/err[] #LODE30
mean(eauxres[])
mean(elode50[])
mean(elode30[])

#MSE<0.95
kk1<-which(eauxres<0.95)
kk2<-which(elode50<0.95)
k<-kk2[which(kk2%in%kk1)]
mean(eauxres[k])
mean(elode50[k])

#SW<5%
wsw_lev<-which(sw_lev[]<0.05)
wsw_y<-which(sw_y[]<0.05)
wsw<-unique(c(wsw_lev,wsw_y))
mean(eauxres[wsw])
mean(elode50[wsw])
mean(elode30[wsw])




