library(KFAS)

soglia<-qnorm(0.995)

#plot
plot(Nile,main="Nilo",lwd=2,xlab="Anno",ylab="Portata fiume",cex.lab=1.3,cex.main=1.5)

y<-as.numeric(Nile)
n<-length(y)

#modello
mod1 <- SSModel(y ~ SSMtrend(1, array(0, c(1, 1, n)), a1 = y[1], P1 = var(y), P1inf = 0),
                H = array(0, c(1, 1, n)))

updtfn <- function(pars, model) {
  model$H[1, 1, ] <- exp(pars[1]) 
  model$Q[1, 1, ] <- exp(pars[2]) 
  model
}
init <- c(log(var(y)/4),log(var(y)/4)) # initial values
fit1 <- fitSSM(mod1, init, updtfn)

# State smoother
smo1 <- KFS(fit1$model)
# Disturbance smoother
dsmo1 <- KFS(fit1$model, smoothing = "disturbance")

# Analisi dei residui ausiliari
auxres_level <- rstandard(dsmo1, "state")
plot(auxres_level,main="Residui ausiliari del trend")
abline(h=c(-soglia,+soglia),lty=3,lwd=2,col="grey")
auxres_y <- rstandard(smo1, "pearson")
plot(auxres_y,main="Residui ausiliari additivi")
abline(h=c(-soglia,+soglia),lty=3,lwd=2,col="grey")

#AO
ndx_ao<- which(abs(auxres_y) > soglia)
if (length(ndx_ao)==0){
  AO=AO<-matrix(0,nrow=length(y),ncol=1)
} else {
  AO<-matrix(0,nrow=length(y),ncol=length(ndx_ao))
  for(l in 1:length(ndx_ao)){
    AO[ndx_ao[l],l]<-1
  }
}

#Cambio di livello
ndx <- which(abs(auxres_level[1:(n-1)]) > soglia) # salti 
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


#Ristimo modello con gli outlier
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

#Non uso le condizioni iniziali diffuse se ridà il warning
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

#Grafico della serie originale e del trend stimato
stima_auxres<-level+shiftot

yauxres<-ts(stima_auxres,start = time(Nile)[1],frequency = 1)
plot(Nile,cex.lab=1.3,cex.main=1.5,main="Nilo",xlab="Anno",ylab="Portata fiume")
lines(yauxres,col="green4",lwd=2)

#Guardo gli outlier trovati
cat("Anno in cui è avvenuto il cambio di livello: ",time(Nile)[max+1],"\n") 
cat("Anno in cui è presente un outlier additivo: ",time(Nile)[ndx_ao],"\n") 
