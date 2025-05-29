####################################################################################################


samplealpha <- function(alpha,beta,N,Tt,nj,c1,d1,x,ff){


  n=ncol(x)

  #FF é feito por nós; Devemos alcançar uma taxa de aceitação dee 40%;


  alphaprop<-rgamma(1,alpha*ff, rate=ff)


  ## Temp e temp1 parecem desnecessários

  #temp=array(NA,dim=c(n,1))
  #temp1=array(NA,dim=c(n,1))

  #Tt é auxiliar;
  #c1, d1 são hiperparâmetros;

  # dados <- na.ommit*(drop(data))
  #Transformar x em data

  palpha=(sum(nj)+c1-1)*log(alpha)+alpha*sum(log(x),na.rm = T)-beta*( sum(x^alpha,na.rm = T)+sum(t(N)*(Tt^alpha)) ) -d1*alpha
  palphaprop=(sum(nj)+c1-1)*log(alphaprop)+alphaprop*sum(log(x),na.rm = T)-beta*( sum(x^alphaprop,na.rm = T)+sum(t(N)*(Tt^alphaprop)) )-d1*alphaprop

  #Tirar o 0.0??01

  logprob<-palphaprop+log(dgamma(alpha,alphaprop*ff, rate=ff)+0.0000001)-(palpha + log(dgamma(alphaprop,alpha*ff, rate=ff)+0.0000001))

  probac<-min(c(1,exp(logprob)))

  u<-runif(1)

  if(u<probac){
    res<-alphaprop
    rejei=1


  }
  else{
    res<-alpha
    rejei=0
  }

  res=list(res,rejei)
  res

}


###################################################################################################



samplebeta <- function(alpha,beta,N,Tt,nj,c,d,x,ff)
{

  n=ncol(x)

  betaprop<-rgamma(1,shape=beta*ff, rate=ff)



  temp=array(NA,dim=c(n,1))
  temp1=array(NA,dim=c(n,1))


  pbeta<-(sum(nj)+c-1)*log(beta)-beta*( sum(x^alpha,na.rm = T)+d+sum(t(N)*(Tt^alpha)) )
  pbetaprop<-(sum(nj)+c-1)*log(betaprop)-betaprop*( sum(x^alpha,na.rm = T)+d+sum(t(N)*(Tt^alpha)) )


  logprob<-pbetaprop+log(dgamma(beta,shape=betaprop*ff, rate=ff)+0.0000001)-(pbeta + log(dgamma(betaprop,shape=beta*ff, rate=ff)+0.0000001))

  probac<-min(c(1,exp(logprob)))

  u<-runif(1)

  if(u<probac){
    res<-betaprop
    rejei=1


  }
  else{
    res<-beta
    rejei=0
  }

  res=list(res,rejei)

}

####################################################################################################

sampleb <- function(W,v,b,loca,ab,bb,X,Psi,u1){

  bprop=rgamma(1,shape=b*u1, rate = u1)

  SSigprop=gSigma(bprop,v,loca)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }

  SSig=gSigma(b,v,loca)
  logp=-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(b)-bb*b

  SSigprop=gSigma(bprop,v,loca)
  logpprop=-0.5*t(W-X%*%Psi)%*%solve(SSigprop)%*%(W-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop

  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1)+1e-17)-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)+1e-17))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop

    rejei=1


  }else{

    bprox=b

    rejei=0;

  }



  res=list(bprox,rejei)
  res
}

##################################################################

#def ? loca

gSigma <- function(b,v,loca){
  n<-nrow(loca)
  R<-exp(-b*(as.matrix(dist(loca))))
  mat<-v*R
  mat
}

###################################################################################################
#usu?rio colocaa ab e bb

samplev <- function(W,v,b,loca,ab,bb,X,Psi,u1){

  vprop=rgamma(1,shape=v*u1, rate = u1)

  SSigprop=gSigma(b,vprop,loca)

  if((det(SSigprop)==0)|(vprop< 0.005)){
    return(list(v,0))
  }

  SSig=gSigma(b,v,loca)
  logp=-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(v)-bb*v

  SSigprop=gSigma(b,vprop,loca)
  logpprop=-0.5*t(W-X%*%Psi)%*%solve(SSigprop)%*%(W-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(vprop)-bb*vprop

  logprob=logpprop+log(dgamma(v,shape=vprop*u1,rate=u1)+1e-17)-(logp+log(dgamma(vprop,shape=v*u1,rate=u1)+1e-17))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    vprox=vprop

    rejei=1


  }else{

    vprox=v

    rejei=0;

  }



  res=list(vprox,rejei)
  res
}

####################################################################################################

sampleW <-function(W,loca,X,Psi,b,v,nj,N,u1){

  n=nrow(W)
  Wprop=mvrnorm(1,W,u1*diag(1,n))
  SSig=gSigma(b,v,loca)

  postW=sum(as.matrix(nj)*W)-sum(exp(W))+sum(t(N)*W)-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)
  postWprop=sum(as.matrix(nj)*Wprop)-sum(exp(Wprop))+sum(t(N)*Wprop)-0.5*t(Wprop-X%*%Psi)%*%solve(SSig)%*%(Wprop-X%*%Psi)

  prob=min(exp((postWprop)-(postW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=Wprop

    rejei=1


  }else{

    Wprox=W
    rejei=0
  }

  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}


###################################################################################################

samplePsi<-function(W,X,Psi,V,M,u1,loca,b,v){

  n=nrow(Psi)
  Psiprop=mvrnorm(1,Psi,u1*solve(t(X)%*%X))
  SSig=gSigma(b,v,loca)

  postPsi=-0.5*t(Psi-M)%*%solve(V)%*%(Psi-M)-0.5*t(W-X%*%Psi)%*%solve(SSig)%*%(W-X%*%Psi)

  postPsiprop=-0.5*t(Psiprop-M)%*%solve(V)%*%(Psiprop-M)-0.5*t(W-X%*%Psiprop)%*%solve(SSig)%*%(W-X%*%Psiprop)


  prob=min(exp((postPsiprop)-(postPsi)),1)


  u=runif(1,0,1)

  if(u<prob){

    Psiprox=Psiprop

    rejei=1

  }else{

    Psiprox=Psi
    rejei=0
  }





  res=as.matrix(Psiprox)
  res=list(res,rejei)
  res

}




####################################################################################################


matrixtransf <- function(matrixuser,limuser){

  matrixuser <- data.matrix(matrixuser,rownames.force = NA)

  maior <- apply(matrixuser,2,function(x) length(x[x>limuser]))
  maior <- max(maior)


  ajustado <- c()


  for(i in 1:ncol(matrixuser)){
    linhaprov <- c()

    linha <- matrixuser[,i]

    ## pensar em outro apply
    for(j in 1:nrow(matrixuser)){
      if(linha[j]>limuser){

        linhaprov <- c(linhaprov,j)

      }

    }
    linhaprov

    if(length(linhaprov)<maior){

      linhaprov <- c(linhaprov,rep(0,(maior-length(linhaprov))))

      ajustado<- rbind(ajustado,linhaprov)

    }else{

      ajustado <- rbind(ajustado,linhaprov)
    }

  }

  ajustado <- t(ajustado)
  ajustado[ajustado==0] <- NA

  return(ajustado)
}
