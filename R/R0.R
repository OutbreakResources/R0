R0<-function(yobs,threshold,pi){
  

  
  Like1<-function(R){
    # with or without the conditional
    if (R<=1){
      f<-f_z(z,1,R)
    }else{
      alpha<-R+lamW::lambertW0(-R*exp(-R))
      f<-f_z(z,1,R*exp(-alpha))*exp(-alpha)
    }
    g<-dbinom(Y,matrix(z,nrow(Y),threshold,byrow=TRUE),pi)
    L<- -(sum(log(g %*% f)))
    return(L)
  }
  
  ExpectedR<-1-pi/yobs
  ExpectedR
  # Observed size
  Y<-matrix(yobs,length(yobs),threshold,byrow=FALSE)
  # possible sizes
  z<-matrix(seq(1,threshold),threshold,1)
  
  # L<-Like1(.75) # check function
  res<-optim(median(ExpectedR),Like1, method ="L-BFGS-B", lower = 1e-5)
  
  # 95%CI
  Like1CI<-function(x){
    res<-(Like1(x)-(Like1(res$par)+qchisq(.95, df=1)/2))^2
    return(res)
  }
#   #check
#   p<-sapply(seq(0,4,length=100),Like1CI)  
#   plot(seq(0,4,length=100),p,ylim=c(0,(qchisq(.95, df=1)/2)^2))
  
  lower<-optim(res$par*.9,Like1CI, method ="L-BFGS-B", lower = 1e-6, upper = res$par)
  upper<-optim(res$par*1.1,Like1CI, method ="L-BFGS-B", lower = res$par)
  
  c(mean=res$par,lower=lower$par,upper=upper$par)
}

f_z<-function(z,s,R){ 
  logf<-(z-s-1)*log( s*z)+(z-s)*log(R)+ (-z*R) - lgamma(z-s +1)
  f<-exp(logf)
}