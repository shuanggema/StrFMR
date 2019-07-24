############################################################################
#####    codes for "Structured Analysis of the High-dimensional FMR Model"
#####   
#####
##############################################################################
FMR_In<- function(X,Y,lambda_1,lambda_2){#############dim(X)=n*p,
  n          <- nrow(X)
  p          <- ncol(X)
  
  b=rep(0.1,p)
  beta=as.matrix(b,nc=1)  ##initial guess
  beta0=as.matrix(b/2,nc=1) 
  gamma=as.matrix(b/3,nc=1)
  gamma0=as.matrix(b/2,nc=1)
  f1=0.3
  f10=0.4
  w1=w2=rep(0.3,n)
  w10=w20=rep(0.4,n)
  r10=2
  r11=1
  r20=2
  r21=1
  tau=2*10^-1
  m=1
  while(abs(sum(w1-w10))>0.01){######E step
    m=m+1
    if (m>10){
      break
    }
    
    w10=w1
    w20=w2
    p12<- (f1/(1-f1))*(r11/r21)*exp(-1/2*(r11*Y-X%*%beta)^2+1/2*(r21*Y-X%*%gamma)^2)
    w1<- 1/(1+1/p12)###n*1   conditional probabilities that individual belongs to the cured or uncured 
    w2<-1/(p12+1)
    f1<- mean(w1)
    diff=1
    u=1
    while(diff>0.1){
      u=u+1
      if (u>200){
        break
      }
      
      r10=r11
      r20=r21
      
      f1<- mean(w1)
      Sb=c()
      Sg=c()
      #############M step#########
      beta0=beta
      gamma0=gamma
      Y_sb=as.vector(sqrt(w1))*Y
      X_sb=as.vector(sqrt(w1))*X
      #################################### for j=1
      r11=(sum(Y_sb*as.vector(X_sb%*%beta))+sqrt(sum(Y_sb*as.vector(X_sb%*%beta))^2+4*sum(Y_sb^2)*sum(w1)))/(2*sum(Y_sb^2))
      sum_sj20=0
      for (k in 2:p){
        sum_sj21=sum_sj20+beta0[k]*sum(X_sb[,1]*X_sb[,k])
        sum_sj20=sum_sj21
      }
      Sb[1]=-r11*sum(Y_sb*X_sb[,1])+sum_sj21
      
      Y_sg=as.vector(sqrt(w2))*Y
      X_sg=as.vector(sqrt(w2))*X
      r21=(sum(Y_sg*as.vector(X_sg%*%gamma))+sqrt(sum(Y_sg*as.vector(X_sg%*%gamma))^2+4*sum(Y_sg^2)*sum(w2)))/(2*sum(Y_sg^2))
      
      sum_sj20=0
      for (k in 2:p){
        sum_sj21=sum_sj20+gamma0[k]*sum(X_sg[,1]*X_sg[,k])
        sum_sj20=sum_sj21
      }
      Sg[1]=-r21*sum(Y_sg*X_sg[,1])+sum_sj21
      beta[1]=0
      gamma[1]=0
      ######beta>0,gamma>0, beta>gamma
      T=exp(-(beta0[1]-gamma0[1])^2/tau)
      l=2*n*lambda_2*T/tau
      
      A_1=-Sb[1]-n*lambda_1*((f1)^0.5)+l*gamma0[1]
      
      A_2=-Sg[1]-n*lambda_1*((1-f1)^0.5)+l*beta0[1]
      A_3=-Sb[1]+n*lambda_1*((f1)^0.5)+l*gamma0[1]
      A_4=-Sg[1]+n*lambda_1*((1-f1)^0.5)+l*beta0[1]
      
      if (A_1/(l+sum(X_sb[,1]^2))>0){  
        
        beta[1]=A_1/(l+sum(X_sb[,1]^2))
      }
      if (A_3/(l+sum(X_sb[,1]^2))<0){  
        
        beta[1]=A_3/(l+sum(X_sb[,1]^2))
      }
      if (A_2/(sum(X_sg[,1]^2)+l)>0){  
        
        gamma[1]=A_2/(sum(X_sg[,1]^2)+l)
      } 
      if (A_4/(sum(X_sg[,1]^2)+l)<0){  
        
        gamma[1]=A_4/(sum(X_sg[,1]^2)+l)
      } 
      ###########################################################  
      #################################### for j=2:(p-1)
      for (j in 2:(p-1)){
        sum_sj10=0
        for (k in 1:(j-1)){
          sum_sj11=sum_sj10+beta[k]*sum(X_sb[,j]*X_sb[,k])
          sum_sj10=sum_sj11
        }
        sum_sj20=0
        for(k in c((j+1):p)){
          sum_sj21=sum_sj20+beta0[k]*sum(X_sb[,j]*X_sb[,k])
          sum_sj20=sum_sj21
        }
        Sb[j]=-r11*sum(Y_sb*X_sb[,j])+sum_sj11+sum_sj21
        
        sum_sj10=0
        for (k in 1:(j-1)){
          sum_sj11=sum_sj10+gamma[k]*sum(X_sg[,j]*X_sg[,k])
          sum_sj10=sum_sj11
        }
        sum_sj20=0
        for(k in c((j+1):p)){
          sum_sj21=sum_sj20+gamma0[k]*sum(X_sg[,j]*X_sg[,k])
          sum_sj20=sum_sj21
        }
        Sg[j]=-r21*sum(Y_sg*X_sg[,j])+sum_sj11+sum_sj21
        beta[j]=0
        gamma[j]=0
        T=exp(-(beta0[j]-gamma0[j])^2/tau)
        l=2*n*lambda_2*T/tau
        
        A_1=-Sb[j]-n*lambda_1*(f1^0.5)+l*gamma0[j]
        
        A_2=-Sg[j]-n*lambda_1*((1-f1)^0.5)+l*beta0[j]
        A_3=-Sb[j]+n*lambda_1*(f1^0.5)+l*gamma0[j]
        A_4=-Sg[j]+n*lambda_1*((1-f1)^0.5)+l*beta0[j]
        
        
        if (A_1/(l+sum(X_sb[,j]^2))>0){  
          
          beta[j]=A_1/(l+sum(X_sb[,j]^2))
        }
        if (A_3/(l+sum(X_sb[,j]^2))<0){  
          
          beta[j]=A_3/(l+sum(X_sb[,j]^2))
        }
        if (A_2/(sum(X_sg[,j]^2)+l)>0){  
          
          gamma[j]=A_2/(sum(X_sg[,j]^2)+l)
        } 
        if (A_4/(sum(X_sg[,j]^2)+l)<0){  
          
          gamma[j]=A_4/(sum(X_sg[,j]^2)+l)
        } 
      }
      
      
      #################################### for j=p
      
      sum_sj10=0
      for (k in 1:(p-1)){
        sum_sj11=sum_sj10+beta[k]*sum(X_sb[,p]*X_sb[,k])
        sum_sj10=sum_sj11
      }
      
      Sb[p]=-r11*sum(Y_sb*X_sb[,p])+sum_sj11
      
      sum_sj10=0
      for (k in 1:(p-1)){
        sum_sj11=sum_sj10+gamma[k]*sum(X_sg[,p]*X_sg[,k])
        sum_sj10=sum_sj11
      }
      
      Sg[p]=-r21*sum(Y_sg*X_sg[,p])+sum_sj11
      
      beta[p]=0
      gamma[p]=0
      ######beta>0,gamma>0, beta>gamma
      T=exp(-(beta0[p]-gamma0[p])^2/tau)
      l=2*n*lambda_2*T/tau
      
      A_1=-Sb[p]-n*lambda_1*((f1)^0.5)+l*gamma0[p]
      
      A_2=-Sg[p]-n*lambda_1*((1-f1)^0.5)+l*beta0[p]
      A_3=-Sb[p]+n*lambda_1*((f1)^0.5)+l*gamma0[p]
      A_4=-Sg[p]+n*lambda_1*((1-f1)^0.5)+l*beta0[p]
      
      if (A_1/(l+sum(X_sb[,p]^2))>0){  
        
        beta[p]=A_1/(l+sum(X_sb[,p]^2))
      }
      if (A_3/(l+sum(X_sb[,p]^2))<0){  
        
        beta[p]=A_3/(l+sum(X_sb[,p]^2))
      }
      if (A_2/(sum(X_sg[,p]^2)+l)>0){  
        
        gamma[p]=A_2/(sum(X_sg[,p]^2)+l)
      } 
      if (A_4/(sum(X_sg[,p]^2)+l)<0){  
        
        gamma[p]=A_4/(sum(X_sg[,p]^2)+l)
      }
      
      diff=mean(abs(gamma-gamma0)+abs(beta-beta0)+abs(r10-r11)+abs(r20-r21))
    }
  }
  b1=beta/r11
  b1[which(abs(b1)<10^-3)]=0
  b2=gamma/r21
  b2[which(abs(b2)<10^-3)]=0
  sigma1=1/r11
  sigma2=1/r21
  
  list(b1,b2,f1,sigma1,sigma2,w1) 
  
}#end
FMR_fused<- function(X,Y,lambda_1,lambda_2){#############dim(X)=n*p,
  n          <- nrow(X)
  p          <- ncol(X)
  
  b=rep(0.1,p)
  beta=as.matrix(b,nc=1)  ##initial guess
  beta0=as.matrix(b/2,nc=1) 
  gamma=as.matrix(b/3,nc=1)
  gamma0=as.matrix(b/2,nc=1)
  f1=0.3
  f10=0.4
  w1=w2=rep(0.3,n)
  w10=w20=rep(0.4,n)
  r10=2
  r11=1
  r20=2
  r21=1
  tau=2*10^-1
  m=1
  while(abs(sum(w1-w10))>0.01){######E step
    m=m+1
    if (m>10){
      break
    }
    
    w10=w1
    w20=w2
    p12<- (f1/(1-f1))*(r11/r21)*exp(-1/2*(r11*Y-X%*%beta)^2+1/2*(r21*Y-X%*%gamma)^2)
    w1<- 1/(1+1/p12)###n*1   conditional probabilities that individual belongs to the cured or uncured 
    w2<-1/(p12+1)
    f1<- mean(w1)
    diff=1
    u=1
    while(diff>0.1){
      u=u+1
      if (u>200){
        break
      }
      
      r10=r11
      r20=r21
      
      f1<- mean(w1)
      Sb=c()
      Sg=c()
      #############M step#########
      beta0=beta
      gamma0=gamma
      Y_sb=as.vector(sqrt(w1))*Y
      X_sb=as.vector(sqrt(w1))*X
      #################################### for j=1
      r11=(sum(Y_sb*as.vector(X_sb%*%beta))+sqrt(sum(Y_sb*as.vector(X_sb%*%beta))^2+4*sum(Y_sb^2)*sum(w1)))/(2*sum(Y_sb^2))
      sum_sj20=0
      for (k in 2:p){
        sum_sj21=sum_sj20+beta0[k]*sum(X_sb[,1]*X_sb[,k])
        sum_sj20=sum_sj21
      }
      Sb[1]=-r11*sum(Y_sb*X_sb[,1])+sum_sj21
      
      Y_sg=as.vector(sqrt(w2))*Y
      X_sg=as.vector(sqrt(w2))*X
      r21=(sum(Y_sg*as.vector(X_sg%*%gamma))+sqrt(sum(Y_sg*as.vector(X_sg%*%gamma))^2+4*sum(Y_sg^2)*sum(w2)))/(2*sum(Y_sg^2))
      
      sum_sj20=0
      for (k in 2:p){
        sum_sj21=sum_sj20+gamma0[k]*sum(X_sg[,1]*X_sg[,k])
        sum_sj20=sum_sj21
      }
      Sg[1]=-r21*sum(Y_sg*X_sg[,1])+sum_sj21
      beta[1]=0
      gamma[1]=0
      ######beta>0,gamma>0, beta>gamma
      index=1
      
      A_1=-Sb[index]-n*lambda_1*((f1)^0.5)   ########mlasso  >0
      
      A_2=-Sg[index]-n*lambda_1*((1-f1)^0.5)  ########mlasso  >0
      A_3=-Sb[index]+n*lambda_1*((f1)^0.5)   ########mlasso  <0
      A_4=-Sg[index]+n*lambda_1*((1-f1)^0.5)  ########mlasso  <0
      if (((A_1)/(sum(X_sb[,index]^2))==(A_2)/(sum(X_sg[,index]^2)))&((A_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3)/(sum(X_sb[,index]^2))==(A_4)/(sum(X_sg[,index]^2)))&((A_4)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4)/(sum(X_sg[,index]^2))
        
      }
      #####################################
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=0
        
      }
      if (((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=0
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=0
        
      }
      if (((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=0
        gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      
      ####################
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>(A_2+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_2+n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_1+n*lambda_2)/(sum(X_sb[,index]^2))<(A_2-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_1+n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<(A_4-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_4-n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3-n*lambda_2)/(sum(X_sb[,index]^2))>(A_4+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_3-n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)&((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)&((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      } 
      ###########################################################  
      #################################### for j=2:(p-1)
      for (j in 2:(p-1)){
        sum_sj10=0
        for (k in 1:(j-1)){
          sum_sj11=sum_sj10+beta[k]*sum(X_sb[,j]*X_sb[,k])
          sum_sj10=sum_sj11
        }
        sum_sj20=0
        for(k in c((j+1):p)){
          sum_sj21=sum_sj20+beta0[k]*sum(X_sb[,j]*X_sb[,k])
          sum_sj20=sum_sj21
        }
        Sb[j]=-r11*sum(Y_sb*X_sb[,j])+sum_sj11+sum_sj21
        
        sum_sj10=0
        for (k in 1:(j-1)){
          sum_sj11=sum_sj10+gamma[k]*sum(X_sg[,j]*X_sg[,k])
          sum_sj10=sum_sj11
        }
        sum_sj20=0
        for(k in c((j+1):p)){
          sum_sj21=sum_sj20+gamma0[k]*sum(X_sg[,j]*X_sg[,k])
          sum_sj20=sum_sj21
        }
        Sg[j]=-r21*sum(Y_sg*X_sg[,j])+sum_sj11+sum_sj21
        beta[j]=0
        gamma[j]=0
        ##########################
        index=j
        
        A_1=-Sb[index]-n*lambda_1*((f1)^0.5)   ########mlasso  >0
        
        A_2=-Sg[index]-n*lambda_1*((1-f1)^0.5)  ########mlasso  >0
        A_3=-Sb[index]+n*lambda_1*((f1)^0.5)   ########mlasso  <0
        A_4=-Sg[index]+n*lambda_1*((1-f1)^0.5)  ########mlasso  <0
        
        
        
        if (((A_1)/(sum(X_sb[,index]^2))==(A_2)/(sum(X_sg[,index]^2)))&((A_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
          
          beta[index]=(A_1)/(sum(X_sb[,index]^2))
          gamma[index]=(A_2)/(sum(X_sg[,index]^2))
          
        }
        if (((A_3)/(sum(X_sb[,index]^2))==(A_4)/(sum(X_sg[,index]^2)))&((A_4)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
          
          beta[index]=(A_3)/(sum(X_sb[,index]^2))
          gamma[index]=(A_4)/(sum(X_sg[,index]^2))
          
        }
        #####################################
        if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
          
          beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=0
          
        }
        if (((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
          
          beta[index]=0
          gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        
        if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
          
          beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=0
          
        }
        if (((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
          
          beta[index]=0
          gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        
        ####################
        if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>(A_2+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_2+n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
          
          beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        if (((A_1+n*lambda_2)/(sum(X_sb[,index]^2))<(A_2-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_1+n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
          
          beta[index]=(A_1+n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        
        
        if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<(A_4-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_4-n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
          
          beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_4-n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        
        
        if (((A_3-n*lambda_2)/(sum(X_sb[,index]^2))>(A_4+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_3-n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
          
          beta[index]=(A_3-n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        
        if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)&((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
          
          beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
          
        }
        if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)&((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
          
          beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
          gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
          
        } 
      }
      
      
      #################################### for j=p
      
      sum_sj10=0
      for (k in 1:(p-1)){
        sum_sj11=sum_sj10+beta[k]*sum(X_sb[,p]*X_sb[,k])
        sum_sj10=sum_sj11
      }
      
      Sb[p]=-r11*sum(Y_sb*X_sb[,p])+sum_sj11
      
      sum_sj10=0
      for (k in 1:(p-1)){
        sum_sj11=sum_sj10+gamma[k]*sum(X_sg[,p]*X_sg[,k])
        sum_sj10=sum_sj11
      }
      
      Sg[p]=-r21*sum(Y_sg*X_sg[,p])+sum_sj11
      
      beta[p]=0
      gamma[p]=0
      ######
      index=p
      
      A_1=-Sb[index]-n*lambda_1*((f1)^0.5)   ########mlasso  >0
      
      A_2=-Sg[index]-n*lambda_1*((1-f1)^0.5)  ########mlasso  >0
      A_3=-Sb[index]+n*lambda_1*((f1)^0.5)   ########mlasso  <0
      A_4=-Sg[index]+n*lambda_1*((1-f1)^0.5)  ########mlasso  <0
      if (((A_1)/(sum(X_sb[,index]^2))==(A_2)/(sum(X_sg[,index]^2)))&((A_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3)/(sum(X_sb[,index]^2))==(A_4)/(sum(X_sg[,index]^2)))&((A_4)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4)/(sum(X_sg[,index]^2))
        
      }
      #####################################
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=0
        
      }
      if (((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=0
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=0
        
      }
      if (((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=0
        gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      
      ####################
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>(A_2+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_2+n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_1+n*lambda_2)/(sum(X_sb[,index]^2))<(A_2-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_1+n*lambda_2)/(sum(X_sb[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_1+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<(A_4-n*lambda_2)/(sum(X_sg[,index]^2)))&((A_4-n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3-n*lambda_2)/(sum(X_sb[,index]^2))>(A_4+n*lambda_2)/(sum(X_sg[,index]^2)))&((A_3-n*lambda_2)/(sum(X_sb[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_3-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_1-n*lambda_2)/(sum(X_sb[,index]^2))>0)&((A_4+n*lambda_2)/(sum(X_sg[,index]^2))<0)=="TRUE"){  
        
        beta[index]=(A_1-n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_4+n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      if (((A_3+n*lambda_2)/(sum(X_sb[,index]^2))<0)&((A_2-n*lambda_2)/(sum(X_sg[,index]^2))>0)=="TRUE"){  
        
        beta[index]=(A_3+n*lambda_2)/(sum(X_sb[,index]^2))
        gamma[index]=(A_2-n*lambda_2)/(sum(X_sg[,index]^2))
        
      }
      
      diff=mean(abs(gamma-gamma0)+abs(beta-beta0)+abs(r10-r11)+abs(r20-r21))
    }
  }
  b1=beta/r11
  b1[which(abs(b1)<10^-3)]=0
  b2=gamma/r21
  b2[which(abs(b2)<10^-3)]=0
  sigma1=1/r11
  sigma2=1/r21
  
  list(b1,b2,f1,sigma1,sigma2,w1) 
  
}#end

#############data generating
library(MASS)
set.seed(1000)
n                <- 200
pt=c(15,55,75,105,300,500)
########### Set the loop and p
loop=1
p                <- pt[floor((loop-1)/30)+1]
##########
f1_real=0.2
f2_real=0.8
############## the construction coeffient of sigma
rho        <- 0  
n_in_group <- 5
n_group    <- p%/%5
sigma      <- matrix(0,nr=p,nc=p) ###the structure between the covariables
sigma1     <- matrix(0,nr=n_in_group,nc=n_in_group)
dia      <- rep(1,p)

for(i in 1:n_in_group)
  for(j in i:n_in_group){
    sigma1[i,j]  <- dia[i]*rho^(j-i)
    sigma1[j,i]  <- sigma1[i,j]
  }
for(i in 1:n_group){
  sigma[(n_in_group*(i-1)+1):(n_in_group*i),(n_in_group*(i-1)+1):(n_in_group*i)]<-sigma1
}
X          <- mvrnorm(n, rep(0,p), sigma, empirical = F)


## regression 

beta_real <- c(0.05,0.2,0.05,0.2,3,3,0,0,rep(0,p-8))

gamma_real  <- c(0.05,0.2,0.05,0.2,0,0,-2,-2,rep(0,p-8))

###random

epi1<-0.5
epi2<- 0.5

Y=c()
g0=c(rep(1,round(n*f1_real)),rep(0,n-round(n*f1_real)))
for (i in 1:n){
  
  ###
  if (g0[i]==1){
    Y[i]=X[i,]%*%matrix(beta_real,nc=1)+as.matrix(rnorm(1,0,epi1),nc=1)}
  else{
    Y[i]=X[i,]%*%matrix(gamma_real,nc=1)+as.matrix(rnorm(1,0,epi2),nc=1)
  }
}
vary=f1_real*(t(as.matrix(beta_real,nc=1))%*%sigma%*%as.matrix(beta_real,nc=1)+epi1^2)+f2_real*(t(as.matrix(gamma_real,nc=1))%*%sigma%*%as.matrix(gamma_real,nc=1)+epi2^2)
SNR=vary/(f1_real*epi1^2+f2_real*epi2^2)

######################################################ROC curves
#########################hard threshholding
hard=seq(0,1,length=20)
V_h=matrix(0,2*p,length(hard))
r0=FMR_In(X,Y,0,0)
beta=r0[[1]]
gamma=r0[[2]]

for(ii in 1:length(hard)){
  
  hardthresh<- function (x) {
    if(abs(x)< hard[ii])
    {r=0}else{r=x}
    return(r)
  }
  
  beta_h=matrix(sapply(beta,hardthresh),nc=1)
  gamma_h=matrix(sapply(gamma,hardthresh),nc=1)
  
  theta_h=rbind(beta_h,gamma_h)
  V_h[,ii]=theta_h
}
#########################KM
Lambda_km           <- seq(0,0.2,length=20)
V_km=matrix(0,2*p,length(Lambda_km))
kl=kmeans(Y,2)
Cl_1=which(kl$cluster==1)
Cl_2=which(kl$cluster==2)
f1_km=length(Cl_1)/n
X_1=X[Cl_1,]
X_2=X[Cl_2,]

Y_1=Y[Cl_1]
Y_2=Y[Cl_2]

library("glmnet")
fit1<- glmnet(X_1,Y_1,family = "gaussian")
fit2<- glmnet(X_2,Y_2,family = "gaussian")

for (ii in 1: length(Lambda_km)){
  beta_km <- coef(fit1,s=Lambda_km[ii])[-1]
  gamma_km <- coef(fit2, s=Lambda_km[ii])[-1]
  theta_km=rbind(matrix(beta_km,nc=1),matrix(gamma_km,nc=1))
  V_km[,ii]=theta_km
}
#########################Lasso

Lambda1           <- seq(0,0.01,length=10)
V_l=matrix(0,2*p,length(Lambda1))
for(ii in 1:length(Lambda1)){
  r_l=FMR_In(X ,Y ,Lambda1[ii],0)
  
  beta_l=r_l[[1]]
  gamma_l=r_l[[2]]
  
  theta_l=rbind(beta_l,gamma_l)
  V_l[,ii]=theta_l
}


##############################  indicator
Lambda2           <- seq(0,0.05,0.005)
tuning=expand.grid(Lambda1 ,Lambda2)

V_in=matrix(0,2*p,length(tuning[[1]]))



for(ii in 1:length(tuning[[1]])){
  
  r_in=FMR_In(X,Y,tuning[[1]][ii],tuning[[2]][ii])
  
  beta_in=r_in[[1]]
  gamma_in=r_in[[2]]
  
  theta_in=rbind(beta_in,gamma_in)
  V_in[,ii]=theta_in
}

V_true=c(beta_real,gamma_real)

####################################fused
Lambda2           <- seq(0,0.01,0.001)
tuning=expand.grid(Lambda1 ,Lambda2)

V_fused=matrix(0,2*p,length(tuning[[1]]))

for(ii in 1:length(tuning[[1]])){
  
  r_fused=FMR_fused(X,Y,tuning[[1]][ii],tuning[[2]][ii])
  
  beta_fused=r_fused[[1]]
  gamma_fused=r_fused[[2]]
  
  theta_fused=rbind(beta_fused,gamma_fused)
  V_fused[,ii]=theta_fused
}

V=data.frame(V_true,V_h,V_km,V_l,V_in,V_fused)
write.csv(V,paste0('mixture',"corr0ROC",loop,'.csv'))


######################Estimation 
#################################### MLE
r0=FMR_In(X,Y,0,0)
beta=r0[[1]]
gamma=r0[[2]]

###########################hard threshholding
bic_h=c()
for(ii in 1:length(hard)){
    r=FMR_In(X,Y,0,0)
    beta_0=r[[1]]
    gamma_0=r[[2]]
    hardthresh<- function (x) {
      if(abs(x)< hard[ii])
      {r=0}else{r=x}
      return(r)
    }
    beta_h=sapply(beta_0,hardthresh)
    gamma_h=sapply(gamma_0,hardthresh)
    r11_h=r[[4]]
    r21_h=r[[5]]
    f1_h=r[[3]]
    de=2+1+length(which(beta_h!=0))+length(which(gamma_h!=0))
    bic_h[ii]=-2*sum(log(f1_h*dnorm(Y,X%*%beta_h,r11_h)+(1-f1_h)*dnorm(Y,X%*%gamma_h,r21_h)))+log(n)*de
}

i_th            <- which.min(bic_h)
Lambda_hard    <-hard[i_th]
hardthresh<- function (x) {
  if(abs(x)< Lambda_hard)
  {r=0}else{r=x}
  return(r)
}

beta_h=sapply(beta,hardthresh)
gamma_h=sapply(gamma,hardthresh)

#######################kmeans+Lasso

kl=kmeans(Y,2)
Cl_1=which(kl$cluster==1)
Cl_2=which(kl$cluster==2)
f1_km=length(Cl_1)/n
X_1=X[Cl_1,]
X_2=X[Cl_2,]

Y_1=Y[Cl_1]
Y_2=Y[Cl_2]

fit1<- glmnet(X_1,Y_1,family = "gaussian",intercept=FALSE)
fit2<- glmnet(X_2,Y_2,family = "gaussian",intercept=FALSE)
cv.fit1=cv.glmnet(X_1,Y_1,family='gaussian',type.measure="deviance")
cv.fit2=cv.glmnet(X_2,Y_2,family='gaussian',type.measure="deviance")

beta_km=coef(fit1,s=cv.fit1$lambda.min)
gamma_km=coef(fit2,s=cv.fit2$lambda.min)

########################################FMRLasso
bic_l=c()
for(ii in 1:length(Lambda1)){

    r_l=FMR_In(X,Y,Lambda1[ii],0)
    beta_l=r_l[[1]]
    gamma_l=r_l[[2]]
    r11_l=r_l[[4]]
    r21_l=r_l[[5]]
    f1_l=r_l[[3]]
    de=2+1+length(which(beta_l!=0))+length(which(gamma_l!=0))

    bic_l[ii]=-2*sum(log(f1_l*dnorm(Y,X%*%beta_l,r11_l)+(1-f1_l)*dnorm(Y,X%*%gamma_l,r21_l)))+log(n)*de

}

i_th            <- which.min(bic_l)
Lambda_lasso    <-Lambda1[i_th]
r_l=FMR_In(X,Y,Lambda_lasso,0)

######################################### FMRLASSO+Indicator
bic_in=c()
tuning=expand.grid(Lambda1 ,Lambda2)

for(ii in 1:length(tuning[[1]])){
  r_in=FMR_In(X,Y,tuning[[1]][ii],tuning[[2]][ii])
  beta_in=r_in[[1]]
  gamma_in=r_in[[2]]
  f1_in=r_in[[3]]
  r11_in=r_in[[4]]
  r21_in=r_in[[5]]
  de=2+1+length(which(beta_in!=0))+length(which(gamma_in!=0))

  bic_in[ii]=-2*sum(log(f1_in*dnorm(Y,X%*%beta_in,r11_in)+(1-f1_in)*dnorm(Y,X%*%gamma_in,r21_in)))+log(n)*de


  }

i_th            <- which.min(bic_in)
L1=tuning[[1]][i_th]
L2=tuning[[2]][i_th]
r_in=FMR_In(X,Y,L1,L2)

########################################FMRfused
bic_fused=c()
tuning=expand.grid(Lambda1 ,Lambda2)

for(ii in 1:length(tuning[[1]])){
  r_fused=FMR_fused(X,Y,tuning[[1]][ii],tuning[[2]][ii])
  beta_fused=r_fused[[1]]
  gamma_fused=r_fused[[2]]
  f1_fused=r_fused[[3]]
  r11_fused=r_fused[[4]]
  r21_fused=r_fused[[5]]
  de=2+1+length(which(beta_fused!=0))+length(which(gamma_fused!=0))

  bic_fused[ii]=-2*sum(log(f1_fused*dnorm(Y,X%*%beta_fused,r11_fused)+(1-f1_fused)*dnorm(Y,X%*%gamma_fused,r21_fused)))+log(n)*de


}

i_th            <- which.min(bic_fused)
L1=tuning[[1]][i_th]
L2=tuning[[2]][i_th]
r_fused=FMR_fused(X,Y,L1,L2)

############################coefficients
beta=r0[[1]]
beta_l=r_l[[1]]
beta_in=r_in[[1]]
beta_fused=r_fused[[1]]

gamma=r0[[2]]
gamma_l=r_l[[2]]
gamma_in=r_in[[2]]
gamma_fused=r_fused[[2]]


##TP
######mle
pred_tp=c(beta_real!=0 & beta!=0)
TP_beta=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma!=0)
TP_gamma=length(pred_tp[which(pred_tp==TRUE)])
TP=TP_beta+TP_gamma
TPR=TP/(sum(beta_real!=0)+sum(gamma_real!=0))
###hard threshholding
pred_tp=c(beta_real!=0 & beta_h!=0)
TP_beta_h=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma_h!=0)
TP_gamma_h=length(pred_tp[which(pred_tp==TRUE)])
TP_h=TP_beta_h+TP_gamma_h
TPR_h=TP_h/(sum(beta_real!=0)+sum(gamma_real!=0))

###kmeans
pred_tp=c(beta_real!=0 & beta_km!=0)
TP_beta_km=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma_km!=0)
TP_gamma_km=length(pred_tp[which(pred_tp==TRUE)])
TP_km=TP_beta_km+TP_gamma_km
TPR_km=TP_km/(sum(beta_real!=0)+sum(gamma_real!=0))

###lasso
pred_tp=c(beta_real!=0 & beta_l!=0)
TP_beta_l=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma_l!=0)
TP_gamma_l=length(pred_tp[which(pred_tp==TRUE)])
TP_l=TP_beta_l+TP_gamma_l
TPR_l=TP_l/(sum(beta_real!=0)+sum(gamma_real!=0))

###indicator
pred_tp=c(beta_real!=0 & beta_in!=0)
TP_beta_in=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma_in!=0)
TP_gamma_in=length(pred_tp[which(pred_tp==TRUE)])
TP_in=TP_beta_in+TP_gamma_in
TPR_in=TP_in/(sum(beta_real!=0)+sum(gamma_real!=0))
###fused
pred_tp=c(beta_real!=0 & beta_fused!=0)
TP_beta_fused=length(pred_tp[which(pred_tp==TRUE)])
pred_tp=c(gamma_real!=0 & gamma_fused!=0)
TP_gamma_fused=length(pred_tp[which(pred_tp==TRUE)])
TP_fused=TP_beta_fused+TP_gamma_fused
TPR_fused=TP_fused/(sum(beta_real!=0)+sum(gamma_real!=0))
##FP
######mle
pred_fp=c(beta_real==0 & beta!=0)
FP_beta=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma!=0)
FP_gamma=length(pred_fp[which(pred_fp==TRUE)])
FP=FP_beta+FP_gamma
FDR=FP/(sum(beta_real==0)+sum(gamma_real==0))

###hard threshholding
pred_fp=c(beta_real==0 & beta_h!=0)
FP_beta_h=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma_h!=0)
FP_gamma_h=length(pred_fp[which(pred_fp==TRUE)])
FP_h=FP_beta_h+FP_gamma_h
FDR_h=FP_h/(sum(beta_real==0)+sum(gamma_real==0))

###Kmeans
pred_fp=c(beta_real==0 & beta_km!=0)
FP_beta_km=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma_km!=0)
FP_gamma_km=length(pred_fp[which(pred_fp==TRUE)])
FP_km=FP_beta_km+FP_gamma_km
FDR_km=FP_km/(sum(beta_real==0)+sum(gamma_real==0))


###lasso
pred_fp=c(beta_real==0 & beta_l!=0)
FP_beta_l=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma_l!=0)
FP_gamma_l=length(pred_fp[which(pred_fp==TRUE)])
FP_l=FP_beta_l+FP_gamma_l
FDR_l=FP_l/(sum(beta_real==0)+sum(gamma_real==0))

###indicator
pred_fp=c(beta_real==0 & beta_in!=0)
FP_beta_in=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma_in!=0)
FP_gamma_in=length(pred_fp[which(pred_fp==TRUE)])
FP_in=FP_beta_in+FP_gamma_in
FDR_in=FP_in/(sum(beta_real==0)+sum(gamma_real==0))

##Fused
pred_fp=c(beta_real==0 & beta_fused!=0)
FP_beta_fused=length(pred_fp[which(pred_fp==TRUE)])
pred_fp=c(gamma_real==0 & gamma_fused!=0)
FP_gamma_fused=length(pred_fp[which(pred_fp==TRUE)])
FP_fused=FP_beta_fused+FP_gamma_fused
FDR_fused=FP_fused/(sum(beta_real==0)+sum(gamma_real==0))


re=data.frame(TPR,TPR_h,TPR_km,TPR_l,TPR_in,TPR_fused,FDR,FDR_h,FDR_km,FDR_l,FDR_in,FDR_fused,beta,beta_h,beta_km,beta_l,beta_in,beta_fused,gamma,gamma_h,gamma_km,gamma_l,gamma_in,gamma_fused)

write.csv(re,paste0('mixture',"corr0",loop,'.csv'))
