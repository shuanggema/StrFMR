######################################Reading cutaneous melanoma (SKCM) data

X_mela=read.csv("Z:\\User\\Documents\\MQ\\Mixture result\\real data analysis\\SelectX_mela.csv")
Y_mela=read.csv("Z:\\User\\Documents\\MQ\\Mixture result\\real data analysis\\Y_mela.csv")

mixture=union(which(Y_mela$AJCC_PATHOLOGIC_TUMOR_STAGE==2),which(Y_mela$AJCC_PATHOLOGIC_TUMOR_STAGE==3))
Y=matrix(Y_mela[mixture,2],nc=1)
X=as.matrix(X_mela[mixture,],nr=length(mixture))
X=t(na.omit(t(apply(X,2,scale))))
Y=apply(Y,2,scale)

mixture_1=which(Y_mela$AJCC_PATHOLOGIC_TUMOR_STAGE==2)
mixture_2=which(Y_mela$AJCC_PATHOLOGIC_TUMOR_STAGE==3)

X1=X_mela[mixture_1,]
X2=X_mela[mixture_2,]

Y1=Y_mela[mixture_1,2]
Y2=Y_mela[mixture_2,2]

###################################screening

screen=function(x){
  #r=cor(x,Y)
  r=cor.test(x,Y)
  r$p.value
}

p=apply(X, 2,screen)
se=c()
for (i in 1:300){
  se[i]=which(names(p)==names(sort(p,decreasing = FALSE)[1:300])[i])
}

X_s=X[,se]
gene=names(p[se])
##########################
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
  
}#e
#############data generating

p=dim(X_s)[2]
n=dim(X_s)[1]
#######################
# Lambda1           <- seq(0,0.3,length=5)
# Lambda2           <- seq(0,0.01,0.002)
Lambda1           <- seq(0.001,0.05,length=5)
Lambda2           <- seq(0,0.0008,length=9)

tuning=expand.grid(Lambda1 ,Lambda2)

###################MLE
r0=FMR_In(X_s,Y,0,0)
###################################################LASSO
bic_l=c()
for(ii in 1:length(Lambda1)){
 
    r_l=FMR_In(X_s,Y,Lambda1[ii],0)
    beta_l=r_l[[1]]
    gamma_l=r_l[[2]]
    r11_l=r_l[[4]]
    r21_l=r_l[[5]]
    f1_l=r_l[[3]]
    de=2+1+length(which(beta_l!=0))+length(which(gamma_l!=0))
    bic_l[ii]=-2*sum(log(f1_l*dnorm(Y,X_s%*%beta_l,r11_l)+(1-f1_l)*dnorm(Y,X_s%*%gamma_l,r21_l)))+log(n)*de

  }

i_th            <- which.min(bic_l)
i_th=1
Lambda_lasso    <-Lambda1[i_th]
r_l=FMR_In(X_s,Y,Lambda_lasso,0)
beta_l=r_l[[1]]
gamma_l=r_l[[2]]
f1_l=r_l[[3]]
#####################################################proposed
bic_in=c()
for(ii in 1:length(tuning[[1]])){
 
    r_in=FMR_In(X_s,Y,tuning[[1]][ii],tuning[[2]][ii])
    beta_in=r_in[[1]]
    gamma_in=r_in[[2]]

    f1_in=r_in[[3]]
    r11_in=r_in[[4]]
    r21_in=r_in[[5]]
    de=2+1+length(which(beta_in!=0))+length(which(gamma_in!=0))
    bic_in[ii]=-2*sum(log(f1_in*dnorm(Y,X_s%*%beta_in,r11_in)+(1-f1_in)*dnorm(Y,X_s%*%gamma_in,r21_in)))+log(n)*de
    
}

i_th            <- which.min(bic_in)
L1=tuning[[1]][i_th]
L2=tuning[[2]][i_th]
r_in=FMR_In(X_s,Y,L1,L2)
beta_in=r_in[[1]]
gamma_in=r_in[[2]]


beta=r0[[1]]
gamma=r0[[2]]

f1=r0[[3]]
f1_l=r_l[[3]]
f1_in=r_in[[3]]
r11=r0[[4]]
r11_l=r_l[[4]]
r11_in=r_in[[4]]
r21=r0[[5]]
r21_l=r_l[[5]]
r21_in=r_in[[5]]


re=data.frame(colnames(X_s),Lambda_lasso,L1,L2,beta,beta_l,beta_in,gamma,gamma_l,gamma_in,f1,f1_l,f1_in,r11,r11_l,r11_in,r21,r21_l,r21_in)

write.csv(re,paste0('realdata',"ME",'.csv'))


