#install.packages("glmnet");install.packages("itertools");install.packages("foreach");install.packages("doParallel");install.packages("RCAL");install.packages("truncnorm");install.packages("ranger")
require(glmnet);require(itertools);require(foreach);require(doParallel);require(truncnorm);require(RCAL);require(ranger)
ptm=proc.time()
n_cores=32
registerDoParallel(n_cores)
setting="a";N=400;d1=100;d2=50;d=d1+d2;tmax=200;K=5;value_trunc=1;M=floor(N/K);num_est=9;num_totest=num_est+1;
logistic=function(x){exp(x)/(1+exp(x))}
func_link=function(x){(abs(x+1)+0.1)/(abs(x+1)+1)}
# Parameters
beta_pi1_0=0
beta_pi1_S1=c(rep(1,2),rep(0,d1-2))
beta_pi2_0_0=0
beta_pi2_0_S1=c(1,0,0,rep(0,d1-3))
beta_pi2_0_S2=c(rep(0.5,4),rep(0,d2-4))
beta_pi2_1_0=0
beta_pi2_1_S1=-c(1,0,0,rep(0,d1-3))
beta_pi2_1_S2=-c(rep(0.5,4),rep(0,d2-4))
beta_g0_0=1
beta_g0_S1=c(1,rep(0,d1-1))
beta_g0_S2=c(rep(0.5,4),rep(0,d2-4))
beta_g1_0=-1
beta_g1_S1=-c(1,rep(0,d1-1))
beta_g1_S2=-c(rep(0.5,4),rep(0,d2-4))
W1=matrix(0,nrow=d1,ncol=d2)
W0=matrix(0,nrow=d1,ncol=d2)
for (i in 1:d1){
  for (j in 1:d2){
    if (i==j){
      W1[i,j]=1
      W0[i,j]=1
    } else if (abs(i-j)==1){
      W1[i,j]=0.8
      W0[i,j]=0.7
    } else if (abs(i-j)==2){
      W0[i,j]=0.7^2
    }
  }
}
Q1=W1*0.5;Q0=W0*0.5
gamma_1_0=beta_g1_0+sum(beta_g1_S2)
gamma_0_0=beta_g0_0
theta=gamma_1_0-gamma_0_0

set.seed(1234)
seeds <- sample.int(1e8, tmax) # Generate random seeds

results_par <- foreach (t=1:tmax, .errorhandling = "pass", .combine=rbind, .packages=c("glmnet","RCAL","ranger","truncnorm")) %dopar% {
  set.seed(seeds[t]) # Set seed inside the foreach loop to ensure reproducibility
  
  # Data generation
  delta1=rtruncnorm(N,-value_trunc, value_trunc)
  delta2=matrix(rtruncnorm(N*d2,-value_trunc, value_trunc),nrow=N)
  S1=matrix(rnorm(N*d1),nrow=N)
  pi1=logistic(beta_pi1_0+S1%*%beta_pi1_S1)
  A1=rbinom(N,1,pi1)
  S2=(S1^2-1)%*%Q1+S1%*%W1+1+delta1+delta2
  S2[A1==0]=(S1[A1==0,]^2-1)%*%Q0+S1[A1==0,]%*%W0+delta2[A1==0]
  S2bar=cbind(S1,S2)
  S2_1=(S1^2-1)%*%Q1+S1%*%W1+1+delta1+delta2
  S2_0=(S1^2-1)%*%Q0+S1%*%W0+delta2
  mu2_1=beta_g1_0+S1%*%beta_g1_S1+S2_1%*%beta_g1_S2
  mu2_0=beta_g0_0+S1%*%beta_g0_S1+S2_0%*%beta_g0_S2
  gamma_1_S1=beta_g1_S1+W1%*%beta_g1_S2
  gamma_0_S1=beta_g0_S1+W0%*%beta_g0_S2
  mu1_1=gamma_1_0+S1%*%gamma_1_S1+(S1^2-1)%*%Q1%*%beta_g1_S2
  mu1_0=gamma_0_0+S1%*%gamma_0_S1+(S1^2-1)%*%Q0%*%beta_g0_S2
  ww1=beta_pi2_0_0+S1%*%beta_pi2_0_S1+S2_0%*%beta_pi2_0_S2
  pi2_0=func_link(ww1)
  ww2=beta_pi2_1_0+S1%*%beta_pi2_1_S1+S2_1%*%beta_pi2_1_S2
  pi2_1=func_link(ww2)
  pi2=(1-A1)*pi2_0+A1*pi2_1
  A2=rbinom(N,1,pi2)
  g=(A1+A2==0)*(beta_g0_0+S1%*%beta_g0_S1+S2%*%beta_g0_S2)+(A1*A2==1)*(beta_g1_0+S1%*%beta_g1_S1+S2%*%beta_g1_S2)
  zeta=rnorm(N)
  Y=g+zeta
  
  # Initialization
  pred_pi1_1=matrix(0,nrow=num_est,ncol=N)
  pred_pi2_0=matrix(0,nrow=num_est,ncol=N);pred_pi2_1=matrix(0,nrow=num_est,ncol=N)
  pred_mu2_1=matrix(0,nrow=num_est,ncol=N);pred_mu1_1=matrix(0,nrow=num_est,ncol=N)
  pred_mu2_0=matrix(0,nrow=num_est,ncol=N);pred_mu1_0=matrix(0,nrow=num_est,ncol=N)
  
  # Oracle
  pred_pi1_1[num_est,]=pi1
  pred_pi2_1[num_est,]=pi2_1;pred_pi2_0[num_est,]=pi2_0
  pred_mu2_1[num_est,]=mu2_1;pred_mu2_0[num_est,]=mu2_0
  pred_mu1_1[num_est,]=mu1_1;pred_mu1_0[num_est,]=mu1_0
  
  time_all=rep(0,num_totest)
  nrho=100
  # SMDR1
  time=proc.time()[3]
  pred_pi1_1_0_SMDR1=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index_train1=index_nk1[1:floor(length(index_nk1)/2)]
    index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
    index_train3=index_nk2[1:floor(length(index_nk2)/2)]
    index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_train1,],y=A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha1=c(-log(1/mean(A1[index_train1]) - 1), rep(0,d1)) # Use a constant estimate if no convergence for all rho values
    pred_pi1_1[1,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_train1,],y=1-A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_train1]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR1[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A1=0|X), a seperate PS estimate to ensure orthogonality for theta0
    #pi2_1
    iw=as.numeric(A1[index_train2]/logistic(alpha1[1]+S1[index_train2,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_train2]) - 1), rep(0,d))
    pred_pi2_1[1,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_train2])/logistic(alpha0[1]+S1[index_train2,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=1-A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_train2])) - 1), rep(0,d))
    pred_pi2_0[1,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=A1[index_train3]*A2[index_train3]*exp(-beta1[1]-S2bar[index_train3,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_train3,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[1,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=(1-A1[index_train3])*(1-A2[index_train3])*exp(-beta0[1]-S2bar[index_train3,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_train3,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[1,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new1+A2[index_train4]*(Y[index_train4]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_train4,]%*%beta1[-1]), weights=A1[index_train4]*exp(-alpha1[1]-S1[index_train4,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[1,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new0+(1-A2[index_train4])*(Y[index_train4]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_train4,]%*%beta0[-1]), weights=(1-A1[index_train4])*exp(-alpha0[1]-S1[index_train4,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[1,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[1]=time_all[1]+proc.time()[3]-time
  
  # SMDR2
  time=proc.time()[3]
  pred_pi1_1_0_SMDR2=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_nk,],y=A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1[2,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_nk,],y=1-A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR2[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A0=1|X)
    #pi2_1
    iw=as.numeric(A1[index_nk]/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_nk]) - 1), rep(0,d))
    pred_pi2_1[2,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_nk])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=1-A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_nk])) - 1), rep(0,d))
    pred_pi2_0[2,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=A1[index_nk]*A2[index_nk]*exp(-beta1[1]-S2bar[index_nk,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[2,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=(1-A1[index_nk])*(1-A2[index_nk])*exp(-beta0[1]-S2bar[index_nk,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[2,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new1+A2[index_nk]*(Y[index_nk]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_nk,]%*%beta1[-1]), weights=A1[index_nk]*exp(-alpha1[1]-S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[2,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new0+(1-A2[index_nk])*(Y[index_nk]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_nk,]%*%beta0[-1]), weights=(1-A1[index_nk])*exp(-alpha0[1]-S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[2,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[2]=time_all[2]+proc.time()[3]-time
  
  nrho=50
  # SMDR1
  time=proc.time()[3]
  pred_pi1_1_0_SMDR3=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index_train1=index_nk1[1:floor(length(index_nk1)/2)]
    index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
    index_train3=index_nk2[1:floor(length(index_nk2)/2)]
    index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_train1,],y=A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha1=c(-log(1/mean(A1[index_train1]) - 1), rep(0,d1)) # Use a constant estimate if no convergence for all rho values
    pred_pi1_1[3,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_train1,],y=1-A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_train1]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR3[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A1=0|X), a seperate PS estimate to ensure orthogonality for theta0
    #pi2_1
    iw=as.numeric(A1[index_train2]/logistic(alpha1[1]+S1[index_train2,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_train2]) - 1), rep(0,d))
    pred_pi2_1[3,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_train2])/logistic(alpha0[1]+S1[index_train2,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=1-A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_train2])) - 1), rep(0,d))
    pred_pi2_0[3,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=A1[index_train3]*A2[index_train3]*exp(-beta1[1]-S2bar[index_train3,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_train3,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[3,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=(1-A1[index_train3])*(1-A2[index_train3])*exp(-beta0[1]-S2bar[index_train3,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_train3,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[3,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new1+A2[index_train4]*(Y[index_train4]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_train4,]%*%beta1[-1]), weights=A1[index_train4]*exp(-alpha1[1]-S1[index_train4,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[3,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new0+(1-A2[index_train4])*(Y[index_train4]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_train4,]%*%beta0[-1]), weights=(1-A1[index_train4])*exp(-alpha0[1]-S1[index_train4,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[3,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[3]=time_all[3]+proc.time()[3]-time
  
  # SMDR2
  time=proc.time()[3]
  pred_pi1_1_0_SMDR4=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_nk,],y=A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1[4,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_nk,],y=1-A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR4[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A0=1|X)
    #pi2_1
    iw=as.numeric(A1[index_nk]/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_nk]) - 1), rep(0,d))
    pred_pi2_1[4,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_nk])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=1-A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_nk])) - 1), rep(0,d))
    pred_pi2_0[4,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=A1[index_nk]*A2[index_nk]*exp(-beta1[1]-S2bar[index_nk,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[4,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=(1-A1[index_nk])*(1-A2[index_nk])*exp(-beta0[1]-S2bar[index_nk,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[4,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new1+A2[index_nk]*(Y[index_nk]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_nk,]%*%beta1[-1]), weights=A1[index_nk]*exp(-alpha1[1]-S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[4,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new0+(1-A2[index_nk])*(Y[index_nk]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_nk,]%*%beta0[-1]), weights=(1-A1[index_nk])*exp(-alpha0[1]-S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[4,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[4]=time_all[4]+proc.time()[3]-time
  
  nrho=10
  # SMDR1
  time=proc.time()[3]
  pred_pi1_1_0_SMDR5=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index_train1=index_nk1[1:floor(length(index_nk1)/2)]
    index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
    index_train3=index_nk2[1:floor(length(index_nk2)/2)]
    index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_train1,],y=A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha1=c(-log(1/mean(A1[index_train1]) - 1), rep(0,d1)) # Use a constant estimate if no convergence for all rho values
    pred_pi1_1[5,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_train1,],y=1-A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_train1]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR5[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A1=0|X), a seperate PS estimate to ensure orthogonality for theta0
    #pi2_1
    iw=as.numeric(A1[index_train2]/logistic(alpha1[1]+S1[index_train2,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_train2]) - 1), rep(0,d))
    pred_pi2_1[5,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_train2])/logistic(alpha0[1]+S1[index_train2,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=1-A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_train2])) - 1), rep(0,d))
    pred_pi2_0[5,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=A1[index_train3]*A2[index_train3]*exp(-beta1[1]-S2bar[index_train3,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_train3,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[5,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=(1-A1[index_train3])*(1-A2[index_train3])*exp(-beta0[1]-S2bar[index_train3,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_train3,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[5,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new1+A2[index_train4]*(Y[index_train4]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_train4,]%*%beta1[-1]), weights=A1[index_train4]*exp(-alpha1[1]-S1[index_train4,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[5,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new0+(1-A2[index_train4])*(Y[index_train4]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_train4,]%*%beta0[-1]), weights=(1-A1[index_train4])*exp(-alpha0[1]-S1[index_train4,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[5,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[5]=time_all[5]+proc.time()[3]-time
  
  # SMDR2
  time=proc.time()[3]
  pred_pi1_1_0_SMDR6=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_nk,],y=A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1[6,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_nk,],y=1-A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR6[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A0=1|X)
    #pi2_1
    iw=as.numeric(A1[index_nk]/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_nk]) - 1), rep(0,d))
    pred_pi2_1[6,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_nk])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=1-A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_nk])) - 1), rep(0,d))
    pred_pi2_0[6,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=A1[index_nk]*A2[index_nk]*exp(-beta1[1]-S2bar[index_nk,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[6,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=(1-A1[index_nk])*(1-A2[index_nk])*exp(-beta0[1]-S2bar[index_nk,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[6,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new1+A2[index_nk]*(Y[index_nk]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_nk,]%*%beta1[-1]), weights=A1[index_nk]*exp(-alpha1[1]-S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[6,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new0+(1-A2[index_nk])*(Y[index_nk]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_nk,]%*%beta0[-1]), weights=(1-A1[index_nk])*exp(-alpha0[1]-S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[6,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[6]=time_all[6]+proc.time()[3]-time
  
  nrho=200
  # SMDR1
  time=proc.time()[3]
  pred_pi1_1_0_SMDR7=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index_train1=index_nk1[1:floor(length(index_nk1)/2)]
    index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
    index_train3=index_nk2[1:floor(length(index_nk2)/2)]
    index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_train1,],y=A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha1=c(-log(1/mean(A1[index_train1]) - 1), rep(0,d1)) # Use a constant estimate if no convergence for all rho values
    pred_pi1_1[7,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_train1,],y=1-A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_train1]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR7[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A1=0|X), a seperate PS estimate to ensure orthogonality for theta0
    #pi2_1
    iw=as.numeric(A1[index_train2]/logistic(alpha1[1]+S1[index_train2,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_train2]) - 1), rep(0,d))
    pred_pi2_1[7,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_train2])/logistic(alpha0[1]+S1[index_train2,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_train2,],y=1-A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_train2])) - 1), rep(0,d))
    pred_pi2_0[7,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=A1[index_train3]*A2[index_train3]*exp(-beta1[1]-S2bar[index_train3,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_train3,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[7,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=(1-A1[index_train3])*(1-A2[index_train3])*exp(-beta0[1]-S2bar[index_train3,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_train3,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[7,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new1+A2[index_train4]*(Y[index_train4]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_train4,]%*%beta1[-1]), weights=A1[index_train4]*exp(-alpha1[1]-S1[index_train4,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[7,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_train4,],y=pred_mu2_new0+(1-A2[index_train4])*(Y[index_train4]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_train4,]%*%beta0[-1]), weights=(1-A1[index_train4])*exp(-alpha0[1]-S1[index_train4,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[7,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[7]=time_all[7]+proc.time()[3]-time
  
  # SMDR2
  time=proc.time()[3]
  pred_pi1_1_0_SMDR8=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    #pi1_1
    fit=glm.regu.cv(x=S1[index_nk,],y=A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha1=fit$sel.bet[,1]
    if (sum(is.na(alpha1))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1[8,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
    #pi1_0
    fit=glm.regu.cv(x=S1[index_nk,],y=1-A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    alpha0=fit$sel.bet[,1]
    if (sum(is.na(alpha0))>0) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
    pred_pi1_1_0_SMDR8[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A0=1|X)
    #pi2_1
    iw=as.numeric(A1[index_nk]/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta1=fit$sel.bet[,1]
    if (sum(is.na(beta1))>0) beta1=c(-log(mean(iw)/mean(iw*A2[index_nk]) - 1), rep(0,d))
    pred_pi2_1[8,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
    #pi2_0
    iw=as.numeric((1-A1[index_nk])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]))
    fit=glm.regu.cv(x=S2bar[index_nk,],y=1-A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1)))
    beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
    if (sum(is.na(beta0))>0) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_nk])) - 1), rep(0,d))
    pred_pi2_0[8,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
    #mu2_1
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=A1[index_nk]*A2[index_nk]*exp(-beta1[1]-S2bar[index_nk,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_1[8,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu2_0
    fit=cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=(1-A1[index_nk])*(1-A2[index_nk])*exp(-beta0[1]-S2bar[index_nk,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu2_0[8,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
    #mu1_1
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new1+A2[index_nk]*(Y[index_nk]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_nk,]%*%beta1[-1]), weights=A1[index_nk]*exp(-alpha1[1]-S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_1[8,index]=predict(fit,newx=S1[index,],s="lambda.min")
    #mu1_0
    fit=cv.glmnet(x=S1[index_nk,],y=pred_mu2_new0+(1-A2[index_nk])*(Y[index_nk]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_nk,]%*%beta0[-1]), weights=(1-A1[index_nk])*exp(-alpha0[1]-S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho)
    pred_mu1_0[8,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  time_all[8]=time_all[8]+proc.time()[3]-time
  
  
  # Estimation of theta
  pred_theta_t=c();pred_sig2_t=c()
  for (j in 1:num_est){
    gamma2_1=A1*A2/(pred_pi1_1[j,]*pred_pi2_1[j,]);gamma1_1=A1/pred_pi1_1[j,]
    psi_1=gamma2_1*Y-(gamma1_1-1)*pred_mu1_1[j,]-(gamma2_1-gamma1_1)*pred_mu2_1[j,]
    gamma2_0=(1-A1)*(1-A2)/((1-pred_pi1_1[j,])*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/(1-pred_pi1_1[j,])
    if (j==1){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR1*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR1
    }
    if (j==2){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR2*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR2
    }
    if (j==3){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR3*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR3
    }
    if (j==4){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR4*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR4
    }
    if (j==5){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR5*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR5
    }
    if (j==6){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR6*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR6
    }
    if (j==7){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR7*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR7
    }
    if (j==8){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR8*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR8
    }
    psi_0=gamma2_0*Y-(gamma1_0-1)*pred_mu1_0[j,]-(gamma2_0-gamma1_0)*pred_mu2_0[j,]
    pred_theta_t=c(pred_theta_t,mean(psi_1-psi_0))
  }
  
  # Asymptotic variance estimator
  for (j in 1:num_est){
    gamma2_1=A1*A2/(pred_pi1_1[j,]*pred_pi2_1[j,]);gamma1_1=A1/pred_pi1_1[j,]
    psi_1=gamma2_1*Y-(gamma1_1-1)*pred_mu1_1[j,]-(gamma2_1-gamma1_1)*pred_mu2_1[j,]
    gamma2_0=(1-A1)*(1-A2)/((1-pred_pi1_1[j,])*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/(1-pred_pi1_1[j,])
    if (j==1){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR1*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR1
    }
    if (j==2){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR2*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR2
    }
    if (j==3){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR3*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR3
    }
    if (j==4){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR4*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR4
    }
    if (j==5){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR5*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR5
    }
    if (j==6){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR6*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR6
    }
    if (j==7){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR7*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR7
    }
    if (j==8){
      gamma2_0=(1-A1)*(1-A2)/(pred_pi1_1_0_SMDR8*(1-pred_pi2_0[j,]));gamma1_0=(1-A1)/pred_pi1_1_0_SMDR8
    }
    psi_0=gamma2_0*Y-(gamma1_0-1)*pred_mu1_0[j,]-(gamma2_0-gamma1_0)*pred_mu2_0[j,]
    psi=psi_1-psi_0-pred_theta_t[j]
    pred_sig2_t=c(pred_sig2_t,mean(psi^2))
  }
  
  # Empdiff
  pred_theta_t=c(pred_theta_t,mean(Y[which(A1*A2==1)])-mean(Y[which(A1+A2==0)]))
  pred_sig2_t=c(pred_sig2_t,var(Y[which(A1*A2==1)])*sum(A1*A2==1)/N+var(Y[which(A1+A2==0)])*sum(A1+A2==0)/N)
  
  result=c(pred_theta=pred_theta_t,pred_sig2=pred_sig2_t)
  message(paste("Iteration ",t," done; Error:",paste(sprintf("%.3f", pred_theta_t-theta),collapse = " "),"; Time:",paste(sprintf("%.3f", time_all),collapse = " ")))
  result
}
pred_theta=results_par[,1:num_totest]
pred_sig2=results_par[,(num_totest+1):(2*num_totest)]

# Save the prediction values
filename1=paste("predtheta_DTE_tunning_a_N",N,"d1",d1,"d2",d2,".RData",sep="")
save(pred_theta,file=filename1)
filename2=paste("predvar_DTE_tunning_a_N",N,"d1",d1,"d2",d2,".RData",sep="")
save(pred_sig2,file=filename2)

# Set NaN prediction values to Inf
pred_theta[is.nan(pred_theta)] <- Inf
pred_sig2[is.nan(pred_sig2)] <- Inf

# Robust (median-type) empirical estimates
bias=apply(pred_theta-theta,2,median)
RMSE=apply((pred_theta-theta)^2,2,median)^0.5
AL=apply(2*sqrt(pred_sig2/N)*qnorm(1-0.05/2),2,median)
AC=apply(abs(pred_theta-theta)<=sqrt(pred_sig2/N)*qnorm(1-0.05/2),2,mean)
ESD=apply(pred_theta,2,function(x){1.4826*median(abs(x-median(x)))})
ASD=apply(pred_sig2^0.5,2,median)/sqrt(N)

# Report results in a Latex-friendly format
est=c("SMDR1-rho100","SMDR2-rho100","SMDR1-rho50","SMDR2-rho50","SMDR1-rho10","SMDR2-rho10","SMDR1-rho200","SMDR2-rho200","oracle","empdiff")
table_result=function(x){
  paste(est[x],"&",format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AL[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(AC[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(ESD[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(ASD[x],digits=3),nsmall=3,scientific=FALSE),sep="")
}
result=apply(t(1:num_totest),2,table_result)
result=list(result=c("estimator&Bias&RMSE&AL&AC&ESD&ASD",result))
result 

proc.time()-ptm


