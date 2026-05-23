require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL);require(ranger);require(readr);require(haven)
ptm=proc.time()
setting="b";n_cores=32;num_est=8;num_totest=num_est+1;tmax=200;K=5;nrho=50;level_signal=2
registerDoParallel(n_cores)

dat <- read_sas("bencost.sas7bdat")
data <- read_csv("data.csv")
tmp=cbind(dat[,34],data[,c(2:2340)])
tmp=tmp[complete.cases(tmp),]
tmp=as.data.frame(tmp)
data=tmp[tmp[,2] !=-1 & tmp[,3] !=-1,]
dat=0
idx=(data[,2] %in% c(1,3)) & (data[,3] %in% c(1,3)) # Consider population under treatments 1 or 3
data=data[idx,]

A1_all=data[,2];A2_all=data[,3];Y_all=data[,1]
indicator1 <- function(x) I(x == 3) + 1 - 1
N=length(A1_all)
A1_all=indicator1(A1_all)
A2_all=indicator1(A2_all)
Y_all=log(Y_all+1)
S_all=data[,5:ncol(data)]

# Remove the unbalance features
threshold=0.1;name_last_S1="FEMALE";idx_last_S1=which(names(S_all)==name_last_S1)
keep=c()
for (i in 1:ncol(S_all)){
  freq=table(factor(S_all[,i]))/nrow(S_all)
  keep=c(keep,max(freq)<=0.9)
}
keep_S1=keep[1:idx_last_S1]
S2bar_all=as.matrix(S_all[,keep])
d1=sum(keep_S1);d=sum(keep==1);d2=d-d1;M=floor(N/K)
N;d1;d

# Hyperparameter range for RF
mtry_range=floor(seq(from=1,to=50,length.out=5))
min.node.size_range=floor(seq(from=10,to=sqrt(N),length.out=3))
paras=numeric()
for (mtry in mtry_range){
  for (min.node.size in min.node.size_range){
    paras=rbind(paras,c(mtry,min.node.size))
  }
}

set.seed(1234)
logistic=function(x){exp(x)/(1+exp(x))}

resample <- sample.int(N)
S2bar=S2bar_all[resample,];S1=S2bar[,1:d1];S2=S2bar[,(d1+1):d]
A1=A1_all[resample];A2=A2_all[resample];Y=Y_all[resample]

# Initialization
pred_pi1_1=matrix(0,nrow=num_est,ncol=N)
pred_pi2_0=matrix(0,nrow=num_est,ncol=N);pred_pi2_1=matrix(0,nrow=num_est,ncol=N)
pred_mu2_1=matrix(0,nrow=num_est,ncol=N);pred_mu1_1=matrix(0,nrow=num_est,ncol=N)
pred_mu2_0=matrix(0,nrow=num_est,ncol=N);pred_mu1_0=matrix(0,nrow=num_est,ncol=N)

# SMDR1
pred_pi1_1_0_SMDR1=rep(0,N)
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  index_nk1=index_nk[1:floor(length(index_nk)/2)]
  index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  index_train1=index_nk1[1:floor(length(index_nk1)/2)]
  index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
  index_train3=index_nk2[1:floor(length(index_nk2)/2)]
  index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
  #pi1_1
  fit=tryCatch(glm.regu.cv(x=S1[index_train1,],y=A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  alpha1=fit$sel.bet[,1]
  if (sum(is.na(alpha1))>0|is.null(fit)) alpha1=c(-log(1/mean(A1[index_train1]) - 1), rep(0,d1)) # Use a constant estimate if no convergence for all rho values
  pred_pi1_1[1,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
  #pi1_0
  fit=tryCatch(glm.regu.cv(x=S1[index_train1,],y=1-A1[index_train1],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  alpha0=fit$sel.bet[,1]
  if (sum(is.na(alpha0))>0|is.null(fit)) alpha0=c(-log(1/mean(1-A1[index_train1]) - 1), rep(0,d1))
  pred_pi1_1_0_SMDR1[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A1=0|X), a seperate PS estimate to ensure orthogonality for theta0
  #pi2_1
  iw=as.numeric(A1[index_train2]/logistic(alpha1[1]+S1[index_train2,]%*%alpha1[-1]))
  fit=tryCatch(glm.regu.cv(x=S2bar[index_train2,],y=A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  beta1=fit$sel.bet[,1]
  if (sum(is.na(beta1))>0|is.null(fit)) beta1=c(-log(mean(iw)/mean(iw*A2[index_train2]) - 1), rep(0,d))
  pred_pi2_1[1,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
  #pi2_0
  iw=as.numeric((1-A1[index_train2])/logistic(alpha0[1]+S1[index_train2,]%*%alpha0[-1]))
  fit=tryCatch(glm.regu.cv(x=S2bar[index_train2,],y=1-A2[index_train2],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
  if (sum(is.na(beta0))>0|is.null(fit)) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_train2])) - 1), rep(0,d))
  pred_pi2_0[1,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
  #mu2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=A1[index_train3]*A2[index_train3]*exp(-beta1[1]-S2bar[index_train3,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_train3,]%*%alpha1[-1]),family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1[1,index]=rep(mean(Y[index_train3][A1[index_train3]*A2[index_train3]==1]),length(index))
    pred_mu2_new1=rep(mean(Y[index_train3][A1[index_train3]*A2[index_train3]==1]),length(index_train4))
  } else {
    pred_mu2_1[1,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
  }
  #mu2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index_train3,],y=Y[index_train3],weights=(1-A1[index_train3])*(1-A2[index_train3])*exp(-beta0[1]-S2bar[index_train3,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_train3,]%*%alpha0[-1]),family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0[1,index]=rep(mean(Y[index_train3][A1[index_train3]+A2[index_train3]==0]),length(index))
    pred_mu2_new0=rep(mean(Y[index_train3][A1[index_train3]*A2[index_train3]==1]),length(index_train4))
  } else {
    pred_mu2_0[1,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_train4,],s="lambda.min")
  }
  #mu1_1
  Y_DR=pred_mu2_new1+A2[index_train4]*(Y[index_train4]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_train4,]%*%beta1[-1])
  fit=tryCatch(cv.glmnet(x=S1[index_train4,],y=Y_DR,weights=A1[index_train4]*exp(-alpha1[1]-S1[index_train4,]%*%alpha1[-1]),family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1[1,index]=rep(mean(Y_DR[A1[index_train4]==1]),length(index))
  } else {
    pred_mu1_1[1,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  #mu1_0
  Y_DR=pred_mu2_new0+(1-A2[index_train4])*(Y[index_train4]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_train4,]%*%beta0[-1])
  fit=tryCatch(cv.glmnet(x=S1[index_train4,],y=Y_DR,weights=(1-A1[index_train4])*exp(-alpha0[1]-S1[index_train4,]%*%alpha0[-1]),family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0[1,index]=rep(mean(Y_DR[A1[index_train4]==0]),length(index))
  } else {
    pred_mu1_0[1,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
}

# SMDR2
time=proc.time()[3]
pred_pi1_1_0_SMDR2=rep(0,N)
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  
  #pi1_1
  fit=tryCatch(glm.regu.cv(x=S1[index_nk,],y=A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  alpha1=fit$sel.bet[,1]
  if (sum(is.na(alpha1))>0|is.null(fit)) alpha1=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
  pred_pi1_1[2,index]=logistic(alpha1[1]+S1[index,]%*%alpha1[-1]) # P(A1=1|X)
  #pi1_0
  fit=tryCatch(glm.regu.cv(x=S1[index_nk,],y=1-A1[index_nk],fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  alpha0=fit$sel.bet[,1]
  if (sum(is.na(alpha0))>0|is.null(fit)) alpha0=c(-log(1/mean(1-A1[index_nk]) - 1), rep(0,d1))
  pred_pi1_1_0_SMDR2[index]=logistic(alpha0[1]+S1[index,]%*%alpha0[-1]) # P(A0=1|X)
  #pi2_1
  iw=as.numeric(A1[index_nk]/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]))
  fit=tryCatch(glm.regu.cv(x=S2bar[index_nk,],y=A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  beta1=fit$sel.bet[,1]
  if (sum(is.na(beta1))>0|is.null(fit)) beta1=c(-log(mean(iw)/mean(iw*A2[index_nk]) - 1), rep(0,d))
  pred_pi2_1[2,index]=logistic(beta1[1]+S2bar[index,]%*%beta1[-1])  # P(A2=1|V,A1=1)
  #pi2_0
  iw=as.numeric((1-A1[index_nk])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]))
  fit=tryCatch(glm.regu.cv(x=S2bar[index_nk,],y=1-A2[index_nk],iw=iw,fold=5,loss="cal",nrho=nrho,tune.fac=0.01^(1/(nrho-1))),error=function(e) NULL)
  beta0=fit$sel.bet[,1]#P(A2=0|V,A1=0)
  if (sum(is.na(beta0))>0|is.null(fit)) beta0=c(-log(mean(iw)/mean(iw*(1-A2[index_nk])) - 1), rep(0,d))
  pred_pi2_0[2,index]=logistic(-beta0[1]-S2bar[index,]%*%beta0[-1]) # P(A2=1|V,A1=0)=1-P(A2=0|V,A1=0)
  #mu2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=A1[index_nk]*A2[index_nk]*exp(-beta1[1]-S2bar[index_nk,]%*%beta1[-1])/logistic(alpha1[1]+S1[index_nk,]%*%alpha1[-1]),family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1[2,index]=rep(mean(Y[index_nk][A1[index_nk]*A2[index_nk]==1]),length(index))
    pred_mu2_new1=rep(mean(Y[index_nk][A1[index_nk]*A2[index_nk]==1]),length(index_nk))
  } else {
    pred_mu2_1[2,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new1=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
  }
  #mu2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index_nk,],y=Y[index_nk],weights=(1-A1[index_nk])*(1-A2[index_nk])*exp(-beta0[1]-S2bar[index_nk,]%*%beta0[-1])/logistic(alpha0[1]+S1[index_nk,]%*%alpha0[-1]),family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0[2,index]=rep(mean(Y[index_nk][A1[index_nk]+A2[index_nk]==0]),length(index))
    pred_mu2_new0=rep(mean(Y[index_nk][A1[index_nk]+A2[index_nk]==0]),length(index_nk))
  } else {
    pred_mu2_0[2,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_new0=predict(fit,newx=S2bar[index_nk,],s="lambda.min")
  }
  #mu1_1
  Y_DR=pred_mu2_new1+A2[index_nk]*(Y[index_nk]-pred_mu2_new1)/logistic(beta1[1]+S2bar[index_nk,]%*%beta1[-1])
  iw=A1[index_nk]*exp(-alpha1[1]-S1[index_nk,]%*%alpha1[-1])
  fit=tryCatch(cv.glmnet(x=S1[index_nk,],y=Y_DR,weights=iw,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1[2,index]=rep(mean(Y_DR[A1[index_nk]==1]),length(index))
  } else {
    pred_mu1_1[2,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  #mu1_0
  Y_DR=pred_mu2_new0+(1-A2[index_nk])*(Y[index_nk]-pred_mu2_new0)/logistic(beta0[1]+S2bar[index_nk,]%*%beta0[-1])
  iw=(1-A1[index_nk])*exp(-alpha0[1]-S1[index_nk,]%*%alpha0[-1])
  fit=tryCatch(cv.glmnet(x=S1[index_nk,],y=Y_DR,weights=iw,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0[2,index]=rep(mean(Y_DR[A1[index_nk]==0]),length(index))
  } else {
    pred_mu1_0[2,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
}

# DTL1
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  index_nk1=index_nk[1:floor(length(index_nk)/2)]
  index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  index_train1=index_nk1[1:floor(length(index_nk1)/2)]
  index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
  index_train3=index_nk2[1:floor(length(index_nk2)/2)]
  index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
  #pi1
  fit=tryCatch(cv.glmnet(x=S1[index_train1,],y=A1[index_train1],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi1_1[3,index]=rep(mean(A1[index_train1]),length(index))
  } else {
    pred_pi1_1[3,index]=predict(fit,newx=S1[index,],type="response",s="lambda.min")
  }
  #pi2_1
  index_train2_1=index_train2[A1[index_train2]==1]
  fit=tryCatch(cv.glmnet(x=S2bar[index_train2_1,],y=A2[index_train2_1],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_1[3,index]=rep(mean(A2[index_train2_1]),length(index))
  } else {
    pred_pi2_1[3,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #pi2_0
  index_train2_0=index_train2[A1[index_train2]==0]
  fit=tryCatch(cv.glmnet(x=S2bar[index_train2_0,],y=A2[index_train2_0],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_0[3,index]=rep(mean(A2[index_train2_0]),length(index))
  } else {
    pred_pi2_0[3,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #mu2_1
  index_train3_11=index_train3[A1[index_train3]*A2[index_train3]==1]
  index_train4_1=index_train4[A1[index_train4]==1]
  fit=tryCatch(cv.glmnet(x=S2bar[index_train3_11,],y=Y[index_train3_11],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1[3,index]=rep(mean(Y[index_train3_11]),length(index))
    pred_mu2_1_nocross=rep(mean(Y[index_train3_11]),length(index_train4_1))
  } else {
    pred_mu2_1[3,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_1_nocross=predict(fit,newx=S2bar[index_train4_1,],s="lambda.min")
  }
  #mu2_0
  index_train3_00=index_train3[A1[index_train3]+A2[index_train3]==0]
  index_train4_0=index_train4[A1[index_train4]==0]
  fit=tryCatch(cv.glmnet(x=S2bar[index_train3_00,],y=Y[index_train3_00],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0[3,index]=rep(mean(Y[index_train3_00]),length(index))
    pred_mu2_0_nocross=rep(mean(Y[index_train3_00]),length(index_train4_0))
  } else {
    pred_mu2_0[3,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_0_nocross=predict(fit,newx=S2bar[index_train4_0,],s="lambda.min")
  }
  #mu1_1
  fit=tryCatch(cv.glmnet(x=S1[index_train4_1,],y=pred_mu2_1_nocross,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1[3,index]=rep(mean(pred_mu2_1_nocross),length(index))
  } else {
    pred_mu1_1[3,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  #mu1_0
  fit=tryCatch(cv.glmnet(x=S1[index_train4_0,],y=pred_mu2_0_nocross,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0[3,index]=rep(mean(pred_mu2_0_nocross),length(index))
  } else {
    pred_mu1_0[3,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
}

# DTL2
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  index1_nk=index_nk[which(A1[index_nk]==1)];index0_nk=index_nk[which(A1[index_nk]==0)]
  index11_nk=index_nk[which(A1[index_nk]*A2[index_nk]==1)];index00_nk=index_nk[which(A1[index_nk]+A2[index_nk]==0)]
  #pi1
  fit=tryCatch(cv.glmnet(x=S1[index_nk,],y=A1[index_nk],family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi1_1[4,index]=rep(mean(A1[index_nk]),length(index))
  } else {
    pred_pi1_1[4,index]=predict(fit,newx=S1[index,],type="response",s="lambda.min")
  }
  #pi2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index1_nk,],y=A2[index1_nk],family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_1[4,index]=rep(mean(A2[index1_nk]),length(index))
  } else {
    pred_pi2_1[4,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #pi2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index0_nk,],y=A2[index0_nk],family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_0[4,index]=rep(mean(A2[index0_nk]),length(index))
  } else {
    pred_pi2_0[4,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #mu2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index11_nk,],y=Y[index11_nk],family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1[4,index]=rep(mean(Y[index11_nk]),length(index))
    pred_mu2_1_nocross=rep(mean(Y[index11_nk]),length(index1_nk))
  } else {
    pred_mu2_1[4,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_1_nocross=predict(fit,newx=S2bar[index1_nk,],s="lambda.min")
  }
  #mu2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index00_nk,],y=Y[index00_nk],family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0[4,index]=rep(mean(Y[index00_nk]),length(index))
    pred_mu2_0_nocross=rep(mean(Y[index00_nk]),length(index0_nk))
  } else {
    pred_mu2_0[4,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
    pred_mu2_0_nocross=predict(fit,newx=S2bar[index0_nk,],s="lambda.min")
  }
  #mu1_1
  fit=tryCatch(cv.glmnet(x=S1[index1_nk,],y=pred_mu2_1_nocross,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1[4,index]=rep(mean(pred_mu2_1_nocross),length(index))
  } else {
    pred_mu1_1[4,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
  #mu1_0
  fit=tryCatch(cv.glmnet(x=S1[index0_nk,],y=pred_mu2_0_nocross,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0[4,index]=rep(mean(pred_mu2_0_nocross),length(index))
  } else {
    pred_mu1_0[4,index]=predict(fit,newx=S1[index,],s="lambda.min")
  }
}

# IPW
pred_pi1_1[5,]=pred_pi1_1[4,];pred_pi2_1[5,]=pred_pi2_1[4,];pred_pi2_0[5,]=pred_pi2_0[4,]

# S-DRL
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  index_nk1=index_nk[1:floor(length(index_nk)/2)]
  index1_nk1=index_nk1[which(A1[index_nk1]==1)];index0_nk1=index_nk1[which(A1[index_nk1]==0)]
  index11_nk1=index_nk1[which(A1[index_nk1]*A2[index_nk1]==1)];index00_nk1=index_nk1[which(A1[index_nk1]+A2[index_nk1]==0)]
  index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  index1_nk2=index_nk2[which(A1[index_nk2]==1)];index0_nk2=index_nk2[which(A1[index_nk2]==0)]
  index11_nk2=index_nk2[which(A1[index_nk2]*A2[index_nk2]==1)];index00_nk2=index_nk2[which(A1[index_nk2]+A2[index_nk2]==0)]
  index1_nk=index_nk[which(A1[index_nk]==1)];index0_nk=index_nk[which(A1[index_nk]==0)]
  index11_nk=index_nk[which(A1[index_nk]*A2[index_nk]==1)];index00_nk=index_nk[which(A1[index_nk]+A2[index_nk]==0)]
  #pi1
  fit=tryCatch(cv.glmnet(x=S1[index_nk,],y=A1[index_nk],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi1_1[6,index]=rep(mean(A1[index_nk]),length(index))
  } else {
    pred_pi1_1[6,index]=predict(fit,newx=S1[index,],type="response",s="lambda.min")
  }
  #pi2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index1_nk,],y=A2[index1_nk],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_1[6,index]=rep(mean(A2[index1_nk]),length(index))
  } else {
    pred_pi2_1[6,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #pi2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index0_nk,],y=A2[index0_nk],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_0[6,index]=rep(mean(A2[index0_nk]),length(index))
  } else {
    pred_pi2_0[6,index]=predict(fit,newx=S2bar[index,],type="response",s="lambda.min")
  }
  #mu2_1
  fit=tryCatch(cv.glmnet(x=S2bar[index11_nk,],y=Y[index11_nk],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1[6,index]=rep(mean(Y[index11_nk]),length(index))
  } else {
    pred_mu2_1[6,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
  }
  #mu2_0
  fit=tryCatch(cv.glmnet(x=S2bar[index00_nk,],y=Y[index00_nk],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0[6,index]=rep(mean(Y[index00_nk]),length(index))
  } else {
    pred_mu2_0[6,index]=predict(fit,newx=S2bar[index,],s="lambda.min")
  }
  #mu1_1
  fit=tryCatch(cv.glmnet(x=S2bar[index1_nk1,],y=A2[index1_nk1],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_1_new1=rep(mean(A2[index1_nk1]),length(index1_nk2))
  } else {
    pred_pi2_1_new1=predict(fit,newx=S2bar[index1_nk2,],type="response",s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index11_nk1,],y=Y[index11_nk1],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1_new1=rep(mean(Y[index11_nk1]),length(index1_nk2))
  } else {
    pred_mu2_1_new1=predict(fit,newx=S2bar[index1_nk2,],s="lambda.min")
  }
  Y1_new1=as.numeric(pred_mu2_1_new1+A2[index1_nk2]*(Y[index1_nk2]-pred_mu2_1_new1)/pred_pi2_1_new1)
  fit=tryCatch(cv.glmnet(x=S1[index1_nk2,],y=Y1_new1,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1_new1=rep(mean(Y1_new1),length(index))
  } else {
    pred_mu1_1_new1=predict(fit,newx=S1[index,],s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index1_nk2,],y=A2[index1_nk2],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_1_new2=rep(mean(A2[index1_nk2]),length(index1_nk1))
  } else {
    pred_pi2_1_new2=predict(fit,newx=S2bar[index1_nk1,],type="response",s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index11_nk2,],y=Y[index11_nk2],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_1_new2=rep(mean(Y[index11_nk2]),length(index1_nk1))
  } else {
    pred_mu2_1_new2=predict(fit,newx=S2bar[index1_nk1,],s="lambda.min")
  }
  Y1_new2=as.numeric(pred_mu2_1_new2+A2[index1_nk1]*(Y[index1_nk1]-pred_mu2_1_new2)/pred_pi2_1_new2)
  fit=tryCatch(cv.glmnet(x=S1[index1_nk1,],y=Y1_new2,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1_new2=rep(mean(Y1_new2),length(index))
  } else {
    pred_mu1_1_new2=predict(fit,newx=S1[index,],s="lambda.min")
  }
  pred_mu1_1[6,index]=(pred_mu1_1_new1+pred_mu1_1_new2)/2
  #mu1_0
  fit=tryCatch(cv.glmnet(x=S2bar[index0_nk1,],y=A2[index0_nk1],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_0_new1=rep(mean(A2[index0_nk1]),length(index0_nk2))
  } else {
    pred_pi2_0_new1=predict(fit,newx=S2bar[index0_nk2,],type="response",s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index00_nk1,],y=Y[index00_nk1],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu2_0_new1=rep(mean(Y[index00_nk1]),length(index0_nk2))
  } else {
    pred_mu2_0_new1=predict(fit,newx=S2bar[index0_nk2,],s="lambda.min")
  }
  Y0_new1=as.numeric(pred_mu2_0_new1+(1-A2[index0_nk2])*(Y[index0_nk2]-pred_mu2_0_new1)/(1-pred_pi2_0_new1))
  fit=tryCatch(cv.glmnet(x=S1[index0_nk2,],y=Y0_new1,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0_new1=rep(mean(Y0_new1),length(index))
  } else {
    pred_mu1_0_new1=predict(fit,newx=S1[index,],s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index0_nk2,],y=A2[index0_nk2],family="binomial",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_pi2_0_new2=rep(mean(A2[index0_nk2]),length(index0_nk1))
  } else {
    pred_pi2_0_new2=predict(fit,newx=S2bar[index0_nk1,],type="response",s="lambda.min")
  }
  fit=tryCatch(cv.glmnet(x=S2bar[index00_nk2,],y=Y[index00_nk2],family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_1_new2=rep(mean(Y[index00_nk2]),length(index0_nk1))
  } else {
    pred_mu2_0_new2=predict(fit,newx=S2bar[index0_nk1,],s="lambda.min")
  }
  Y0_new2=as.numeric(pred_mu2_0_new2+(1-A2[index0_nk1])*(Y[index0_nk1]-pred_mu2_0_new2)/(1-pred_pi2_0_new2))
  fit=tryCatch(cv.glmnet(x=S1[index0_nk1,],y=Y0_new2,family="gaussian",nfolds=5),error=function(e) NULL)
  if (is.null(fit)) {
    pred_mu1_0_new2=rep(mean(Y0_new2),length(index))
  } else {
    pred_mu1_0_new2=predict(fit,newx=S1[index,],s="lambda.min")
  }
  pred_mu1_0[6,index]=(pred_mu1_0_new1+pred_mu1_0_new2)/2
}

# SDR+
#pi1
fit=tryCatch(cv.glmnet(x=S1,y=A1,family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_pi2_1[7,]=rep(mean(A1),N)
} else {
  pred_pi1_1[7,]=predict(fit,newx=S1,type="response",s="lambda.min")
}
#pi2_1
fit=tryCatch(cv.glmnet(x=S2bar[A1==1,],y=A2[A1==1],family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_pi2_1[7,]=rep(mean(A2[A1==1]),N)
} else {
  pred_pi2_1[7,]=predict(fit,newx=S2bar,type="response",s="lambda.min")
}
#pi2_0
fit=tryCatch(cv.glmnet(x=S2bar[A1==0,],y=A2[A1==0],family="binomial",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_pi2_0[7,]=rep(mean(A2[A1==0]),N)
} else {
  pred_pi2_0[7,]=predict(fit,newx=S2bar,type="response",s="lambda.min")
}
#mu2_1
fit=tryCatch(cv.glmnet(x=S2bar[A1*A2==1,],y=Y[A1*A2==1],family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_mu2_1[7,]=rep(mean(Y[A1*A2==1]),N)
} else {
  pred_mu2_1[7,]=predict(fit,newx=S2bar,s="lambda.min")
}
#mu2_0
fit=tryCatch(cv.glmnet(x=S2bar[A1+A2==0,],y=Y[A1+A2==0],family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_mu2_0[7,]=rep(mean(Y[A1+A2==0]),N)
} else {
  pred_mu2_0[7,]=predict(fit,newx=S2bar,s="lambda.min")
}
#mu1_1
Y_DR=pred_mu2_1[7,A1==1]+A2[A1==1]*(Y[A1==1]-pred_mu2_1[7,A1==1])/pred_pi2_1[7,A1==1]
fit=tryCatch(cv.glmnet(x=S1[A1==1,],y=Y_DR,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_mu1_1[7,]=rep(mean(Y_DR),N)
} else {
  pred_mu1_1[7,]=predict(fit,newx=S1,s="lambda.min")
}
#mu1_0
Y_DR=pred_mu2_0[7,A1==0]+(1-A2[A1==0])*(Y[A1==0]-pred_mu2_0[7,A1==0])/(1-pred_pi2_0[7,A1==0])
fit=tryCatch(cv.glmnet(x=S1[A1==0,],y=Y_DR,family="gaussian",nfolds=5,nlambda=nrho),error=function(e) NULL)
if (is.null(fit)) {
  pred_mu1_0[7,]=rep(mean(Y_DR),N)
} else {
  pred_mu1_0[7,]=predict(fit,newx=S1,s="lambda.min")
}

# SDR-RF
time=proc.time()[3]
for (k in 1:K){
  if (k==K) {
    index=((k-1)*M+1):N
  } else {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  index1_nk=index_nk[which(A1[index_nk]==1)];index0_nk=index_nk[which(A1[index_nk]==0)]
  index11_nk=index_nk[which(A1[index_nk]*A2[index_nk]==1)];index00_nk=index_nk[which(A1[index_nk]+A2[index_nk]==0)]
  #pi1
  errs=numeric() # Tunning through OOB Errors
  for (j in 1:nrow(paras)){
    fit=ranger(y = A1[index_nk], x = data.frame(X=S1[index_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = A1[index_nk], x = data.frame(X=S1[index_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_pi=predict(fit, data = data.frame(X=S1[index,]))$predictions[,1]
  pred_pi[which(pred_pi==0)]=0.01;pred_pi[which(pred_pi==1)]=0.99 # Truncate RF PS estimates to [0.01,0.99]
  pred_pi1_1[8,index]=pred_pi
  #pi2_1
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = A2[index1_nk], x = data.frame(X=S2bar[index1_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = A2[index1_nk], x = data.frame(X=S2bar[index1_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_pi=predict(fit, data = data.frame(X=S2bar[index,]))$predictions[,1]
  pred_pi[which(pred_pi==0)]=0.01;pred_pi[which(pred_pi==1)]=0.99
  pred_pi2_1[8,index]=pred_pi
  pred_pi=predict(fit, data = data.frame(X=S2bar[index1_nk,]))$predictions[,1]
  pred_pi[which(pred_pi==0)]=0.01;pred_pi[which(pred_pi==1)]=0.99
  pred_pi2_1_nocross=pred_pi
  #pi2_0
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = A2[index0_nk], x = data.frame(X=S2bar[index0_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = A2[index0_nk], x = data.frame(X=S2bar[index0_nk,]), num.trees = 200, verbose = FALSE, probability = TRUE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_pi=predict(fit, data = data.frame(X=S2bar[index,]))$predictions[,1]
  pred_pi[which(pred_pi==0)]=0.01;pred_pi[which(pred_pi==1)]=0.99
  pred_pi2_0[8,index]=pred_pi
  pred_pi=predict(fit, data = data.frame(X=S2bar[index0_nk,]))$predictions[,1]
  pred_pi[which(pred_pi==0)]=0.01;pred_pi[which(pred_pi==1)]=0.99
  pred_pi2_0_nocross=pred_pi
  #mu2_1
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = Y[index11_nk], x = data.frame(X=S2bar[index11_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = Y[index11_nk], x = data.frame(X=S2bar[index11_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_mu2_1[8,index]=predict(fit, data = data.frame(X=S2bar[index,]))$predictions
  pred_mu2_1_nocross=predict(fit, data = data.frame(X=S2bar[index1_nk,]))$predictions
  #mu2_0
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = Y[index00_nk], x = data.frame(X=S2bar[index00_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = Y[index00_nk], x = data.frame(X=S2bar[index00_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_mu2_0[8,index]=predict(fit, data = data.frame(X=S2bar[index,]))$predictions
  pred_mu2_0_nocross=predict(fit, data = data.frame(X=S2bar[index0_nk,]))$predictions
  #mu1_1
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = as.numeric(pred_mu2_1_nocross+A2[index1_nk]*(Y[index1_nk]-pred_mu2_1_nocross)/pred_pi2_1_nocross), x = data.frame(X=S1[index1_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = as.numeric(pred_mu2_1_nocross+A2[index1_nk]*(Y[index1_nk]-pred_mu2_1_nocross)/pred_pi2_1_nocross), x = data.frame(X=S1[index1_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_mu1_1[8,index]=predict(fit, data = data.frame(X=S1[index,]))$predictions
  #mu1_0
  errs=numeric()
  for (j in 1:nrow(paras)){
    fit=ranger(y = as.numeric(pred_mu2_0_nocross+(1-A2[index0_nk])*(Y[index0_nk]-pred_mu2_0_nocross)/(1-pred_pi2_0_nocross)), x = data.frame(X=S1[index0_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j,1], min.node.size = paras[j,2])
    errs=c(errs,fit$prediction.error)
  }
  j_opt=which.min(errs)
  fit=ranger(y = as.numeric(pred_mu2_0_nocross+(1-A2[index0_nk])*(Y[index0_nk]-pred_mu2_0_nocross)/(1-pred_pi2_0_nocross)), x = data.frame(X=S1[index0_nk,]), num.trees = 200, verbose = FALSE, mtry = paras[j_opt,1], min.node.size = paras[j_opt,2])
  pred_mu1_0[8,index]=predict(fit, data = data.frame(X=S1[index,]))$predictions
}

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
  psi_0=gamma2_0*Y-(gamma1_0-1)*pred_mu1_0[j,]-(gamma2_0-gamma1_0)*pred_mu2_0[j,]
  psi=psi_1-psi_0-pred_theta_t[j]
  pred_sig2_t=c(pred_sig2_t,mean(psi^2))
}

# Empdiff
pred_theta_t=c(pred_theta_t,mean(Y[which(A1*A2==1)])-mean(Y[which(A1+A2==0)]))
pred_sig2_t=c(pred_sig2_t,var(Y[which(A1*A2==1)])*N/sum(A1*A2==1)+var(Y[which(A1+A2==0)])*N/sum(A1+A2==0))

pred_theta=pred_theta_t;pred_sig2=pred_sig2_t

SE=sqrt(pred_sig2/N)
p_value=2*pnorm(-abs(pred_theta)/SE)
CI_lower=pred_theta-1.96*SE
CI_upper=pred_theta+1.96*SE

# Report results in a Latex-friendly format
est=c("SMDR1","SMDR2","DTL1","DTL2","IPW","S-DRL","SDR+","SDR-RF","empdiff")
table_result=function(x){
  paste(est[x],"&",format(round(pred_theta[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(SE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(CI_lower[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(CI_upper[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(p_value[x],digits=3),nsmall=3,scientific=FALSE),sep="")
}
result=apply(t(1:num_totest),2,table_result)
result=list(result=c("estimator&theta&SE&CI_lower&CI_upper&p_value",result))
result 

proc.time()-ptm
