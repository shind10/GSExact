#' GSExact function
#' This is Testing. It's my very first time uploading package onto Github


GSExact_basic=function(input_data, J,n, burn, iter){

  final_data = input_data

  x1 = final_data$x1[seq(1, J*n, n)]
  x2 = final_data$x2[seq(1, J*n, n)]
  x3 = final_data$x3[seq(1, J*n, n)]
  vij = final_data$vij

  x1_missing=ifelse(is.na(x1)==TRUE, 1, 0 )
  x2_missing=ifelse(is.na(x2)==TRUE, 1, 0 )
  vm=ifelse(is.na(vij)==TRUE, 1, 0)
  ##########################

  init_lmer=lmer(vij~x1+x2+x3+x1:x2+(1|id), data=final_data)
  lme=summary(init_lmer)

  init_tau = lme$varcor$id[1]; Tau=init_tau
  init_sig = lme$sigma^2; Sig=init_sig

  B0=lme$coefficients[1,1]  #initial values for beta
  B1=lme$coefficients[2,1]
  B2=lme$coefficients[3,1]
  B3=lme$coefficients[4,1]
  B4=lme$coefficients[5,1]

  Beta=matrix(c(B0,B1,B2,B3,B4),5,1)

  ## initial Alpha and initial Tx.star . for bivariate X1, X2 | X3 ~

  a0=summary(lm(x1~x3))$coefficient[1,1]
  a1=summary(lm(x1~x3))$coefficient[2,1]
  a2=summary(lm(x2~x3))$coefficient[1,1]
  a3=summary(lm(x2~x3))$coefficient[2,1]

  Alpha=matrix(c(a0,a1,a2,a3),4,1); Alpha # alpha will be used for mean.

  vx1_3=var(resid(lm(x1~x3, data=final_data)))
  vx2_3=var(resid(lm(x2~x3, data=final_data)))
  covx12_3=summary(lm(x1~x2+x3, data=final_data))$coefficients[2,1]*var(resid(lm(x2~x3, data=final_data)))

  Tx.star=S0=matrix(c(vx1_3, covx12_3, covx12_3, vx2_3), 2, 2)

  ##initial value how?

  final_data$x1=ifelse(is.na(final_data$x1)==TRUE, mean(x1, na.rm=T), final_data$x1)
  final_data$x2=ifelse(is.na(final_data$x2)==TRUE, mean(x2, na.rm=T), final_data$x2)

  ybar=aggregate(final_data$vij, by=list(final_data$id), mean, na.rm=T)$x
  ybar2=ifelse(is.na(ybar)==TRUE, mean(final_data$vij, na.rm=TRUE), ybar)
  yvec=rep(ybar2, each=4)
  final_data$vij=ifelse(is.na(final_data$vij)==TRUE, yvec, final_data$vij)

  x1=final_data$x1[seq(1, dim(final_data)[1],4)]
  x2=final_data$x2[seq(1, dim(final_data)[1],4)]
  vij=final_data$vij

  ###variance covariance matrix check
  #True_Tmat

  ## Priors

  sa=1
  sb=0.5

  ta=1
  tb=0.5


  save.beta0=save.beta1=save.beta2=save.beta3=save.beta4=save.beta5=save.beta6=NULL
  save.tau=save.sig=NULL
  save.Tmat=array(0, dim=c(2,2,iter))
  save.Alpha=array(0, dim=c(4,1,iter))


  for (k in 1:iter){

    ## step 1 : Uj

    data_list = split(final_data$vij, f = final_data$id) #vij = final_data$vij. make sure.
    vobs_mat = vector(mode = "list", length = J)

    for (i in 1:J){
      vobs_mat[[i]]=matrix(data_list[[i]],n,1)
    }

    Bmat = matrix(c(B0,B1,B2,B3,B4),5,1)
    one = rep(1,J)
    Xmat = rbind(one, x1, x2, x3, x1*x2)

    # joint distribution of uj and vj

    dj_mat=NULL
    uj_mat=NULL
    delj=n*(1/Sig)+(1/Tau); varj=1/delj

    for (j in 1:J){
      one=matrix(1,4,1)
      dj_mat=vobs_mat[[j]] - one%*%t(Xmat[,j])%*%Bmat
      meanj=varj*(1/Sig)*sum(dj_mat)
      uj_mat[j]= rnorm(1,meanj, sqrt(varj))
    }

    uj=uj_mat

    final_data$uj_vec=rep(uj, each=4)

    ##########################

    ## step 2 : Tau  ## Lets do it with Dr.Shin's derivation. 1031 there is little discrapency

    #  new_ta=ta+0.5*J;
    #  new_tb=tb+0.5*sum(uj^2);

    #  tauinv=rgamma(1,new_ta, new_tb)
    #  Tau=1/tauinv;

    new_ta=ta+0.5*J;
    new_tb=solve(sum(uj^2)/2+(1/tb))
    Tau=rinvgamma(1, new_ta, scale=new_tb)
    save.tau=rbind(save.tau, Tau)

    ## step 3 : Draw Beta Directly

    Xj_mat = vector(mode = "list", length = J)

    for (i in 1:J){
      Xj_mat[[i]] = matrix(matrix(1,4,1)%*%Xmat[,i], 4, 5)
    }

    sum_xx=matrix(0,5,5)
    sum_xm=matrix(0,5,1)

    for (j in 1:J){
      sum_xx=sum_xx+t(Xj_mat[[j]])%*%Xj_mat[[j]]
      sum_xm=sum_xm+t(Xj_mat[[j]])%*%(vobs_mat[[j]]-matrix(1,4,1)%*%uj[j])
    }

    Beta_mean=solve(sum_xx)%*%sum_xm  # beta hat in SY2002; estimated Beta.
    Beta_var=Sig*solve(sum_xx)

    Beta=mvrnorm(1, Beta_mean, Beta_var)

    B0=Beta[1]; B1=Beta[2]; B2=Beta[3]; B3=Beta[4];B4=Beta[5]

    save.beta0=rbind(save.beta0, B0)
    save.beta1=rbind(save.beta1, B1)
    save.beta2=rbind(save.beta2, B2)
    save.beta3=rbind(save.beta3, B3)
    save.beta4=rbind(save.beta4, B4)

    ### Step 4 : Sigma
    ### calculate eij ###

    XB = NULL

    for (j in 1:J){
      XB[j]=t(Xmat[,j])%*%Beta
    }

    final_data$XB_vec=rep(XB, each=4)
    final_data$eij=final_data$vij-final_data$XB_vec-final_data$uj_vec

    #new_sa=sa+(n*J)/2
    #new_sb=sb+0.5*sum(final_data$eij^2)
    new_sa=sa+(n*J)/2
    new_sb=solve((1/sb)+0.5*sum(final_data$eij^2))
    Sig=rinvgamma(1,new_sa,scale=new_sb)
    save.sig=rbind(save.sig, Sig)

    ################################

    ### Step 5 : making yij
    ### Not use Mvj but we are using updated x1 and x2
    e_ij=rnorm((n*J),0,sd=sqrt(Sig))
    vij=ifelse(vm==1, final_data$XB_vec+final_data$uj_vec+e_ij, final_data$vij)
    final_data$vij=vij

    ##### Bivariate distribution of X1,X2|X3  #####

    Y=matrix(cbind(x1,x2), J, 2)
    ax3=as.matrix(cbind(rep(1,J),x3,rep(0,J),rep(0,J)),4,J)
    x3a=as.matrix(cbind(rep(0,J),rep(0,J),rep(1,J),x3),4,J)

    AX_mat=vector(mode = "list", length = J)

    for (i in 1:J){
      AX_mat[[i]]=matrix(rbind(ax3[i,], x3a[i,]),2,4)
    }

    sum_xx=matrix(0,4,4)
    sum_xm=matrix(0,4,1)

    for (j in 1:J){
      sum_xx=sum_xx+t(AX_mat[[j]])%*%solve(Tx.star)%*%AX_mat[[j]]
      sum_xm=sum_xm+t(AX_mat[[j]])%*%t(solve(Tx.star))%*%matrix(t(Y[j,]),2,1) #
    }

    alpha_mean=solve(sum_xx)%*%sum_xm
    alpha_var=solve(sum_xx)
    Alpha=mvrnorm(1, alpha_mean, alpha_var);
    save.Alpha[,,k]=matrix(Alpha,4,1)


    a0=Alpha[1];a1=Alpha[2];a2=Alpha[3];a3=Alpha[4];

    ## check with the simulated value.

    Mu_x1=matrix(0,2,J)
    for (j in 1:J){
      Mu_x1[,j]=AX_mat[[j]]%*%matrix(Alpha,4,1)
    }

    nu0=4
    Sn=S0+(t(Y)-Mu_x1)%*%t(t(Y)-Mu_x1)
    Tmat=solve(rwish(nu0+J,solve(Sn)))
    Tx.star=Tmat
    save.Tmat[,,k]=Tmat

    ### Now, We have Mmat and Tmat. 2*1 and 2*2 8/25
    ### Now, We have Mmat and Tmat. 2*1 and 2*2 8/25
    Mx1=Mu_x1[1,]; Mx2=Mu_x1[2,];  # simply a0+a1x3 and a2+a3x3  hence, E(x1|x2 x3), E(x2|x1 x3)

    Tii=Tmat[1,1];Tss=Tmat[2,2];Tis=Tmat[1,2];Tsi=Tmat[2,1]

    Mu_x1x2=Mx1+Tis*(1/Tss)*(x2-Mx2)
    T_x1x2=Tii-Tis*(1/Tss)*Tsi  # var(x1|x2x3)

    Mu_x2x1=Mx2+Tsi*(1/Tii)*(x1-Mx1)
    T_x2x1=Tss-Tsi*(1/Tii)*Tis        # var(x2|x1x3)
    ##### 1017 now, new thing is implemented. yj,x1 | x2 x3 uj
    ## First Draw the Conditional Bivariate Distribution of Yj and X1j | X2j X3j Uj
    # mean x1|x2 x3 we already have it
    # and we are using updated yj
    data_list = split(final_data$vij, f = final_data$id) #vij = final_data$vij. make sure.
    vobs_mat = vector(mode = "list", length = J)

    for (i in 1:J){
      vobs_mat[[i]]=matrix(data_list[[i]],n,1)
    }

    meany=B0+B2*x2+B3*x3+uj+(B1+B4*x2)*Mu_x1x2
    meanx=Mu_x1x2

    v11=(B1+B4*x2)^2*T_x1x2; v12=(B1+B4*x2)*T_x1x2; v22=T_x1x2

    one4=matrix(1,4,1)
    jmat=matrix(1,4,4)

    vmat=array(0, dim=c(n+1,n+1,J))

    for (i in 1:J){
      vmat[,,i][1:n,1:n]=jmat*v11[i]+diag(n)*Sig ## 4 x 4
      vmat[,,i][1:n,1+n]=one4*v12[i]   # 4 * 1
      vmat[,,i][n+1,1:n]=t(one4)*v12[i] # 1 * 4
      vmat[,,i][n+1,n+1]=v22 # 1 * 1
    }   # each Vmat is different for cluster to cluster

    ## impute x1 based on x1|y x2 x3 uj
    imp_x1=NULL

    for (j in 1:J){
      if(x1_missing[j]==0){imp_x1[j]=x1[j]}
      else{
        mean_j=meanx[j]+vmat[,,j][n+1,1:n]%*%solve(vmat[,,j][1:n,1:n])%*%(vobs_mat[[j]]-one4*meany[j])
        var_j=vmat[,,j][n+1,n+1]-vmat[,,j][n+1,1:n]%*%solve(vmat[,,j][1:n,1:n])%*%vmat[,,j][1:n,1+n]
        imp_x1[j]=rnorm(1, mean=mean_j, sd=sqrt(var_j))
      }
    }

    x1=imp_x1 #update x1

    final_data$x1=rep(x1, each=4) #update in data_setn


    ############Now is the distribution of x2|x1 y uj ####################################

    meany=B0+B1*x1+B3*x3+uj+(B2+B4*x1)*Mu_x2x1
    meanx=Mu_x2x1

    v11=(B2+B4*x1)^2*T_x2x1; v12=(B2+B4*x1)*T_x2x1; v22=T_x2x1
    one4=matrix(1,4,1)
    jmat=matrix(1,4,4)

    vmat=array(0, dim=c(n+1,n+1,J))

    for (i in 1:J){
      vmat[,,i][1:n,1:n]=jmat*v11[i]+diag(n)*Sig ## 4 x 4
      vmat[,,i][1:n,1+n]=one4*v12[i]   # 4 * 1
      vmat[,,i][n+1,1:n]=t(one4)*v12[i] # 1 * 4
      vmat[,,i][n+1,n+1]=v22 # 1 * 1
    }   # each Vmat is different for cluster to cluster

    ## impute x1 based on x1|y x2 x3 uj
    imp_x2=NULL

    for (j in 1:J){
      if(x2_missing[j]==0){imp_x2[j]=x2[j]}
      else{
        mean_j=meanx[j]+vmat[,,j][n+1,1:n]%*%solve(vmat[,,j][1:n,1:n])%*%(vobs_mat[[j]]-one4*meany[j])
        var_j=vmat[,,j][n+1,n+1]-vmat[,,j][n+1,1:n]%*%solve(vmat[,,j][1:n,1:n])%*%vmat[,,j][1:n,1+n]
        imp_x2[j]=rnorm(1, mean=mean_j, sd=sqrt(var_j))
      }
    }

    x2=imp_x2 #update x1

    final_data$x2=rep(x2, each=4) #update in data_setn

  }

  mc.tau=save.tau
  mc.sig=save.sig
  mc.beta0=save.beta0
  mc.beta1=save.beta1
  mc.beta2=save.beta2
  mc.beta3=save.beta3
  mc.beta4=save.beta4

  mc.mat=matrix(rbind(mc.tau, mc.sig, mc.beta0, mc.beta1, mc.beta2, mc.beta3, mc.beta4), 7, iter)

  mean.tau=mean(save.tau[burn:iter])
  mean.sig=mean(save.sig[burn:iter])
  mean.beta0=mean(save.beta0[burn:iter])
  mean.beta1=mean(save.beta1[burn:iter])
  mean.beta2=mean(save.beta2[burn:iter])
  mean.beta3=mean(save.beta3[burn:iter])
  mean.beta4=mean(save.beta4[burn:iter])

  mean_mat=matrix(c(mean.tau,mean.sig,mean.beta0,mean.beta1,mean.beta2,mean.beta3,mean.beta4),7,1)

  median.tau=median(save.tau[burn:iter])
  median.sig=median(save.sig[burn:iter])
  median.beta0=median(save.beta0[burn:iter])
  median.beta1=median(save.beta1[burn:iter])
  median.beta2=median(save.beta2[burn:iter])
  median.beta3=median(save.beta3[burn:iter])
  median.beta4=median(save.beta4[burn:iter])

  mean_mat=matrix(c(median.tau,median.sig,median.beta0,median.beta1,median.beta2,median.beta3,median.beta4),7,1)

  mean.Tmat=apply(save.Tmat[,,burn:iter], c(1,2), mean)
  mean.Alpha=apply(save.Alpha[,,burn:iter], c(1), mean)

  var.tau=var(save.tau[burn:iter])
  var.sig=var(save.sig[burn:iter])
  var.beta0=var(save.beta0[burn:iter])
  var.beta1=var(save.beta1[burn:iter])
  var.beta2=var(save.beta2[burn:iter])
  var.beta3=var(save.beta3[burn:iter])
  var.beta4=var(save.beta4[burn:iter])

  var_mat=matrix(c(var.tau,var.sig,var.beta0,var.beta1,var.beta2,var.beta3,var.beta4),7,1)

  cp.tau=ifelse(quantile(save.tau[burn:iter], c(0.025,0.975))[1]<4 & 4<quantile(save.tau[burn:iter], c(0.025,0.975))[2],1,0)
  cp.sig=ifelse(quantile(save.sig[burn:iter], c(0.025,0.975))[1]<16 & 16<quantile(save.sig[burn:iter], c(0.025,0.975))[2],1,0)
  cp.beta0=ifelse(quantile(save.beta0[burn:iter], c(0.025,0.975))[1]<1 & 1<quantile(save.beta0[burn:iter], c(0.025,0.975))[2],1,0)
  cp.beta1=ifelse(quantile(save.beta1[burn:iter], c(0.025,0.975))[1]<1 & 1<quantile(save.beta1[burn:iter], c(0.025,0.975))[2],1,0)
  cp.beta2=ifelse(quantile(save.beta2[burn:iter], c(0.025,0.975))[1]<1 & 1<quantile(save.beta2[burn:iter], c(0.025,0.975))[2],1,0)
  cp.beta3=ifelse(quantile(save.beta3[burn:iter], c(0.025,0.975))[1]<1 & 1<quantile(save.beta3[burn:iter], c(0.025,0.975))[2],1,0)
  cp.beta4=ifelse(quantile(save.beta4[burn:iter], c(0.025,0.975))[1]<1 & 1<quantile(save.beta4[burn:iter], c(0.025,0.975))[2],1,0)

  cp_mat=matrix(c(cp.tau, cp.sig, cp.beta0, cp.beta1, cp.beta2, cp.beta3, cp.beta4),7,1)

  ci.tau=quantile(save.tau[burn:iter], c(0.025,0.975))
  ci.sig=quantile(save.sig[burn:iter], c(0.025,0.975))
  ci.beta0=quantile(save.beta0[burn:iter], c(0.025,0.975))
  ci.beta1=quantile(save.beta1[burn:iter], c(0.025,0.975))
  ci.beta2=quantile(save.beta2[burn:iter], c(0.025,0.975))
  ci.beta3=quantile(save.beta3[burn:iter], c(0.025,0.975))
  ci.beta4=quantile(save.beta4[burn:iter], c(0.025,0.975))

  ci_mat=matrix( rbind(ci.tau, ci.sig, ci.beta0, ci.beta1, ci.beta2, ci.beta3, ci.beta4), 7,2)

  return(list(mean_mat, var_mat, mean.Tmat, mean.Alpha, cp_mat, ci_mat, mc.mat))
}

