# Mixture multigroup factor analysis for loadings and intercept (non-)variance
# ----------------------------------------------------------------
# Code written by Kim De Roover
# This version builds on configural invariance across all groups (same zero loadings in case of CFA, same number of factors in case of EFA) and deals with differences in loadings and intercepts simultaneously
# (finds clusters of groups based on similarity of their loadings and intercepts, given the user-specified number of clusters)
# For model selection, it is advised to use BIC_G (number of groups as sample size) in combination with CHull (see preprint)
# Please cite publications: https://www.tandfonline.com/doi/full/10.1080/10705511.2020.1866577
# and https://doi.apa.org/doi/10.1037/met0000355

# INPUT:
# Xsup = data matrix for all groups (rows are subjects nested within groups, columns are the variables to be factor-analyzed)
# N_gs = vector with number of subjects for each group (in the same order as they appear in the data matrix)
# nclust = user-specified number of clusters
# nfactors = user-specified number of factors
# Maxiter = maximum number of iterations
# nruns = number of starts (based on pre-selected random partitions when start = 1)
# preselect = percentage of best starts taken in pre-selection (increase to speed up startprocedure)
# design = matrix indicating position of zero loadings with '0' and non-zero loadings with '1' (specify for CFA, leave unspecified for EFA)
#          (using different design matrices for different clusters is currently not supported)
# startpartition = partition of groups to start from (use with start = 2 and nruns = 1)

# OUTPUT:
# z_gks = cluster memberships of groups (posterior classification probabilities)
# pi_ks= mixing proportions (prior classification probabilities)
# Lambda_ks = cluster-specific loadings, access loadings of cluster k via Lambda_ks[[k]]
# Psi_gs = group-specific unique variances, access loadings of group g via Psi_gs[[g]]
# Phi_gks = group- and cluster-specific factor (co)variances, access (co)variances of group g in cluster k via Phi_gks[[g,k]]
# tau_ks = group-specific means, access intercepts of cluster k via tau_ks[k,]
# alpha_gks = group- and cluster-specific factor means, access factor means of group g in cluster k via alpha_gks[[g,k]]
# bestloglik = loglikelihood of best start
# logliks = loglikelihoods of all starts
# nrpars = number of free parameters, to be used for model selection in combination with bestloglik
# convergence = 2 if converged on loglikelihood, 1 if converged on parameter changes, 0 if not converged
# nractivatedconstraints = number of constraints on the unique variances (across groups) to avoid unique variances approaching zero

MixtureMG_FA_loadingsandintercepts <- function(Xsup,N_gs,nclust,nfactors,Maxiter = 1000,start = 1,nruns = 50,design = 0,preselect = 10,startpartition){
  
  Xsup=as.matrix(Xsup)
  ngroup <- length(N_gs)
  if(nrow(N_gs)!=ngroup || is.null(nrow(N_gs))){ # make sure N_gs is a column vector
    N_gs_colvec=matrix(0,ngroup,1)
    for (g in 1:ngroup){
      N_gs_colvec[g,]=N_gs[g]
    }
    N_gs <- N_gs_colvec
  }
  nvar <- ncol(Xsup)
  N <- sum(N_gs);
  IM <- diag(nclust)
  Ncum <- matrix(0,ngroup,2)
  Ncum[1,1]=1 
  Ncum[1,2]=N_gs[1]
  for(g in 2:ngroup){
    Ncum[g,1]=sum(N_gs[1:(g-1)])+1 # Ncum[g,1]: first row in Xsup for group g
    Ncum[g,2]=sum(N_gs[1:g]) # Ncum[g,2]: last row in Xsup for group g
  }
  
  if (sum(design)==0){ # if design is unspecified, EFA is used
    design=matrix(1,nvar,nfactors)
    EFA=1
  } else {
    EFA=0
  }
  
  
  # compute group-specific means
  mu_gs=matrix(0,ngroup,nvar)
  for(g in 1:ngroup){
    X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
    mu_gs[g,] <- apply(X_g,2,mean)
  }
  
  
  if(start==1){
    if(nclust>1){
      # pre-selection of random partitions
      nrtrialstarts=nruns*round(100/preselect) # generate 'nruns'*(100/preselect) different random partitions
      randpartvecs=matrix(0,nrtrialstarts,ngroup);
      for (trialstart in 1:nrtrialstarts){
        aris=1;
        iterrnd=0;
        while (sum(aris==1)>0 && iterrnd<5){
          cl=0;
          while(length(cl)<nclust){
            randpartvec <- sample(1:nclust,ngroup,replace=TRUE) # generate random partition
            cl=unique(randpartvec)
          }
          nstartprev=trialstart-1
          aris=matrix(0,nstartprev,1)
          if (nstartprev>0){
            for (r in 1:nstartprev){
              prevpartvec=randpartvecs[r,]
              aris[r]<-adjrandindex(prevpartvec,randpartvec)
            }
            iterrnd=iterrnd+1;
          }
        }
        randpartvecs[trialstart,]=randpartvec
      }
      ODLLs_trialstarts=rep(0,nrtrialstarts,1)
      for (trialstart in 1:nrtrialstarts){ # select 'preselect'% (default 10%) best fitting random partitions
        randpartvec=randpartvecs[trialstart,];
        if(nclust>1){
          z_gks=IM[randpartvec,]
          pi_ks=(1/ngroup)*apply(z_gks,2,sum)
        }
        else {
          z_gks=t(randpartvec)
          pi_ks=1
        }
        N_gks=diag(N_gs[,1])%*%z_gks
        N_ks=apply(N_gks,2,sum)
        
        tau_ks <- matrix(0,nclust,nvar)
        for(k in 1:nclust){
          Xsup_k=matrix(0,N_ks[k],nvar)
          for(g in 1:ngroup){
            if(randpartvec[g]==k){
              X=Xsup[Ncum[g,1]:Ncum[g,2],]
              Xsup_k[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=X
            }
          }
          tau_ks[k,] <- apply(Xsup_k,2,mean)
        }
        # center data with initialized intercepts
        Xsupcent <- matrix(0,N,nvar)
        for(g in 1:ngroup){
          k=randpartvec[g]
          X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
          Xsupcent[Ncum[g,1]:Ncum[g,2],]=X_g-t(matrix(tau_ks[k,],ncol=N_gs[g],nrow=nvar))
        }
        
	      Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
        uniq_ks <- matrix(0,nclust,nvar)
	      for(k in 1:nclust){
          Xsup_kc=matrix(0,N_ks[k],nvar)
          for(g in 1:ngroup){
            if(randpartvec[g]==k){
              Xc=Xsupcent[Ncum[g,1]:Ncum[g,2],]
              Xsup_kc[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=Xc
            }
          }
          S_k <- (1/N_ks[k])*(t(Xsup_kc)%*%Xsup_kc)
          ed<-eigen(S_k, symmetric=TRUE, only.values = FALSE)
          val<-ed$values
          u<-ed$vectors
          totalerror=sum((val[-seq_len(nfactors)]))
          meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
          Uniq=rep(meanerror,nvar)
          lambda_k=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
          if (EFA==0){
            lambda_k <- procr(lambda_k,design)
            lambda_k=lambda_k*design # non-zero loadings should be indicated with '1' for this to work properly
          }
          Lambda_ks[[k]]=lambda_k
          uniq_ks[k,]=Uniq
        }
        
        
        
        Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
        Phi_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # initialize group- and cluster-specific factor covariances
        alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
        for(g in 1:ngroup){
	        k=randpartvec[g]
          Psi_gs[[g]]=diag(uniq_ks[k,])
          X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
          for(k2 in 1:nclust){
	          lambda_k=Lambda_ks[[k2]]
	          Phi_gks[[g,k2]]=diag(nfactors)
            X_c=X_g-t(matrix(tau_ks[k2,],ncol=N_gs[g],nrow=nvar))
            udv=svd(X_c%*%lambda_k)
            U=udv$u
            V=udv$v
            Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
            Fvar=apply(Fscores,2,var)
            Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
            alpha_gks[[g,k2]] <- apply(Fscores,2,mean)
          }
        }

        
        Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          for(g in 1:ngroup){
            psi_g=Psi_gs[[g]]
            invPsi_g=diag(1/diag(psi_g))
            phi_gk=Phi_gks[[g,k]]
            invPhi_gk=phi_gk # phi_gk is still identity matrix
            sigma_gk=lambda_k %*% phi_gk %*% t(lambda_k) + psi_g
            Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
            invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_gk+t(lambda_k)%*%invPsi_g%*%lambda_k)
            invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
            invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(lambda_k)%*%invPsi_g; # Woodbury identity
          }
        }
        
        S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
        S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
        for(g in 1:ngroup){
          S_g=matrix(0,nvar,nvar)
          for(k in 1:nclust){
            lambda_k=Lambda_ks[[k]]
            # if(N_gks[g,k]!=0){
              X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
              Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar))
              S_gk=(1/N_gs[g])*(t(Xc_gk)%*%Xc_gk)
              S_gks[[g,k]] <- S_gk
              S_g=S_g+N_gks[g,k]*S_gk
            # } else {
            #   S_gks[[g,k]] <- matrix(0,nvar,nvar)
            # }
          }
          S_gs[[g]] <- (1/N_gs[g])*S_g
        }
        
        # compute Beta_gks and theta_gks
        Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
        #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
        for(g in 1:ngroup){
          for(k in 1:nclust){
            lambda_k=Lambda_ks[[k]]
            phi_gk=Phi_gks[[g,k]]
            invsigma_gk=invSigma_gks[[g,k]]
            beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
            Beta_gks[[g,k]]=beta_gk;
            S_gk=S_gks[[g,k]]
            theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_gk%*%t(beta_gk)
            Theta_gks[[g,k]]=theta_gk
            meanexpEta_gks[[g,k]]=(mu_gs[g,]-tau_ks[k,]-alpha_gks[[g,k]]%*%t(lambda_k))%*%t(beta_gk)
            #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
          }
        }
        
        Output_Mstep <- MixtureMG_FA_loadingsandintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mu_gs,tau_ks,alpha_gks)
        Lambda_ks=Output_Mstep$Lambda_ks
        Psi_gs=Output_Mstep$Psi_gs
        Phi_gks=Output_Mstep$Phi_gks
        tau_ks=Output_Mstep$tau_ks
        alpha_gks=Output_Mstep$alpha_gks
        Sigma_gks=Output_Mstep$Sigma_gks
        invSigma_gks=Output_Mstep$invSigma_gks
        nractivatedconstraints=Output_Mstep$nractivatedconstraints
        
        
        # compute observed-data log-likelihood for start
        ODLL_trialstart=0;
        loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
        loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
        for(g in 1:ngroup){
          X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
          for(k in 1:nclust){
            logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
            invSigma_gk=invSigma_gks[[g,k]]
            lambda_k=Lambda_ks[[k]]
            Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar)) # centered data per group
            loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
            for (n in 1:N_gs[g]){
              Xc_n=Xc_gk[n, ,drop=FALSE]
              loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_gk,Xc_n))
            }
            loglik_gks[g,k]=loglik_gk
            loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
          }
          m_i=max(loglik_gksw[g,]);
          for(k in 1:nclust){
            loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
          }
          ODLL_trialstart=ODLL_trialstart+log(sum(loglik_gksw[g,]))+m_i;
        }
        ODLLs_trialstarts[trialstart]=ODLL_trialstart;
      }
      sortedstarts <- sort(ODLLs_trialstarts,decreasing=TRUE,index.return=TRUE);
      index_beststarts=sortedstarts$ix[1:nruns]
      if (nrtrialstarts>nruns){
        randpartvecs=randpartvecs[index_beststarts,]
        if(nruns==1){
          randpartvecs=matrix(randpartvecs)
        }
      }
      else {
        nruns=nrtrialstarts
      }
    }
    else {
      nruns=1;
      randpartvecs=matrix(1,1,ngroup)
    }
  }
  
  # start of loop of multiple starts
  convergence <- 1
  logliks <- matrix(0,nruns,2)
  for(run in 1:nruns){
    nractivatedconstraints <- 0
    if(start==1){
      if(nruns>1){
        randpartvec <- randpartvecs[run,]
      }
      else {
        randpartvec <- randpartvecs
      }
    }
    if(start==2){
      randpartvec <- startpartition
    }
    
    if(nclust>1){
      z_gks=IM[randpartvec,]
      pi_ks=(1/ngroup)*apply(z_gks,2,sum)
    }
    else {
      z_gks=t(randpartvec)
      pi_ks=1
    }
    N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=apply(N_gks,2,sum)
    
    tau_ks <- matrix(0,nclust,nvar)
    for(k in 1:nclust){
      Xsup_k=matrix(0,N_ks[k],nvar)
      for(g in 1:ngroup){
        if(randpartvec[g]==k){
          X=Xsup[Ncum[g,1]:Ncum[g,2],]
          Xsup_k[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=X
        }
      }
      tau_ks[k,] <- apply(Xsup_k,2,mean)
    }
    # center data with initialized intercepts
    Xsupcent <- matrix(0,N,nvar)
    for(g in 1:ngroup){
      k=randpartvec[g]
      X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      Xsupcent[Ncum[g,1]:Ncum[g,2],]=X_g-t(matrix(tau_ks[k,],ncol=N_gs[g],nrow=nvar))
    }
    
    Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
    uniq_ks <- matrix(0,nclust,nvar)
    for(k in 1:nclust){
      Xsup_kc=matrix(0,N_ks[k],nvar)
      for(g in 1:ngroup){
        if(randpartvec[g]==k){
          Xc=Xsupcent[Ncum[g,1]:Ncum[g,2],]
          Xsup_kc[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=Xc
        }
      }
      S_k <- (1/N_ks[k])*(t(Xsup_kc)%*%Xsup_kc)
      ed<-eigen(S_k, symmetric=TRUE, only.values = FALSE)
      val<-ed$values
      u<-ed$vectors
      totalerror=sum((val[-seq_len(nfactors)]))
      meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
      Uniq=rep(meanerror,nvar)
      lambda_k=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
      if (EFA==0){
        lambda_k <- procr(lambda_k,design)
        lambda_k=lambda_k*design # non-zero loadings should be indicated with '1' for this to work properly
      }
      Lambda_ks[[k]]=lambda_k
      uniq_ks[k,]=Uniq
    }
    
    
    
    Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
    Phi_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # initialize group- and cluster-specific factor covariances
    alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
    for(g in 1:ngroup){
      k=randpartvec[g]
      Psi_gs[[g]]=diag(uniq_ks[k,])
      X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      for(k2 in 1:nclust){
        lambda_k=Lambda_ks[[k2]]
        Phi_gks[[g,k2]]=diag(nfactors)
        X_c=X_g-t(matrix(tau_ks[k2,],ncol=N_gs[g],nrow=nvar))
        udv=svd(X_c%*%lambda_k)
        U=udv$u
        V=udv$v
        Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
        Fvar=apply(Fscores,2,var)
        Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
        alpha_gks[[g,k2]] <- apply(Fscores,2,mean)
      }
    }
    
    # # initialize prior and posterior classification probabilities
    # if(nclust>1){
    #   z_gks=IM[randpartvec,]
    #   pi_ks=(1/ngroup)*apply(z_gks,2,sum)
    # }
    # else {
    #   z_gks=t(randpartvec)
    #   pi_ks=1
    # }
    
    
    Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      for(g in 1:ngroup){
        psi_g=Psi_gs[[g]]
        invPsi_g=diag(1/diag(psi_g))
        phi_gk=Phi_gks[[g,k]]
        invPhi_gk=solve(phi_gk)
        sigma_gk=lambda_k %*% phi_gk %*% t(lambda_k) + psi_g
        Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
        invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_gk+t(lambda_k)%*%invPsi_g%*%lambda_k)
        invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
        invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(lambda_k)%*%invPsi_g; # Woodbury identity
      }
    }
    
    # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
    for(g in 1:ngroup){
      X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        invSigma_gk=invSigma_gks[[g,k]]
        lambda_k=Lambda_ks[[k]]
        Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar)) # centered data per group
        loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
        for (n in 1:N_gs[g]){
          Xc_n=Xc_gk[n, ,drop=FALSE]
          loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_gk,Xc_n))
        }
        loglik_gks[g,k]=loglik_gk
      }
    }
    
    iter=0;
    conv1=1;
    conv2=1;
    ODLL=-Inf;
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      if (k==1){
        pars=lambda_k[design==1]
      }
      else {
        pars=c(pars,lambda_k[design==1])
      }
    }
    pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
    pars=unlist(pars)
    while(min(conv1,conv2)>1e-4 && iter<101){
      prev_ODLL=ODLL
      prev_Lambda_ks=Lambda_ks
      prev_Psi_gs=Psi_gs
      prev_Phi_gks=Phi_gks
      prev_tau_ks=tau_ks
      prev_alpha_gks=alpha_gks
      prev_pars=pars
      iter=iter+1
      
      # **E-step**: compute the posterior classification probabilities
      z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust, nfactors)
      
      
      N_gks=diag(N_gs[,1])%*%z_gks
      N_ks=apply(N_gks,2,sum)
      
      
      # update mixing proportions
      pi_ks=(1/ngroup)*apply(z_gks,2,sum)
      
      S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          # if(N_gks[g,k]!=0){
          X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
          Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar))
          S_gk=(1/N_gs[g])*(t(Xc_gk)%*%Xc_gk)
          S_gks[[g,k]] <- S_gk
          S_g=S_g+N_gks[g,k]*S_gk
          # } else {
          #   S_gks[[g,k]] <- matrix(0,nvar,nvar)
          # }
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }
      
      # compute Beta_gks and theta_gks
      Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
      #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
      for(g in 1:ngroup){
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          phi_gk=Phi_gks[[g,k]]
          invsigma_gk=invSigma_gks[[g,k]]
          beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
          Beta_gks[[g,k]]=beta_gk;
          S_gk=S_gks[[g,k]]
          theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_gk%*%t(beta_gk)
          Theta_gks[[g,k]]=theta_gk
          meanexpEta_gks[[g,k]]=(mu_gs[g,]-tau_ks[k,]-alpha_gks[[g,k]]%*%t(lambda_k))%*%t(beta_gk)
          #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
        }
      }
      
      Output_Mstep <- MixtureMG_FA_loadingsandintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mu_gs,tau_ks,alpha_gks)
      Lambda_ks=Output_Mstep$Lambda_ks
      Psi_gs=Output_Mstep$Psi_gs
      Phi_gks=Output_Mstep$Phi_gks
      tau_ks=Output_Mstep$tau_ks
      alpha_gks=Output_Mstep$alpha_gks
      Sigma_gks=Output_Mstep$Sigma_gks
      invSigma_gks=Output_Mstep$invSigma_gks
      nractivatedconstraints=Output_Mstep$nractivatedconstraints
      
      
      # check on change in observed-data log-likelihood
      ODLL=0;
      loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
      loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
      for(g in 1:ngroup){
        X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
        for(k in 1:nclust){
          logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
          invSigma_gk=invSigma_gks[[g,k]]
          lambda_k=Lambda_ks[[k]]
          Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar)) # centered data per group
          loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
          for (n in 1:N_gs[g]){
            Xc_n=Xc_gk[n, ,drop=FALSE]
            loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_gk,Xc_n))
          }
          loglik_gks[g,k]=loglik_gk
          loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
        }
        m_i=max(loglik_gksw[g,]);
        for(k in 1:nclust){
          loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
        }
        ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
      }
      
      
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        if (k==1){
          pars=lambda_k[design==1]
        }
        else {
          pars=c(pars,lambda_k[design==1])
        }
      }
      pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
      pars=unlist(pars)
      ind=which(pars!=0 & prev_pars!=0)
      parsdiff <- pars[ind]-prev_pars[ind] #mapply("-",pars,prev_pars)
      parsdiffdiv <- parsdiff/prev_pars[ind] #mapply("/",parsdiff,prev_pars)
      parsabsdiffdiv <- abs(parsdiffdiv)#lapply(parsdiffdiv,function(x){abs(x)})
      conv1 <- sum(parsabsdiffdiv)#do.call(sum,parsabsdiffdiv)
      
      conv2=ODLL-prev_ODLL
      
      #if(ODLL-prev_ODLL<0){
      #  ODLL-prev_ODLL
      #}
      
      
    } # end while-loop till convergence
    logliks[run,]=c(ODLL,nractivatedconstraints);
     
    
    if (run==1) {
      bestz_gks=z_gks
      bestpi_ks=pi_ks
      bestLambda_ks=Lambda_ks
      bestPsi_gs=Psi_gs
      bestPhi_gks=Phi_gks
      besttau_ks=tau_ks
      bestalpha_gks=alpha_gks
      bestSigma_gks=Sigma_gks
      bestinvSigma_gks=invSigma_gks
      bestloglik=ODLL
      bestloglik_gks=loglik_gks
      bestiter=iter
      bestconv1=conv1
      bestconv2=conv2
    }
    else {
      if (ODLL>bestloglik){
        bestz_gks=z_gks
        bestpi_ks=pi_ks
        bestLambda_ks=Lambda_ks
        bestPsi_gs=Psi_gs
        bestPhi_gks=Phi_gks
        besttau_ks=tau_ks
        bestalpha_gks=alpha_gks
        bestSigma_gks=Sigma_gks
        bestinvSigma_gks=invSigma_gks
        bestloglik=ODLL
        bestloglik_gks=loglik_gks
        bestiter=iter
        bestconv1=conv1
        bestconv2=conv2
      }
    }
    
    
  } # end for-loop over multiple starts
  
  z_gks=bestz_gks
  pi_ks=bestpi_ks
  Lambda_ks=bestLambda_ks
  Psi_gs=bestPsi_gs
  Phi_gks=bestPhi_gks
  tau_ks=besttau_ks
  alpha_gks=bestalpha_gks
  Sigma_gks=bestSigma_gks
  invSigma_gks=bestinvSigma_gks
  ODLL=bestloglik
  loglik_gks=bestloglik_gks
  iter=bestiter
  conv1=bestconv1
  conv2=bestconv2
  
  for(k in 1:nclust){
    lambda_k=Lambda_ks[[k]]
    if (k==1){
      pars=lambda_k[design==1]
    }
    else {
      pars=c(pars,lambda_k[design==1])
    }
  }
  pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
  pars=unlist(pars)
  while(min(conv1,conv2)>1e-6 && iter<Maxiter+1){ # iterate till convergence for best start
    prev_ODLL=ODLL
    prev_Lambda_ks=Lambda_ks
    prev_Psi_gs=Psi_gs
    prev_Phi_gks=Phi_gks
    prev_tau_ks=tau_ks
    prev_alpha_gks=alpha_gks
    prev_pars=pars
    iter=iter+1
    
    # **E-step**: compute the posterior classification probabilities
    z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust, nfactors)
    
    
    N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=apply(N_gks,2,sum)
    
    # update mixing proportions
    pi_ks=(1/ngroup)*apply(z_gks,2,sum)
    
    
    S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        # if(N_gks[g,k]!=0){
        X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
        Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar))
        S_gk=(1/N_gs[g])*(t(Xc_gk)%*%Xc_gk)
        S_gks[[g,k]] <- S_gk
        S_g=S_g+N_gks[g,k]*S_gk
        # } else {
        #   S_gks[[g,k]] <- matrix(0,nvar,nvar)
        # }
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }
    
    # compute Beta_gks and theta_gks
    Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
    #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
    for(g in 1:ngroup){
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        phi_gk=Phi_gks[[g,k]]
        invsigma_gk=invSigma_gks[[g,k]]
        beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
        Beta_gks[[g,k]]=beta_gk;
        S_gk=S_gks[[g,k]]
        theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_gk%*%t(beta_gk)
        Theta_gks[[g,k]]=theta_gk
        meanexpEta_gks[[g,k]]=(mu_gs[g,]-tau_ks[k,]-alpha_gks[[g,k]]%*%t(lambda_k))%*%t(beta_gk)
        #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
      }
    }
    
    Output_Mstep <- MixtureMG_FA_loadingsandintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mu_gs,tau_ks,alpha_gks)
    Lambda_ks=Output_Mstep$Lambda_ks
    Psi_gs=Output_Mstep$Psi_gs
    Phi_gks=Output_Mstep$Phi_gks
    tau_ks=Output_Mstep$tau_ks
    alpha_gks=Output_Mstep$alpha_gks
    Sigma_gks=Output_Mstep$Sigma_gks
    invSigma_gks=Output_Mstep$invSigma_gks
    nractivatedconstraints=Output_Mstep$nractivatedconstraints
    
    ODLL=0;
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
    loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
    for(g in 1:ngroup){
      X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        invSigma_gk=invSigma_gks[[g,k]]
        lambda_k=Lambda_ks[[k]]
        Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar)) # centered data per group
        loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
        for (n in 1:N_gs[g]){
          Xc_n=Xc_gk[n, ,drop=FALSE]
          loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_gk,Xc_n))
        }
        loglik_gks[g,k]=loglik_gk
        loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
      }
      m_i=max(loglik_gksw[g,]);
      for(k in 1:nclust){
        loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
      }
      ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
    }
    
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      if (k==1){
        pars=lambda_k[design==1]
      }
      else {
        pars=c(pars,lambda_k[design==1])
      }
    }
    pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
    pars=unlist(pars)
    ind=which(pars!=0 & prev_pars!=0)
    parsdiff <- pars[ind]-prev_pars[ind] #mapply("-",pars,prev_pars)
    parsdiffdiv <- parsdiff/prev_pars[ind] #mapply("/",parsdiff,prev_pars)
    parsabsdiffdiv <- abs(parsdiffdiv)#lapply(parsdiffdiv,function(x){abs(x)})
    conv1 <- sum(parsabsdiffdiv)#do.call(sum,parsabsdiffdiv)
    
    conv2=ODLL-prev_ODLL
    #if(ODLL-prev_ODLL<0){
    #  ODLL-prev_ODLL
    #}
    bestloglik=ODLL
    
    
  } # end while-loop till convergence
  
  if (conv2<1e-6){
    convergence=2 # convergence in terms of loglikelihood
  }
  else {
    if (conv1<1e-6) {
      convergence=1 # convergence in terms of parameters
    }
    else {
      convergence=0 # no convergence
    }
  }
  
  
  # set scale of factors across groups PER CLUSTER (and, in case of EFA, make them orthogonal)
  for(k in 1:nclust){
    if(N_ks[k]>0){
      theta_k=matrix(0,nfactors,nfactors)
      for(g in 1:ngroup){
        theta_gk=Theta_gks[[g,k]]
        theta_k=theta_k+(N_gks[g,k]/N_ks[k])*theta_gk;
      }
      if(nfactors>1){
        if(EFA==1){
          # find matrix square root via eigenvalue decomposition
          ed=eigen(theta_k)
          sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
          invsqrtFscale=solve(sqrtFscale)
        }
        else {
          sqrtFscale=diag(diag(theta_k^(1/2)))
          invsqrtFscale=diag(diag((1/theta_k)^(1/2)))
        }
      }
      else {
        sqrtFscale=sqrt(theta_k)
        invsqrtFscale=1/sqrtFscale
      }
      
      for(g in 1:ngroup){
        phi_gk=Phi_gks[[g,k]]
        phi_gk=invsqrtFscale%*%phi_gk%*%invsqrtFscale;
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry
      }
      Lambda_ks[[k]]=Lambda_ks[[k]]%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
      for(g in 1:ngroup){ # transform the factor means accordingly
        alpha_gks[[g,k]]=alpha_gks[[g,k]]%*%invsqrtFscale
      }
    }
  }
  
  # identification of factor means and intercepts
  mean_alpha_ks=matrix(0,nclust,nfactors)
  for(g in 1:ngroup){ 
    for(k in 1:nclust){
      mean_alpha_ks[k,]=mean_alpha_ks[k,]+N_gks[g,k]/N_ks[k]*alpha_gks[[g,k]]
    }
  }
  # Translation of factor means per cluster and cycle 1 update indicator intercepts (see De Roover, 2021)
  for(k in 1:nclust){
    for(g in 1:ngroup){
      alpha_gks[[g,k]]=alpha_gks[[g,k]]-mean_alpha_ks[k,]
    }
    if(N_ks[k]>0){
      lambda_k=Lambda_ks[[k]]
      suminvSigma=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvSigma=matrix(0,1,nvar)
      for(g in 1:ngroup){
        invSigma_gk=invSigma_gks[[g,k]]
        summeansminusalphaLambdainvSigma=summeansminusalphaLambdainvSigma+N_gks[g,k]*(mu_gs[g,]-alpha_gks[[g,k]]%*%t(lambda_k))%*%invSigma_gk
        suminvSigma=suminvSigma+N_gks[g,k]*invSigma_gk
      }
      tau_ks[k,]=t(solve(suminvSigma,t(summeansminusalphaLambdainvSigma)))
    }
  }
  
  
  
  if(EFA==1){
    nrpars=nclust-1+(nvar*nfactors-(nfactors*(nfactors-1)*(1/2)))*nclust+(nfactors*(nfactors+1)/2)*(ngroup-nclust)+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }
  else {
    nrpars=nclust-1+sum(design)*nclust+(nfactors*(nfactors+1)/2)*ngroup-nclust*nfactors+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }
  
 
  output_list <- list(z_gks=z_gks,pi_ks=pi_ks,Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,tau_ks=tau_ks,alpha_gks=alpha_gks,bestloglik=bestloglik,logliks=logliks,nrpars=nrpars,convergence=convergence,nractivatedconstraints=nractivatedconstraints)
  
  return(output_list)
} # end main function



# functions for E-step/posterior classification probabilities and M-step


# Update the cluster-membership probabilities z_gk
# Reuses the loglik_gks to save time

UpdPostProb <- function(pi_ks, loglik_gks, ngroup, nclust, nfact){
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow=ngroup,ncol=nclust)
  
  for(g in 1:ngroup){
    for(k in 1:nclust){
      z_gks[g,k] <- log(pi_ks[k])+loglik_gks[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow 
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclust))
  }
  
  # divide by the rowwise sum of the above calculated part 
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  z_gks <- round(z_gks,digits=16)
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks
  
  return(z_gks)
}



MixtureMG_FA_loadingsandintercepts_Mstep <- function(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mu_gs,tau_ks,alpha_gks){
  nractivatedconstraints <- 0
  ngroup <- length(N_gs)
  N_ks=apply(N_gks,2,sum)
  
  invPsi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
  for(g in 1:ngroup){
    psi_g=Psi_gs[[g]]
    invPsi_g=diag(1/diag(psi_g))
    invPsi_gs[[g]] <- invPsi_g
  }
  
  # update indicator intercepts
  for(k in 1:nclust){
    if(N_ks[k]>0){
      lambda_k=Lambda_ks[[k]]
      suminvPsi=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvPsi=matrix(0,1,nvar)
      for(g in 1:ngroup){
        invPsi_g=invPsi_gs[[g]]
        summeansminusalphaLambdainvPsi=summeansminusalphaLambdainvPsi+N_gks[g,k]*(mu_gs[g,]-alpha_gks[[g,k]]%*%t(lambda_k)-meanexpEta_gks[[g,k]]%*%t(lambda_k))%*%invPsi_g
        suminvPsi=suminvPsi+N_gks[g,k]*invPsi_g
      }
      tau_ks[k,]=t(solve(suminvPsi,t(summeansminusalphaLambdainvPsi)))
    }
  }
  
  # update factor means
  for(g in 1:ngroup){
    invPsi_g=invPsi_gs[[g]]
    for(k in 1:nclust){
      if(N_ks[k]>0){
        lambda_k=Lambda_ks[[k]]
        alpha_gks[[g,k]]=(mu_gs[g,]-tau_ks[k,]-meanexpEta_gks[[g,k]]%*%t(lambda_k))%*%t(solve(((t(lambda_k)%*%invPsi_g)%*%lambda_k),t(invPsi_g%*%lambda_k)))
      }
    }
  }
  
  # update factor loadings
  for(k in 1:nclust){
    if(N_ks[k]>0){
      lambda_k=matrix(0,nvar,nfactors)
      for(j in 1:nvar){
        nfactors_j=sum(design[j,])
        sumSbeta=matrix(0,1,nfactors_j)
        sumthetaalpha=matrix(0,nfactors_j,nfactors_j)
        summeansalpha=matrix(0,1,nfactors_j)
        for(g in 1:ngroup){
          psi_g=Psi_gs[[g]]
          S_gk=S_gks[[g,k]]
          beta_gk=Beta_gks[[g,k]]
          beta_gk=beta_gk[design[j,]==1, ,drop=FALSE]
          sumSbeta=sumSbeta+(N_gks[g,k]/psi_g[j,j])*S_gk[j,]%*%t(beta_gk)
          theta_gk=Theta_gks[[g,k]]
          theta_gk=theta_gk[design[j,]==1,design[j,]==1]
          alpha_gk=alpha_gks[[g,k]]
          alpha_gk=alpha_gk[design[j,]==1]
          meanexpeta_gk=meanexpEta_gks[[g,k]]
          meanexpeta_gk=meanexpeta_gk[design[j,]==1]
          sumthetaalpha=sumthetaalpha+(N_gks[g,k]/psi_g[j,j])*(theta_gk+alpha_gk%*%t(alpha_gk)+meanexpeta_gk%*%t(alpha_gk))
          summeansalpha=summeansalpha+(N_gks[g,k]/psi_g[j,j])*((mu_gs[g,j]-tau_ks[k,j])%*%alpha_gk)
        }
        lambda_k[j,design[j,]==1]= t(solve(sumthetaalpha,t(sumSbeta+summeansalpha)))
      }
      Lambda_ks[[k]]=lambda_k
    }
  }
  
  
  # update unique variances
  nractivatedconstraints=0
  for(g in 1:ngroup){
    S_g=S_gs[[g]]
    sum2SbetaB_BthetaB=0
    for(k in 1:nclust){
      if(N_ks[k]>0){
        lambda_k=matrix(0,nvar,nfactors)
        beta_gk=Beta_gks[[g,k]]
        theta_gk=Theta_gks[[g,k]]
        sum2SbetaB_BthetaB=sum2SbetaB_BthetaB+(N_gks[g,k]/N_gs[g])*(2*lambda_k%*%beta_gk%*%S_gk-lambda_k%*%theta_gk%*%t(lambda_k)) # modelimplied reduced covariance matrix on sample level, based on old structure matrix and sigma_gk, weighting based on new z_gks
      }
    }
    psi_g=diag(diag(S_g-sum2SbetaB_BthetaB))
    if (sum(diag(psi_g)<.0001)>0){ # track "heywood" cases
      ind=diag(psi_g)<.0001
      d=diag(psi_g);
      d[ind]=0.0001;
      psi_g=diag(d);
      nractivatedconstraints=nractivatedconstraints+sum(ind)
    }
    Psi_gs[[g]]=psi_g
  }
  
  # update factor (co)variances
  for(g in 1:ngroup){
    for(k in 1:nclust){
      if(N_ks[k]>0){
        theta_gk=Theta_gks[[g,k]]
        phi_gk=theta_gk
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry to avoid accumulation of asymmetry over iterations
      }
    }
  }
  
  
  # update (inv)Sigma_gks
  Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
  invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
  for(k in 1:nclust){
    lambda_k=Lambda_ks[[k]]
    for(g in 1:ngroup){
      psi_g=Psi_gs[[g]]
      invPsi_g=diag(1/diag(psi_g))
      phi_gk=Phi_gks[[g,k]]
      invPhi_gk=solve(phi_gk)
      sigma_gk=lambda_k %*% phi_gk %*% t(lambda_k) + psi_g
      Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
      invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_gk+t(lambda_k)%*%invPsi_g%*%lambda_k)
      invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
      invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(lambda_k)%*%invPsi_g; # Woodbury identity
    }
  }
  
  output_list <- list(Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,tau_ks=tau_ks,alpha_gks=alpha_gks,Sigma_gks=Sigma_gks,invSigma_gks=invSigma_gks,nractivatedconstraints=nractivatedconstraints)
  
  return(output_list)
}

# computation of adjusted rand index
adjrandindex <- function(part1,part2){
  
  IM1=diag(max(part1))
  IM2=diag(max(part2))
  A=IM1[part1,]
  B=IM2[part2,]
  
  T = t(A)%*%B
  N = sum(T)
  Tc = apply(T,2,sum)
  Tr = apply(T,1,sum)
  a = (sum(T^2) - N)/2
  b = (sum(Tr^2) - sum(T^2))/2
  c = (sum(Tc^2) - sum(T^2))/2
  d = (sum(T^2) + N^2 - sum(Tr^2) - sum(Tc^2))/2
  ARI = (choose(N,2)*(a + d) - ((a+b)*(a+c)+(c+d)*(b+d)))/(choose(N,2)^2 - ((a+b)*(a+c)+(c+d)*(b+d)))
  
  return(ARI)
}

# Procrustes rotation (orthogonal)
procr <- function(x,y){
  s <- svd(t(x)%*%y)
  U <- s$u # X = U D V'
  D <- s$d
  V <- s$v
  R <- U%*%t(V)
  yhat <- x%*%R # rotated x that approximates y
  
  return(yhat)
}


