########################################
## Mixture Multigroup Factor Analysis ##
########################################
# Main function for performing MMG-FA with selection of the number of clusters
# Code written by Kim De Roover
# Please cite publications: https://doi.apa.org/doi/10.1037/met0000355
# and https://www.tandfonline.com/doi/full/10.1080/10705511.2020.1866577
# See these papers for recommendations on sample size, number of starts, model selection and how to move on from the results of these analyses.
# The version that finds clusters of groups based on similarity of their factor loadings and intercepts at the same time will be added soon

# INPUT:
# data = data file you want to analyze: FIRST COLUMN SHOULD CONTAIN NUMERICAL GROUP IDs, all remaining variables are factor-analyzed
# cluster.spec = either "loadings" or "intercepts" to find clusters of groups based on equivalence of factor loadings or intercepts (imposing equal loadings), respectively,
# or c("loadings","intercepts") to find clusters of groups based on equivalence of loadings AND intercepts
# nsclust = vector of length two, indicating the minimal and maximal number of clusters (it is recommended to set the minimal number to one)
# nfactors = number of factors
# Maxiter = maximum number of iterations used in each MMG-FA analysis
# nruns = number of (preselected) random starts (important setting to avoid local maxima in case of few groups and/or small groups)
# design = matrix indicating position of zero loadings with '0' and non-zero loadings with '1' (for CFA, leave unspecified for EFA)
#          (using different design matrices for different clusters is not supported)
# rotation = rotation criterion, currently either "oblimin" or "varimax" (0 = no rotation)
# preselect = percentage of best starts taken in pre-selection (increase to speed up startprocedure) 

# OUTPUT:
# overview = overview of fitted MMG-FA solutions with loglikelihood (loglik), number of parameters (nrpars), BIC_G (i.e., using number of groups as sample size), 
# convergence (1 = converged) and number of activated constraints on the unique variances
# MMGFA_solutions = list of all MMG-FA solutions with different numbers of clusters
#          Access parameter values of solution with preferred number of clusters as, for example, OutputObject$MMGFA_solutions$"2.clusters"


MixtureMG_FA <- function(data,cluster.spec,nsclust,nfactors,Maxiter = 5000,nruns = 50,design=0,rotation=0,preselect = 10){
  if(rotation!=0){
    library("GPArotation")
  }
  nvar <- ncol(data)-1
  Xsup=data[,2:(nvar+1)]
  groups=unique(data[,1])
  N_gs<-tabulate(data[,1])
  N_gs=N_gs[N_gs!=0]
  N=sum(N_gs)
  G=length(N_gs)
  
  overview=matrix(0,nsclust[2]-nsclust[1]+1,6)
  MMGFA_solutions=matrix(list(NA),nrow = 1, ncol=nsclust[2]-nsclust[1]+1)
  if (length(cluster.spec)==1){
    if (cluster.spec=="loadings"){
      for(nclust in nsclust[1]:nsclust[2]){
        if(nclust==1){
          print(paste("Fitting MMG-FA with",nclust,"cluster ..."),quote=FALSE)
          output_nclust <- MixtureMG_FA_loadings(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else {
          print(paste("Fitting MMG-FA with",nclust,"clusters ..."),quote=FALSE)
          if (nclust==G){
            output_nclust <- MixtureMG_FA_loadings(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
          }
          else {
            output_nclust <- MixtureMG_FA_loadings(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
          }
        }
        loglik=output_nclust$bestloglik
        nrpars=output_nclust$nrpars
        convergence=output_nclust$convergence>0
        overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
        if(design==0 && rotation!=0){
          for(k in 1:nclust){
            if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
              rot<-GPForth(output_nclust$Lambda_ks[[k]], method="varimax")
            }
            if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
              rot<-GPFoblq(output_nclust$Lambda_ks[[k]], method="oblimin")
            }
            rotatedloadings=rot$loadings
            output_nclust$Lambda_ks[[k]]=rotatedloadings
            T_matrix=rot$Th
            invT=solve(T_matrix)
            for(g in 1:G){ # counter-rotate all corresponding sets of factor (co)variances
              output_nclust$Phi_gks[[g,k]]=invT%*%output_nclust$Phi_gks[[g,k]]%*%t(invT)
            }
          }
        }
        prefix="Cluster"
        suffix=seq(1:nclust)
        names(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
        colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
        output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,groupspecific.uniquevariances=output_nclust$Psi_gs,groupspecific.means=output_nclust$mu_gs)
        
        MMGFA_solutions[[nclust-nsclust[1]+1]]=output_nclust2
      }
    }
    if (cluster.spec=="intercepts"){
      for(nclust in nsclust[1]:nsclust[2]){
        if(nclust==1){
          print(paste("Fitting MMG-FA with",nclust,"cluster ..."),quote=FALSE)
          output_nclust <- MixtureMG_FA_intercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else {
          print(paste("Fitting MMG-FA with",nclust,"clusters ..."),quote=FALSE)
          if(nclust==G){
            output_nclust <- MixtureMG_FA_intercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
          }
          else{
            output_nclust <- MixtureMG_FA_intercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
          }
        }
        loglik=output_nclust$bestloglik
        nrpars=output_nclust$nrpars
        convergence=output_nclust$convergence>0
        overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
        if(design==0 && rotation!=0){
          if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
            rot<-GPForth(output_nclust$Lambda, method="varimax")
          }
          if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
            rot<-GPFoblq(output_nclust$Lambda, method="oblimin")
          }
          rotatedloadings=rot$loadings
          output_nclust$Lambda=rotatedloadings
          T_matrix=rot$Th
          invT=solve(T_matrix)
          for(g in 1:G){ # counter-rotate all sets of factor (co)variances
            output_nclust$Phi_gs[[g]]=invT%*%output_nclust$Phi_gs[[g]]%*%t(invT)
            for(k in 1:nclust){ # counter-rotate all sets of factor means
              output_nclust$alpha_gks[[g,k]]=output_nclust$alpha_gks[[g,k]]%*%t(invT)
            }
          }
        }
        prefix="Cluster"
        suffix=seq(1:nclust)
        rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
        colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
        colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
        output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,invariant.loadings=output_nclust$Lambda,groupspecific.factorcovariances=output_nclust$Phi_gs,groupspecific.uniquevariances=output_nclust$Psi_gs,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)
        
        MMGFA_solutions[[nclust-nsclust[1]+1]]=output_nclust2
      }
    }
  } else {
    if (is.element("loadings",cluster.spec) & is.element("intercepts",cluster.spec)){
      for(nclust in nsclust[1]:nsclust[2]){
        if(nclust==1){
          print(paste("Fitting MMG-FA with",nclust,"cluster ..."),quote=FALSE)
          output_nclust <- MixtureMG_FA_loadingsandintercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else {
          print(paste("Fitting MMG-FA with",nclust,"clusters ..."),quote=FALSE)
          if(nclust==G){
            output_nclust <- MixtureMG_FA_loadingsandintercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
          }
          else{
            output_nclust <- MixtureMG_FA_loadingsandintercepts(Xsup,N_gs,nclust,nfactors,Maxiter = Maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
          }
        }
        loglik=output_nclust$bestloglik
        nrpars=output_nclust$nrpars
        convergence=output_nclust$convergence>0
        overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
        if(design==0 && rotation!=0){
          for(k in 1:nclust){
            if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
              rot<-GPForth(output_nclust$Lambda_ks[[k]], method="varimax")
            }
            if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
              rot<-GPFoblq(output_nclust$Lambda_ks[[k]], method="oblimin")
            }
            rotatedloadings=rot$loadings
            output_nclust$Lambda_ks[[k]]=rotatedloadings
            T_matrix=rot$Th
            invT=solve(T_matrix)
            for(g in 1:G){ # counter-rotate all sets of factor (co)variances and factor means
              output_nclust$Phi_gks[[g,k]]=invT%*%output_nclust$Phi_gks[[g,k]]%*%t(invT)
              output_nclust$alpha_gks[[g,k]]=output_nclust$alpha_gks[[g,k]]%*%t(invT)
            }
          }
        }
        prefix="Cluster"
        suffix=seq(1:nclust)
        rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
        colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
        colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
        output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,groupspecific.uniquevariances=output_nclust$Psi_gs,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)
        
        MMGFA_solutions[[nclust-nsclust[1]+1]]=output_nclust2
      }
    } 
  }
  if (is.nan(max(overview))==FALSE){
    screeratios=matrix(0,nsclust[2]-nsclust[1]+1,1)
    for(nclust in (nsclust[1]+1):(nsclust[2]-1)){
      LL_nclust=overview[nclust-nsclust[1]+1,2]
      npar_nclust=overview[nclust-nsclust[1]+1,3]
      LL_nclustmin1=overview[nclust-nsclust[1],2]
      npar_nclustmin1=overview[nclust-nsclust[1],3]
      LL_nclustplus1=overview[nclust-nsclust[1]+2,2]
      npar_nclustplus1=overview[nclust-nsclust[1]+2,3]
      screeratios[nclust-nsclust[1]+1]=((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
    }
    overview=cbind(overview[,1:4],screeratios,overview[,5:6])
    prefix=seq(nsclust[1]:nsclust[2])
    suffix="clusters"
    sollistnames<-c(paste(prefix,suffix,sep="."))
    names(MMGFA_solutions)<-noquote(sollistnames)
    overview=as.data.frame(overview,row.names=FALSE)
    colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_G","screeratios","convergence","nr.activated.constraints")
    par(mfrow=c(2,1))
    par(mar=c(4.1,4.1,2.1,2.1))
    plot(overview[,1],overview[,4],type="b",xlab="number of clusters",ylab="BIC_G")
    plot(overview[,1],overview[,2],type="b",xlab="number of clusters",ylab="loglik")
    
    cat("\n")
    print(overview)
    cat("\n")
    cat(">> Choose the best number of clusters ('K_best') based on the plots and access the corresponding cluster memberships and parameter estimates by using, for example, OutputObject$clustermemberships[[K_best]] and OutputObject$clusterspecific.loadings[[K_best]] or OutputObject$clusterspecific.intercepts[[K_best]].")
    cat("\n")
    cat(">> The parameter sets are further subdivided in group- and/or cluster-specific parameter sets.")
    output_MS<-list(overview=overview,MMGFA_solutions=MMGFA_solutions)
  } else {
    cat("You seem to have made a typo in 'cluster.spec' (make sure you use lowercase letters) or have requested a model that is not (or not yet) supported.")
    cat("\n")
    cat("Please try again.")
    cat("\n")
  }
}
