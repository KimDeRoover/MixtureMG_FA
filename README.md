# MixtureMG_FA: R-code for mixture multigroup factor analysis

For now, please cite preprint: https://psyarxiv.com/7fdwv

This version deals with factor loading differences (finds clusters of groups based on similarity of their factor loadings, given the user-specified number of clusters) and uses EFA

For model selection, it is advised to use BIC_G (number of groups as sample size) in combination with CHull (see preprint, https://www.rdocumentation.org/packages/multichull/versions/1.0.0)

# INPUT:
Xsup = data matrix for all groups (rows are subjects nested within groups, columns are the variables to be factor-analyzed)

N_gs = vector with number of subjects for each group (in the same order as they appear in the data matrix)

nclust = user-specified number of clusters

nfactors = user-specified number of factors

Maxiter = maximum number of iterations

nruns = number of starts (based on pre-selected random partitions, when start = 1)

startpartition = partition of groups to start from (use with start = 2 and nruns = 1)

# OUTPUT:
z_gks = cluster memberships of groups (posterior classification probabilities)

pi_ks= mixing proportions (prior classification probabilities)

Lambda_ks = cluster-specific loadings, access loadings of cluster k via Lambda_ks[[k]], rotate as desired

Psi_gs = group-specific unique variances, access loadings of group g via Psi_gs[[g]]

Phi_gks = group- and cluster-specific factor (co)variances, access (co)variances of group g in cluster k via Phi_gks[[g,k]]

mu_gs = group-specific means, access means of group g via mu_gs[g,]

bestloglik = loglikelihood of best start (iterated till complete convergence)

logliks = loglikelihoods of all starts (iterated till preliminary convergence, see paper)

nrpars = number of free parameters, to be used for model selection in combination with bestloglik

convergence = 2 if converged on loglikelihood, 1 if converged on parameter changes, 0 if not converged code for mixture multigroup factor analysis

If you labeled the output of MixtureMG_FA as 'Output', extract the different parameters as follows: e.g., Lambda_ks=Output$Lambda_ks.
