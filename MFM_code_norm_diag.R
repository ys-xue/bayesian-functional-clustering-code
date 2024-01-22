# rm(list=ls())
## We consider undirected graph with selfloops

cut = function(d,c)
{as.numeric(d<=c)}

## Function for log-likelihood related to Jth observation
loglike <- function(clusterassign,param_mu,param_tau, data,J,n) #here J means Jth observation
{
  ################################################################
  
  ## Input: clusterassign = clustering configuration, a n by 1 vector ##
  ##        param_mu, param_tau = probability matrix, a k by k matrix ##
  ##        data = the adjacency matrix, a n by n matrix ##
  ##        J = observation index ##
  ##        n = number of observations ##
  
  ## Output: log-likelihood related to Jth observation ##
  
  #################################################################
  clustersize = max(clusterassign)
  param_mu = as.matrix(param_mu)
  param_tau = as.matrix(param_tau)
  if (J==1) {result2 = 0
  for (ii in c((J):n))
  {
    result2 = result2 + dnorm(data[J,ii],mean = param_mu[clusterassign[J],clusterassign[ii]],sd = sqrt(1/param_tau[clusterassign[J],clusterassign[ii]]),log = T)
    #data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
  }
  output = sum(result2)} else if (J==n){
    result = 0
    for (ii in c(1:(J)))
    {
      result = result + dnorm(data[ii,J],mean = param_mu[clusterassign[ii],clusterassign[J]],sd = sqrt(1/param_tau[clusterassign[ii],clusterassign[J]]),log = T)
      #data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    output = sum(result)
  } else {
    result = 0
    for (ii in c(1:(J)))
    {
      result = result + dnorm(data[ii,J],mean = param_mu[clusterassign[ii],clusterassign[J]],sd = sqrt(1/param_tau[clusterassign[ii],clusterassign[J]]),log = T)
    }
    
    result2 = 0
    for (ii in c((J+1):n))
      
    {
      result2 = result2 + dnorm(data[J,ii],mean = param_mu[clusterassign[J],clusterassign[ii]],sd = sqrt(1/param_tau[clusterassign[J],clusterassign[ii]]),log = T)
    }
    output = sum(result)+sum(result2)}
  output
}

#function for getting m(Aj)
logmargs <- function(clusterassign,data,J,mu_0, t_0, alpha, beta) #here J means Jth observation
{
  ################################################################
  
  ## Input: clusterassign = clustering configuration, a n by 1 vector ##
  ##        data = the adjacency matrix, a n by n matrix ##
  ##        J = observation index ##
  ##        n = number of observations ##
  ##        beta.a, beta.b = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  
  ## Output: m(A_j) in Algotithm 1 (collapsed sampler for MFM-SBM) ##
  
  #################################################################
  clustersize = max(clusterassign)-1
  result = NULL
  for (ii in 1:clustersize)
  {
    if (length(which(clusterassign==ii))==0) {result[ii] = 0} else{
      meanA = mean(data[J,which(clusterassign==ii)])
      S = length(which(clusterassign==ii))
      alphan = alpha + S/2
      tn = t_0 + S
      betan = beta + (meanA-mu_0)^2*t_0*S/(2*(S+t_0)) + 1/2*ifelse(is.na(var(data[J,which(clusterassign==ii)])*(S-1)),0,var(data[J,which(clusterassign==ii)])*(S-1))
      
      result[ii] = lgamma(alphan)-lgamma(alpha) + alpha*log(beta) - alphan*log(betan) + 0.5*(log(t_0)-log(tn)) - S/2*log(2*pi)
      
      # sumA =  sum(data[J,which(clusterassign==ii)[which(clusterassign==ii)>J]]) + sum(data[which(clusterassign==ii)[which(clusterassign==ii)<J],J])
      # S = length(which(clusterassign==ii)[which(clusterassign==ii)>J]) + length(which(clusterassign==ii)[which(clusterassign==ii)<J])
      # result[ii] = lbeta(sumA+beta.a,S-sumA+beta.b)-lbeta(beta.a,beta.b)
    }}
  sum(result)
}


## function for Collapsed sampler for MFM-SBM (main algorithm)
CDMFM <- function(data, data1, lambda1, neighbour,distance, niterations, mu_0,mu_0off, t_0, alpha, beta, GAMMA, LAMBDA, initNClusters)
{
  ## Model: A_{ij}|z,Q \sim Normal(mu_{z_i,z_j},tau_{z_i,z_j}) ##
  ##        (\mu_{r,s},\sigma_{r,s}) \sim NIG with (0,1,1)
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  
  
  ################################################################
  
  ## Input: data = the connection matrix, a n by n matrix ##
  ##        data1 = the upper traiangle for the adjacency matrix, a n by n matrix ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        beta.a, beta.b = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distrition ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         Qout = probability matrix, a k by k matrix ##
  
  #################################################################
  n = dim(data)[1]
  #precomputation for prespecified coefficient VN
  lambda <- LAMBDA
  gamma <- GAMMA
  N=n ## n is the number of oberservations
  VN<-0
  tmax = n+10
  for (t in 1:tmax)
  {
    r = log(0)
    for (k in t:500)
    {
      b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
      m = max(b,r)
      r = log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] = r
  }
  # initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  tau<-matrix(0, initNClusters,initNClusters)
  mu<-matrix(0, initNClusters,initNClusters)
  for (i in 1:initNClusters){
    for (j in i:initNClusters){
      tau[i,j] <- rgamma(1,alpha, beta);
      tau[j,i] = tau[i,j]
      if(i == j) {
        mu[i,j] <- rnorm(1,mu_0,sqrt(1/(t_0*tau)));
        mu[j,i] = mu[i,j]} else {
          mu[i,j] <- rnorm(1,mu_0off,sqrt(1/(t_0*tau)))
          mu[j,i] = mu[i,j];
        }
    }
  }
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          cost = exp(lambda1*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*exp(loglike(clusterAssign_temp,mu,tau,data,i,n))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(logmargs(clusterAssign_1,data,i,mu_0off, t_0, alpha, beta))*exp(VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          tau1 = matrix(0,nClusters+1,nClusters+1)
          tau1[1:nClusters,1:nClusters] = tau
          tau1[nClusters+1,1:(nClusters+1)] = rgamma(nClusters+1,alpha, beta)
          tau1[1:(nClusters+1),nClusters+1] = tau1[nClusters+1,1:(nClusters+1)]
          tau = tau1
          mu1 = matrix(0,nClusters+1,nClusters+1)
          mu1[1:nClusters,1:nClusters] = mu
          mu1[nClusters+1,1:(nClusters+1)] = c(rnorm(nClusters,mu_0off,sqrt(1/(t_0*tau1[nClusters+1,1:(nClusters+1)]))),rnorm(1,mu_0off,sqrt(1/(t_0*tau1[nClusters+1,1:(nClusters+1)])))) #rgamma(nClusters+1,alpha, beta)
          mu1[1:(nClusters+1),nClusters+1] = mu1[nClusters+1,1:(nClusters+1)]
          mu = mu1
          
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {tau = tau
          mu = mu
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          cost = exp(lambda1*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*exp(loglike(clusterAssign_temp,mu,tau,data,i,n))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(logmargs(clusterAssign_1,data,i,mu_0off, t_0, alpha, beta))*exp(VN[nClusters]-VN[nClusters-1])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else
        {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          if (nClusters > 1) {mu = mu[-cur.cluster.i,][,-cur.cluster.i];tau = tau[-cur.cluster.i,][,-cur.cluster.i]} else {mu = mu[-cur.cluster.i,][-cur.cluster.i];tau = tau[-cur.cluster.i,][-cur.cluster.i]}}
      }
    }
    # end for loop over subjects i
    ## update Q ##
    tau = matrix(0, nClusters,nClusters)
    mu = matrix(0, nClusters,nClusters)
    AA = matrix(0,nClusters,nClusters)
    NN = matrix(0,nClusters,nClusters)
    for (r in 1:nClusters){
      for (s in r:nClusters)
      {
        # med = matrix(0,n,n)
        # med[which(clusterAssign==r),which(clusterAssign==s)] = 1
        # med1 = matrix(0,n,n)
        # med1[which(clusterAssign==s),which(clusterAssign==r)] = 1
        # nr = sum(med*lower.tri(med)) + sum(med1*lower.tri(med1))-(r==s)*sum(med1*lower.tri(med1))
        # sumh = sum(data1[clusterAssign==r,clusterAssign==s]) + sum(data1[clusterAssign==s,clusterAssign==r]) - (r==s)*sum(data1[clusterAssign==s,clusterAssign==r])
        if (r != s) {datapoint = as.vector(data[clusterAssign==r,clusterAssign==s])} else {
          datapoint = as.vector(data[clusterAssign==r,clusterAssign==s]*lower.tri(data[clusterAssign==r,clusterAssign==s],diag = T))
          datapoint = datapoint[datapoint!=0]
        }
        
        if (length(datapoint)==0){
          tau[r,s] = rgamma(1,shape=alpha,rate= beta)
          mu[r,s] = rnorm(1,mean = mu_0,sd = sqrt(1/(t_0*tau[r,s])))
        } else {
          nr = length(datapoint)
          sumh = sum(datapoint)
          thetah = sumh/nr
          kr = 1/(t_0+nr)
          tau[r,s] = rgamma(1,shape=(alpha+nr/2),rate= beta + 1/2*(nr*t_0/(t_0+nr)*(thetah-mu_0)^2+ifelse(is.na(var(datapoint)*(nr-1)),0,var(datapoint)*(nr-1)))/2)
          if(r==s) {
            mu[r,s] = rnorm(1,mean = kr*nr*thetah + t_0*kr*mu_0,sd = sqrt(kr*1/tau[r,s]))} else {
              mu[r,s] = rnorm(1,mean = kr*nr*thetah + t_0*kr*mu_0off,sd = sqrt(kr*1/tau[r,s]))
            }
        }
        
        
        mu[s,r] = mu[r,s]
        tau[s,r] = tau[r,s]
        
      }
    }
    History[[iter]] <- list(zout = clusterAssign,tauout = tau,muout = mu)
    # cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

## Dahl's method to summarize the samples from the MCMC
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}
# 
# distance <- readRDS("USgdist.rds")
# 
# n=51
# A = log((inner_prod_matrix-min(inner_prod_matrix))/(max(inner_prod_matrix)-inner_prod_matrix))
# A[which(A==Inf)] = max(A[which(is.finite(A)==T)]) 
# A[which(A==-Inf)] = min(A[which(is.finite(A)==T)]) 
# #diag(A) = 0
# #A = matrix(scale(as.vector(A)),n,n)
# AAA = matrix(0,n,n) ##the upper traiangle for the adjacency matrix
# for (i in 1:n){
#   for (j in i:n){
#     #A[j,i] = A[i,j]
#     AAA[i,j] = A[i,j]
#   }
# }
# # diag(AAA) = 0 ## make it without-selfloop network
# # diag(A) = 0 ## make it without-selfloop networkk
# 
# ## taking the data into the MFM-SBM algorithm
# set.seed(4)
# fit1 = CDMFM_new(data = A, data1 = AAA, lambda1 = 2, neighbour = 4,distance = distance, niterations = 1000, mu_0 = 5, mu_0off = -5, t_0 = 2, alpha = 1, beta = 1, GAMMA=1, LAMBDA = 1, initNClusters = 9)
# 
# result1 = getDahl(fit1, 500)
# 
# save(result1, file = "result_2_4.RData")



calcDB <- function(srvf, final_cluster) {
  #########
  ## inputs:
  ## srvf = the matrix of srvf's for all 51 states
  ## final_cluster = the inference result, result$zout
  K <- unique(final_cluster)
  DB <- 0
  meansrvf <- matrix(0, nrow = nrow(srvf), ncol = length(K))
  for (i in seq_along(K)) {
    meansrvf[, i] <- rowMeans(srvf[, which(final_cluster == K[i]), drop = FALSE])
  }
  ## a matrix version of the term inside Equation (14)
  ratioMatrix <- matrix(0, length(K), length(K)) 
  for (i in 1:length(K)) {
    for (j in 1:length(K)) {
      if (i != j) {
        ratioMatrix[i, j] <- (calcL2(srvf[, which(final_cluster == K[i]),
                                          drop = FALSE]) + 
                                calcL2(srvf[, which(final_cluster == K[j]),
                                            drop = FALSE])) / 
          sqrt(sum( (meansrvf[,i] - meansrvf[,j])^2 ))
      }
    }
  }
  # return(ratioMatrix)
  DB <- mean(sum(matrixStats::rowMaxs(ratioMatrix)))
  return(DB)
}

calcSSE <- function(srvf, final_cluster) {
  K <- unique(final_cluster)
  tosum <- 0
  meansrvf <- matrix(0, nrow = nrow(srvf), ncol = length(K))
  for (i in seq_along(K)) {
    meansrvf[, i] <- rowMeans(srvf[, which(final_cluster == K[i]), drop = FALSE])
  }
  for (i in seq_along(K)) {
    tosum <- tosum +  sum(sweep(srvf[, which(final_cluster == K[i]), drop = FALSE], 1, meansrvf[, i]) ^ 2)
  }
  return(tosum)
}

## calculate pairwise difference
calcL2 <- function(srvf) {
  coll <- ncol(srvf)
  outmat <- matrix(0, nrow = coll, ncol = coll)
  for (i in 1:coll) {
    for (j in i:coll) {
      outmat[i, j] <- outmat[j, i] <- sqrt(sum((srvf[,i] - srvf[,j])^2))
    }
  }
  if (coll == 1) {
    return(0)
  } else {
    ## exclude distances from itself
    return(sum(outmat) / (coll - 1)^2)
  }
  
}


calcRand <- function(fit, ind) {
  ## history = an item from fitList
  # RIs <- matrix(0, nrow = length(history), ncol = length(history[[1]]$Iterates))
  # for (i in 1:nrow(RIs)) {
  #   for (j in 1:ncol(RIs)) {
  #     RIs[i, j] <- fossil::rand.index(history[[i]]$Iterates[[j]]$zout, ind)
  #   }
  # }
  # return(RIs)
  ## history = resultLust[[1]]$fit1$Iterates
  
  RIs <- map_dbl(fit$Iterates, ~fossil::rand.index(.x$zout, ind))
  return(RIs)
}
