#' tloho_lm: Tree-based Low rank Horseshoe regularization prior in linear model
#' 
#' tloho_lm: Tree-based Low rank Horseshoe regularization prior in linear model
#' 
#' The main R function of the paper 
#' 
#' T-LoHo: A Bayesian Regularization Model for Structured Sparsity and Smoothness on Graphs
#' by C. J. Lee, Z. T. Luo, and H. Sang, NeurIPS 2021
#'
#'
#' @param Y n by 1 scalar, real-valued response variables. 
#' @param X n by p matrix, real-valued predictors. If not standardized, it will be standardized so that all columns have L2 norm 1
#' @param graph0 igraph object, reflecting the structure of beta. 
#' @param init_val list with name 'beta' and/or 'trees', Initial value of beta and/or spanning forest. 
#' @param c model size penalization hyperparameter, value between 0 and 1. c = 1   
#' @param tau0 shrinkage strength hyperparameter, positive real.
#' @param MCMC number of total MCMC iteration.
#' @param BURNIN number of burn-in iteration which will be discarded
#' @param THIN thin-in rate. The final number of posterior sample is nsamples = (MCMC-BURNIN)/THIN.
#' @param loss "binder" or "VI", whether to calculate cluster estimate based on Binder loss (default) or VI loss
#' @param hsplus default F, use horseshoe+ prior instead of horseshoe? (experimental)
#' @param seed seed.
#'
#' @import igraph salso mgcv gsl
#'
#' @return  A list containing:\tabular{ll}{
#'    \code{beta_out} \tab nsamples by p matrix, posterior samples of beta\cr
#'    \tab \cr
#'    \code{lambda2_out} \tab list of length nsamples, posterior samples of lambda^2 \cr
#'    \tab \cr
#'    \code{tau2_out} \tab length nsamples vector, posterior samples of tau^2 \cr
#'    \tab \cr
#'    \code{sigmasq_y_out} \tab length nsamples vector, posterior samples of sigma^2 \cr
#'    \tab \cr
#'    \code{cluster_out} \tab nsamples by p matrix, posterior samples of cluster label \cr
#'    \tab \cr
#'    \code{cluster_map} \tab length p vector, maximum a posteriori cluster estimate \cr
#'    \tab \cr
#'    \code{cluster_est} \tab length p vector, cluster estimate based on Bayes estimator that minimizes expected loss \cr
#'    \tab \cr
#'    \code{log_post_out} \tab length nsamples vector, log posterior likelihood(up to a constant) \cr
#'    \tab \cr
#'    \code{acc} \tab list with split/merge/change/hyper move proposal counts and acceptance counts \cr
#'    \tab \cr
#'    \code{map_beta_est} \tab length p vector, maximum a posteriori beta estimate\cr
#'    \tab \cr
#'    \code{mean_beta_est} \tab length p vector, posterior mean beta estimate \cr
#'    \tab \cr
#'    \code{median_beta_est} \tab length p vector, posterior median beta estimate \cr
#'    \tab \cr
#'    \code{map_MSF_est} \tab igraph object, a sample of minimum spanning forest when posterior likelihood is maximized\cr
#'    \tab \cr
#' }
#' @export
#' 
#' 
tloho_lm <- function(Y, X, graph0, init_val=NULL, c = 0.5, tau0 = 1, MCMC = 50000, BURNIN = 40000, THIN = 10, loss = "binder", hsplus = F, seed=NULL){
  ## sanity check ----
  set.seed(seed)
  p = ncol(X) # = vcount(graph0)
  n = nrow(X)
  X = unname(X)
  if(length(Y)!=n) stop("length of Y is not match with nrow of X")
  if(length(V(graph0))!=p) stop("number of verticies in a graph does not match with ncol of X")
  
  if(any(abs(apply(X, 2, function(x) sum(x^2))-1) > 1e-15)){
    cat("design matrix X will be standardized using scale(X)/sqrt(n-1) so that all columns have L2 norm 1.\n")
    X = scale(X)/sqrt(n-1)
  }
  if(!(c >= 0 && c < 1)) stop("c should be between 0(inclusive) and 1(exclusive)")
  if(!(tau0 > 0)) stop("tau0 should be positive")
  
  if(loss == "binder"){
    cat("point estimate will based on binder distance \n")
  }else if(loss == "VI"){
    cat("point estimate will based on Variation of information distance\n")
  }else{
    stop("provide 'binder' or 'VI' for the input of 'loss' argument \n")
  }
  
  if(all(X[lower.tri(X)] == 0, X[upper.tri(X)] == 0, n==p) ){
    cat("X = I after column standardization, i.e. it is normal means model...\n")
    normalmeans = T
  }else{
    normalmeans = F
  }
  
  # sanity check for graph0
  if(any(clusters(graph0)$csize==1)) cat(paste("note: graph contains",sum(clusters(graph0)$csize==1),"isolated nodes in the graph\n"))
  inc_mat = get.edgelist(graph0, names = F) 
  graph0 = graph_from_edgelist(inc_mat, directed = F)
  
  vertex_idx = 1:p 
  E(graph0)$eid = c(1:ecount(graph0))  # edge id
  V(graph0)$vid = c(1:vcount(graph0))  # vertex id

  Comp = clusters(graph0)$no
  move_cnt <- move_acc <- rep(0,4)

  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = delete_vertex_attr(graph0, 'name')
  }
  cat("generating adjacent edge lists from given graph... ")
  adj_list = unname(lapply(as_adj_list(graph0), function(x) x$vid )) 
  adj_edge_list = unname(lapply(as_adj_edge_list(graph0), function(x) x$eid )) 
  cat(" complete! \n")

  if(!is.null(init_val[['beta']])){
    est.Beta <- init_val[['beta']]
    if(length(est.Beta)!=p) stop("invalid initial beta. should be a vector with length vcount(graph)")
    cluster = c(factor(est.Beta, labels = c(1:length(unique(est.Beta)))))
    cluster = as.integer(cluster)
    beta_tilde = unique(est.Beta) # contiguity check needed - see below
  }else{
    cat("Initial beta not provided. Using Clauset-Newman-Moore graph clustering method to get initial clusters of beta...\n")
    temp = cluster_fast_greedy(graph0)
    cluster = membership(temp)
    k = length(sizes(temp))
    beta_tilde = rnorm(k)
  }

  if(!is.null(init_val[['trees']])){
    mstgraph = init_val[['trees']]
  }else{
    edge_status = getEdgeStatus(cluster, inc_mat)
    mstgraph = proposeMST(graph0, edge_status)
  }

  # update eid_btw_mst
  inc_mat_mst = get.edgelist(mstgraph, names = F)
  idx_btw = which(cluster[inc_mat_mst[, 1]] != cluster[inc_mat_mst[, 2]])
  eid_btw_mst = (E(mstgraph)$eid)[idx_btw]
  
  k = max(cluster)
  clust_vid = unname(split(1:p, cluster)) #clust_vid[[k]] is vid of kth cluster

  Phi = getPhi(clust_vid, p)
  
  beta = t(Phi)%*%beta_tilde
  Xtilde = tcrossprod(X, Phi) 
  Xtilde = as.matrix(Xtilde)
  #We avoid Phi calculation during MCMC, we directly update Xtilde = X*t(Phi), not Phi
  
  csize = as.numeric(table(cluster))# cluster size
 
  #contiguity check
  for(iii in 1:k){ 
    mstsubgraph <- induced_subgraph(mstgraph,clust_vid[[iii]])
    if((vcount(mstsubgraph)-1)!=ecount(mstsubgraph)){
      betaduplicate = beta_tilde[iii]
      cat(paste("initial cluster (induced by beta) is not contiguous. check beta value",betaduplicate,"\n"," and vertex id",clust_vid[[iii]],"\n"))
      stop("stop")
    }
  }
  cat("initial cluster contiguity check complete\n")
  
  #### initialize ####
  
  a0 <- b0 <- 0; #cat("set sigma^2 hyperparameter to be a0 = b0 = 0(Jeffreys)\n");
  cat(paste0("set tau0 = ",tau0,"\n"));
  cat(paste0("set c = ",c,"\n"));
  sigMH_tau = 1; cat("set M-H proposal sd sigMH_tau = 1 \n")
  hyper = c(a0, b0, tau0, c)
  tau02 = tau0^2
  lambda2 = rep(1, length(beta_tilde)) #local shrinkage parameter. initialize with 1
  tau2 = tau02 # global shrinkage parameter
  PRECISION = 1/(lambda2*tau2) # length K precision vector

  yty = sum(Y ^ 2) # t(Y) %*% Y
  # cholesky factor of lambda*I + t(Xtilde) %*% Xtilde
  XtX = crossprod(Xtilde)
  diag(XtX) = diag(XtX) + PRECISION # initial local parameter is 1
  R = chol(XtX)
  
  # t(Xtilde) %*% Y
  Xty = as.numeric(crossprod(Xtilde, Y))

  log_like_res = evalmLogLike_lambdavec(R, Xty, PRECISION, a0, b0, yty, n)
  log_like = log_like_res$mlog_like
  bsol = log_like_res$bsol  # solution to t(R) %*% z = t(Xtilde) %*% Y

  # whether an edge in graph0 is within a cluster or bewteen two clusters
  edge_status = apply(matrix(cluster), 2, FUN = getEdgeStatus, inc_mat)

  # MCMC saving objects
  sigmasq_y_out = numeric((MCMC-BURNIN)/THIN)
  lambda2_out = list()
  tau2_out = numeric((MCMC-BURNIN)/THIN)

  cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, p))
  MST_out = list()
  log_post_out = matrix(0,(MCMC-BURNIN)/THIN)
  beta_est = array(0, dim = c((MCMC-BURNIN)/THIN, p))
  
  tau_acc = numeric(MCMC)
  
  # RJMCMC starts ----------------------------
  
  pd = 0.05 # hyper move probability
  for(iter in 1:MCMC) {
      
      ## Step 1 ----------------
      
      if(k == Comp) {pa = 0.95; pb = 0; pc = 0 # only birth/hyper step allowed
      } else if(k == p) {pa = 0; pb = 0.95; pc = 0 # death/change/hyper allowed
      } else {pa = 0.425; pb = 0.425; pc = 0.1}
      move = sample(4, 1, prob = c(pa, pb, pc, pd))
      move_cnt[move] = move_cnt[move] + 1
     
      #cat(iter, move, log_like,'\n', sep = ' ')
      if(move == 1) { 
        ### Step 1(a), split -------
        
        # split an existing cluster - (modified) old cluster size is always greater or equal to the new cluster
        split_res = splitCluster.loho(mstgraph, k, clust_vid, csize)
        vid_new = split_res$vid_new; vid_old = split_res$vid_old
        clust_old = split_res$clust_old;
        clust_new = k + 1
        #PRECISION inherits to the larger splitted cluster
        # prior proposal
        if(!hsplus){
          lambda_star = abs(rcauchy(1)) 
        }else{
          lambda_star = abs(rcauchy(1, scale = abs(rcauchy(1))))
        }
        PRECISION_star <- 1/(lambda_star^2*tau2)
        PRECISION_starvec = PRECISION
        PRECISION_starvec <- append(PRECISION_starvec, PRECISION_star, after = clust_new-1)

        # compute log-prior ratio - prior proposal
        log_A = log(1-c) #+ log(k+1 - Comp) - log(n-k) + 2*dcauchy(lambda_new2,log=T)
        # compute log-proposal ratio
        if(k == p-1) {
          pb_new = 0.85
        } else {pb_new = 0.425}

        log_P = log(pb_new) - log(pa)  #- log(k+1 - Comp)+ log(n-k) - 2*dcauchy(lambda_new2, log = T)

        # update Cholesky factor
        oldcol = rowSums(X[,vid_old, drop = F])/sqrt(length(vid_old))
        newcol = rowSums(X[,vid_new, drop = F])/sqrt(length(vid_new)) # size of vid_new is less than vid_old
        
        if(normalmeans){
          R_new = diag(sqrt(1+PRECISION_starvec), nrow = k+1) # for normal means
        }else{
          R_new = cholSplit.loho(R, Xtilde, lambda_old = PRECISION[clust_old], PRECISION_star,
                                 oldcol, newcol, clust_old = clust_old)
        }
        #
        # update Xty
        Xty_new = updateXty.loho('split', Xty, Y, idx_old = clust_old, idx_new = clust_new,
                                oldcol = oldcol, newcol = newcol, csize = NULL)
        
        # compute log-likelihood ratio
        log_like_res = evalmLogLike_lambdavec(R_new, Xty_new, PRECISION_starvec, a0, b0, yty, n)
        log_like_new = log_like_res$mlog_like
        log_L = log_like_new - log_like

        #acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          move_acc[1] = move_acc[1] + 1

          Xtilde = updateXtilde.loho('split', Xtilde, clust_old, clust_new, oldcol, newcol, csize = NULL)
          update_res = updateSplit(split_res, clust_vid, k, csize, eid_btw_mst,
                                   cluster, edge_status, adj_list, adj_edge_list)
          clust_vid = update_res$clust_vid
          csize = update_res$csize
          eid_btw_mst = update_res$eid_btw_mst
          cluster = update_res$cluster
          k = k + 1

          log_like = log_like_new; bsol = log_like_res$bsol
          R = R_new; Xty = Xty_new
          PRECISION = PRECISION_starvec
          lambda2 = 1/(PRECISION*tau2)
          edge_status = update_res$estatus
          }
      }

      if(move == 2) { 
        ### Step 1(b), merge-------------------
        # merge two existing clusters
        merge_res = mergeCluster.loho(mstgraph, eid_btw_mst, clust_vid, csize,
                                          cluster, inc_mat)
        vid_old = merge_res$vid_old
        clust_old = merge_res$clust_old; clust_new = merge_res$clust_new

        PRECISION_starvec <- PRECISION[-(clust_old)] # erase clust_old th PRECISION
        #if(merge_res$oldinherit){
        #  PRECISION_star <- PRECISION[clust_old] # merging cluster precision - inherits old
        #}else{
          PRECISION_star <- PRECISION[clust_new] #merging cluster precision
        #}
        PRECISION_starvec[clust_new] <- PRECISION_star

        # compute log-prior ratio
        log_A = -log(1-c) #- log(k - Comp) + log(n-k+1) - 2*dcauchy(lambda_old,log=T)
        # # compute log-proposal ratio
        if(k == 1 + Comp) {pa_new = 0.95
        }else {pa_new = 0.425}
        log_P = -(log(pb) - log(pa_new))  #+ log(k - Comp) - log(n-k+1) + 2*dcauchy(lambda_old, log = T)

        #browser()
        if(normalmeans){
          R_new = diag(sqrt(1+PRECISION_starvec), nrow = k-1) # for normal means
        }else{
          R_new = cholMerge.loho(R, Xtilde, PRECISION_star, clust_old, clust_new, csize)#, #merge_res$oldinherit)
        }
        # browser()
        Xty_new = updateXty.loho('merge', Xty, Y, idx_old = clust_old, idx_new = clust_new,
                               oldcol = NULL, newcol = NULL, csize = csize)
        
        # compute log-likelihood ratio
        log_like_res = evalmLogLike_lambdavec(R_new, Xty_new, PRECISION_starvec, a0, b0, yty, n)
        log_like_new = log_like_res$mlog_like
        log_L = log_like_new - log_like

        #acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          move_acc[2] = move_acc[2] + 1
          Xtilde = updateXtilde.loho('merge', Xtilde, clust_old, clust_new, oldcol=NULL, newcol=NULL, csize)
          
          update_res = updateMerge(merge_res, clust_vid, csize, eid_btw_mst,
                                   cluster, edge_status, adj_list, adj_edge_list)
          clust_vid = update_res$clust_vid

          csize = update_res$csize
          eid_btw_mst = update_res$eid_btw_mst
          cluster = update_res$cluster
          k = k - 1

          #XtX = crossprod(Xtilde)
          log_like = log_like_new; bsol = log_like_res$bsol
          R = R_new; Xty = Xty_new
          PRECISION = PRECISION_starvec
          lambda2 = 1/(PRECISION*tau2)
          edge_status = update_res$estatus
        }
      }
  
      if(move == 3){ 
        ### step 1(c), change ------
        
        # first perform death move:
        merge_res = mergeCluster.loho(mstgraph, eid_btw_mst, clust_vid, csize,
                                          cluster, inc_mat, change = T)
        vid_old = merge_res$vid_old
        clust_old = merge_res$clust_old; clust_new = merge_res$clust_new

        vid_new = merge_res$vid_new
        
        PRECISION_starvec <- PRECISION[-(clust_old)] # erase clust_old th PRECISION
        # if(merge_res$oldinherit){
        #   PRECISION_star <- PRECISION[clust_old] # merging cluster precision - inherits old
        # }else{
          PRECISION_star <- PRECISION[clust_new] # merging cluster precision
        #}
        PRECISION_starvec[clust_new] <- PRECISION_star 

        # update Cholesky factor and Xty
        if(normalmeans){
          R_new = diag(sqrt(1+PRECISION_starvec), nrow = k-1) # for normal means
        }else{
          R_new = cholMerge.loho(R, Xtilde, PRECISION_star, clust_old, clust_new, csize)#, #merge_res$oldinherit)
        }
        Xty_new = updateXty.loho('merge', Xty, Y, idx_old = clust_old, idx_new = clust_new,
                               oldcol = NULL, newcol = NULL, csize = csize)
        # update Xtilde
        Xtilde_new = updateXtilde.loho('merge', Xtilde, clust_old, clust_new, oldcol=NULL, newcol=NULL, csize)
        
        # then perform birth move
        split_res = splitCluster.loho(mstgraph, k-1, merge_res$clust_vid, merge_res$csize)

        vid_new = split_res$vid_new; vid_old = split_res$vid_old
        clust_old = split_res$clust_old; clust_new = k
        
        if(!hsplus){
          lambda_starstar = abs(rcauchy(1)) 
        }else{
          lambda_starstar = abs(rcauchy(1, scale = abs(rcauchy(1))))
        }
        PRECISION_starstar <- 1/(lambda_starstar^2*tau2) 
        PRECISION_starstarvec <- append(PRECISION_starvec, PRECISION_starstar, after = clust_new-1)
        
        # update Cholesky factor
        oldcol = rowSums(X[,vid_old, drop = F])/sqrt(length(vid_old))
        newcol = rowSums(X[,vid_new, drop = F])/sqrt(length(vid_new)) # size of vid_new is less than vid_old
        
        R_new =  cholSplit.loho(R_new, Xtilde_new, lambda_old = PRECISION_starvec[clust_old], PRECISION_starstar,
                              oldcol, newcol, clust_old = clust_old)
        
        Xty_new = updateXty.loho('split', Xty_new, Y, idx_old = clust_old, idx_new = clust_new,
                               oldcol = oldcol, newcol = newcol, csize = NULL)
        log_A = 0 #- 2*dcauchy(eta_min,log=T) +  2*dcauchy(eta_new_new2,log=T)
        log_P = 0 #2*dnorm(eta_min, log=T) - 2*dnorm(eta_new_new2, log=T)

        # compute log-likelihood ratio
        log_like_res = evalmLogLike_lambdavec(R_new, Xty_new, PRECISION_starstarvec, a0, b0, yty, n)
        log_like_new = log_like_res$mlog_like
        log_L = log_like_new - log_like

        # acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          move_acc[3] = move_acc[3] + 1
          Xtilde = updateXtilde.loho('split', Xtilde_new, clust_old, clust_new, oldcol, newcol, csize = NULL)
          
          update_res = updateMerge(merge_res, clust_vid, csize, eid_btw_mst,
                                   cluster, edge_status, adj_list, adj_edge_list)
          update_res = updateSplit(split_res, update_res$clust_vid, k-1, update_res$csize,
                                   update_res$eid_btw_mst, update_res$cluster,
                                   update_res$estatus, adj_list, adj_edge_list)
          clust_vid = update_res$clust_vid
          csize = update_res$csize
          eid_btw_mst = update_res$eid_btw_mst
          cluster = update_res$cluster

          log_like = log_like_new; bsol = log_like_res$bsol
          R = R_new; Xty = Xty_new
          PRECISION = PRECISION_starstarvec
          lambda2 = 1/(PRECISION*tau2)
          edge_status = update_res$estatus
        }
      }

      if(move == 4) { 
        ### Step 1(d), hyper ####
        # update MST
        move_acc[4] = move_acc[4] + 1
        mstgraph = proposeMST(graph0, edge_status)
        # update eid_btw_mst
        inc_mat_mst = get.edgelist(mstgraph, names = F)
        idx_btw = which(cluster[inc_mat_mst[, 1]] != cluster[inc_mat_mst[, 2]])
        eid_btw_mst = (E(mstgraph)$eid)[idx_btw]
      }
      
      # Step 2 -----------------------------
      ## Step 2-1, update tau ------------
      #update tau via MH
      #propose new tau2 using lognormal proposal
      tau_new = exp(rnorm(1, 0.5*log(tau2), sigMH_tau))
      tau2_new = tau_new^2
      PRECISION_newtau = 1/(lambda2*tau2_new)
      
      R_new <- cholesky_diagonalcpp(R, diagadd = PRECISION_newtau - PRECISION)
      
      log_like_res = evalmLogLike_lambdavec(R_new, Xty, PRECISION_newtau, a0, b0, yty, n)
      log_like_new = log_like_res$mlog_like
      log_L = log_like_new - log_like

      log_A = log(1+tau2/tau02) - log(1+tau2_new/tau02) # cauchy prior
      
      # acceptance probability
      acc_prob = min(0, log_A + log_L + log(tau_new) - 0.5*log(tau2))
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        tau2 = tau2_new
        PRECISION = PRECISION_newtau
        log_like = log_like_new; bsol = log_like_res$bsol
        R = R_new
        tau_acc[iter] = 1
      }
      # adaptive M-H to maintain acceptance ratio close to 0.35
      if (iter %% 1000 == 0) {
        if(mean(tau_acc[(iter - 999) : iter]) >= 0.35) {
          scale.adj <- exp(min(0.01, 1 / sqrt(iter / 1000)))
        } else if (mean(tau_acc[(iter - 999) : iter]) < 0.35) {
          scale.adj <- exp(-min(0.01, 1 / sqrt(iter / 1000)))
        # } else if (mean(tau_acc[(iter - 99) : iter]) >= 0) { # too low acceptance rate. stuck, escape by increasing scale.adj
        #   scale.adj <- exp(min(0.1, 1 / sqrt(iter / 100)))
        } else {
          scale.adj <- 1
        }
        sigMH_tau <- sigMH_tau * scale.adj
      }

     ## Step 2-2, update sigmasq ---------------
     YPY = yty - drop(sum(bsol^2))  # Y'P(lambda)^{-1}Y
     sigmasq_y = 1/rgamma(1, shape = (n+a0)/2, rate = (b0+YPY)/2)
     
     ## Step 2-3, update tilde(beta) ----------
     beta_tilde = drop(backsolve(R/sqrt(sigmasq_y), bsol/sqrt(sigmasq_y) + rnorm(k)))# rcMVN(R/sqrt(sigmasq_y), drop(bsol)/sqrt(sigmasq_y))

     # Step 3 --------------------------
     # code reference: R package horseshoe
     # update local shrinkage parameters: update lambda_j's in a block using slice sampling ##
     #browser()
     eta = 1/(lambda2)
     if(!hsplus){
       upsi = stats::runif(k,0,1/(1+eta))
     }else{
       upperbound1 = log(eta)/(eta-1)
       upperbound1[which(is.nan(upperbound1))] = 1 # when eta = 1
       upsi = stats::runif(k,0, upperbound1)
     }
     tempps = beta_tilde^2/(2*sigmasq_y*tau2)
     if(!hsplus){
       ub = (1-upsi)/upsi
     }else{
       ub = -gsl::lambert_Wm1(-upsi*exp(-upsi))/upsi #not W0
     }
     # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
     Fub = 1 - exp(-tempps*ub) # exp cdf at ub
     Fub[Fub < (1e-8)] = 1e-8;  # for numerical stability
     up = stats::runif(k,0,Fub)
     eta = -log(1-up)/tempps
     lambda2 = 1/eta;
     
     PRECISION_new = 1/(lambda2*tau2) 
     R_new <- cholesky_diagonalcpp(R, diagadd = PRECISION_new - PRECISION)
     
     log_like_res = evalmLogLike_lambdavec(R_new, Xty, PRECISION_new, a0, b0, yty, n)
     log_like = log_like_res$mlog_like
     bsol = log_like_res$bsol
     R = R_new
     PRECISION = PRECISION_new
     
     # sort cluster indices in the decreasing size order
   
     
     
    # save result ----------------
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      idx = (iter-BURNIN)/THIN
      
      beta_est[idx, ] = as.vector(crossprod(getPhi(clust_vid, p),beta_tilde))
      lambda2_out[[idx]] = lambda2
      tau2_out[idx] = tau2
      sigmasq_y_out[idx] = sigmasq_y
      MST_out[[idx]] = mstgraph
      cluster_out[idx, ] = cluster

      log_post_out[idx] = evalLogPost_HS(beta_tilde, sigmasq_y, lambda2, tau2, k, Y, Xtilde, hyper, Comp)
    }
     
    if(iter %% 1000 == 0) {
      cat('Iteration', iter, 'done, ncluster:', k,', log_likelihood(up to const, beta/sigma integrated out):',log_like, '\n')
    }

  }# end mcmc
  
  # find MAP
  map_idx = which.max(log_post_out)
  cluster_map = cluster_out[map_idx,]
  map_beta_est = beta_est[map_idx,]
  
  # find posterior mean
  mean_beta_est = apply(beta_est, 2, FUN = mean)
  median_beta_est = apply(beta_est, 2, FUN = median)
  
  cat("move count(split/merge/change/hyper)=", move_cnt,"\n")
  cat("move accepted(split/merge/change/hyper)=", move_acc,"\n")
  
  mode(cluster_out) <- "integer"
  mode(cluster_map) <- "integer"
  mode(tau_acc) <- "integer"
  
  cat("get Bayes estimator of partition (clusters) that minimizes loss ")
  cluster_est = salso::dlso(cluster_out, loss = loss)
  
  return(list('beta_out' = beta_est,
              'lambda2_out' = lambda2_out,
              'tau2_out' = tau2_out,
              'sigmasq_y_out' = sigmasq_y_out, 
              'cluster_out' = cluster_out,
              'cluster_map' = cluster_map,
              'cluster_est' = cluster_est,
              'log_post_out' = log_post_out,
              'acc' = list(move_cnt = move_cnt, move_acc = move_acc),
              'map_beta_est' = map_beta_est, 
              'mean_beta_est' = mean_beta_est,
              'median_beta_est' = median_beta_est,
              'map_MSF_est' = MST_out[[map_idx]]))
}