# Specifing dependency packages... now included in imports
# require(igraph)
# require(mgcv) # for mgcv::choldrop and mgcv::cholup
# require(salso) # for cluster point estimate using Dahl's method


# Auxiliary functions for final likelihood calculation  --------------

# function to update log likelihood - marginlized out sigmasq_y

#' Helper function for evaluating log-likelihood
#' 
#' This function gives two output: 
#' 1. mlog_like, -0.5log(|Sigma|)- n/2log(t(y)Sigma^{-1}y/2) 
#' 2. bsol, (t(R))^{-1}*t(Xtilde)y
#'
#' @param R R matrix, K by K right triangular matrix, chol(tau^{-2}Lambda^{-1} + crossprod(Xtilde)). 
#' @param Xty t(Xtilde)*y
#' @param PRECISION tau^{-2}Lambda^{-1}
#' @param a0 (default 0) first hyperparameter for sigma2
#' @param b0 (default 0) second hyperparameter for sigma2
#' @param yty crossprod(y) value.
#' @param n nobs.
#' @param bs_CVAC (default NULL) bsol if precalculated 
#' 
#' @return
#'
#' @keywords internal
#' @noRd
evalmLogLike_lambdavec <- function(R, Xty, PRECISION, a0 = 0, b0 = 0, yty, n, bs_CVAC = NULL) {
  logdet = (1/2) * sum(log(PRECISION))
  logdet = logdet - sum(log(diag(R)))
  if(is.null(bs_CVAC)) {bs_CVAC = forwardsolve(t(R), Xty)}
  distval_CVAC = sum(bs_CVAC^2)
  mlog_like = logdet -(n+a0)/2*log(b0/2 + yty/2 - distval_CVAC/2)#+ distval_CVAC / (2*sigmasq_y)
  return(list(mlog_like = as.numeric(mlog_like), bsol = bs_CVAC))
}

# function to get log posterior density (up to a constant)
evalLogPost_HS <- function(beta_all, sigmasq_y, lambda2, tau2, k, Y, Xtrans, hyper, Comp) {
  lambda = sqrt(lambda2)
  tau = sqrt(tau2)
  n = length(Y); p = length(k)
  a0 = hyper[1]; b0 = hyper[2]; tau0 = hyper[3]; c = hyper[4]
  log_prior =  -(a0/2+1)*log(sigmasq_y) - b0/(2*sigmasq_y) # prior for gamma
  #log_prior = log_prior + (c0/2-1)*log(lambda) - d0*lambda/2
  log_prior = log_prior + sum(dcauchy(lambda, log = T)) + dcauchy(tau, scale = tau0, log = T) ####### prior for lambda and tau
  
  log_prior = log_prior + sum(-lchoose(n-Comp, k-Comp) + k*log(1-c))
  
  #log_prior = log_prior - sum(k)/2*(log(sigmasq_y)-log(lambda)) - lambda/(2*sigmasq_y)*sum(beta_all^2)
  log_prior = log_prior + sum(dnorm(beta_all, mean = 0, sd = lambda*tau*sqrt(sigmasq_y), log =T))
  log_like = -n/2*log(sigmasq_y) - sum((Y - Xtrans%*%beta_all)^2) / (2*sigmasq_y)
  return(log_like + log_prior)
}





# Auxiliary functions for cholesky rank-1 update --------------

# function for updating R when removing j-th column and row of the original matrix t(R)%*%R
# input: n by n matrix R, index j(1 to n)
# output: (n-1) by (n-1) matrix R_new
cholDel <- function(R, j){
  n = nrow(R)
  if(j == n){
    R_new = R[-n, -n, drop = F]
  }else{
    R_new = mgcv::choldrop(R, j)
  }
  R_new
}

# function for updating R when adding [A1, A2, A3] as j-th column and row of the original matrix t(R)%*%R
# input: n by n matrix R , index j(1 to n+1), 
# vector [A1 A2 A3] where A1 has length (j-1), A2 has length 1, A3 has length n+1-j
# output: (n+1) by (n+1) matrix R_new
cholAdd <- function(R, j, A2, A1 = NULL, A3 = NULL) {
  n = nrow(R)
  if(j == n+1) {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(n > 0){
      R_new[1:n,1:n] = R
      S12 = drop(backsolve(R, A1, transpose = T))
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[1:(n+1), n+1] = c(S12, S22)
    }else{
      R_new[1,1] = A2
    }
  } else {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(j > 1) {
      R11 = R[1:(j-1), 1:(j-1), drop = F]
      R_new[1:(j-1), 1:(j-1)] = R11
      S12 = backsolve(R11, A1, transpose = T)
      R_new[1:(j-1), j] = S12
      S13 = R[1:(j-1), j:n, drop = F]
      R_new[1:(j-1), (j+1):(n+1)] = S13
      
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[j, j] = S22
      S23 = (t(A3) - crossprod(S12, S13)) / S22
      R_new[j, (j+1):(n+1)] = S23
      S33 = mgcv::cholup(R[j:n, j:n, drop = F], S23, up = FALSE) # downdate
      R_new[(j+1):(n+1), (j+1):(n+1)] = S33
    } else {
      S22 = sqrt(as.numeric(A2))
      R_new[1, 1] = S22
      S23 = as.numeric(A3) / S22
      R_new[1, 2:(n+1)] = S23
      S33 = mgcv::cholup(R, S23, up = FALSE) # downdate
      R_new[2:(n+1), 2:(n+1)] = S33
    }
  }
  return(R_new)
}



# function for updating R when original matrix is t(R)%*%R =  diag(PRECISION) + t(Xtilde) %*% Xtilde
# input: k by k matrix R, Xtilde,
# lambda_old and lambda_new are scalar values corresponds to old precision and new precision
# oldcol and newcol corresponds to the columns of Xtilde
# clust_old: replacing index. new column corresponding to newcol is added at the end
# output: (k+1) by (k+1) matrix R_new
cholSplit.loho <- function(R, Xtilde, lambda_old, lambda_new, oldcol, newcol, clust_old){
  K = ncol(Xtilde)
  # remove a column
  idx_rm = clust_old
  R_new = cholDel(R, idx_rm)
  # replace with oldcol column
  col_new = oldcol
  A1 = NULL; A3 = NULL
  if(idx_rm > 1) {A1 = crossprod(Xtilde[, 1:(idx_rm-1), drop=F], col_new)}
  if(idx_rm < K) {A3 = crossprod(Xtilde[, (idx_rm+1):K, drop=F], col_new)}
  if((idx_rm == 1)&&(K==1)){
    R_new = sqrt(crossprod(col_new) + lambda_old)  
  }else{
    R_new = cholAdd(R_new, idx_rm, A2 = sum(col_new^2)+lambda_old, A1 = A1, A3 = A3)
  }
  # add a new column 
  idx_add = K + 1 # last position
  col_new = newcol
  
  A3 = NULL # last column so no A_3
  Xtilde[,idx_rm] <- oldcol # this is second adding step.
  A1 = crossprod(Xtilde, col_new)
  #if(change ==T) A3 = crossprod(Xtilde[, idx_add:k_sum, drop=F], col_new)
  R_new = cholAdd(R_new, idx_add, A2 = sum(col_new^2)+lambda_new, A1 = A1, A3 = A3)
  
  return(R_new)
}



# function for updating R when original matrix is t(R)%*%R =  diag(PRECISION) + t(Xtilde) %*% Xtilde
# input: k by k matrix R, Xtilde,
# lambda_new is scalar values corresponds to merging cluster
# clust_old: removing index.
# clust_new: merging index. 
# csize: cluster size, needed for normalization
# output: (k-1) by (k-1) matrix R_new
cholMerge.loho <- function(R, Xtilde, lambda_new, clust_old, clust_new, csize){
  K = ncol(Xtilde)
  
  idx_add = clust_new + 1
  idx_rm_old = clust_old
  col_new = (Xtilde[,clust_old]*sqrt(csize[clust_old]) + Xtilde[,clust_new]*sqrt(csize[clust_new]))/sqrt(csize[clust_old]+csize[clust_new])  # new column for Xtilde
  A1 = NULL; A3 = NULL
  if(idx_add > 1) {A1 = crossprod(Xtilde[, 1:(idx_add-1), drop=F], col_new)} # A_12 matrix
  if(idx_add < K+1) {A3 = crossprod(Xtilde[, idx_add:(K), drop=F], col_new)} # t(A_23) matrix
  
  R_new = cholAdd(R, idx_add, A2 = sum(col_new^2)+lambda_new, A1 = A1, A3 = A3)
  
  # remove a column for clust_new
  idx_rm_new = idx_add - 1
  R_new = cholDel(R_new, idx_rm_new)
  
  # remove a column for clust_old
  R_new = cholDel(R_new, idx_rm_old)
  
  return(R_new)
}

# Auxiliary functions during proposals --------------------


# function to propose a new MST
proposeMST <- function(graph0, edge_status) {
  nedge = length(edge_status) # edge_status is TRUE if crossing(between cluster), FALSE if within cluster
  nb = sum(edge_status)
  nw = nedge - nb
  weight = numeric(nedge)
  weight[!edge_status] = runif(nw, 0, 1/2)
  weight[edge_status] = runif(nb, 1/2, 1)
  mstgraph = mst(graph0, weights = weight)
  return(mstgraph)
}

# function to update Xty
updateXty.loho <- function(move, Xty, y, idx_old, idx_new, oldcol=NULL, newcol = NULL, csize=NULL){
  K = length(Xty)
  if(move == 'split') {
    if(is.null(newcol)) stop("newcol is empty")
    if(is.null(oldcol)) stop("oldcol is empty")
    val_old = sum(oldcol * y)
    val_new = sum(newcol * y)
    Xty[idx_old] = val_old
    if(idx_new <= K) {
      Xty = c(Xty[1:(idx_new-1)], val_new, Xty[idx_new:K])
    } else {
      Xty = c(Xty,val_new)
    }
  }
  if(move == 'merge') {
    if(is.null(csize)) stop("csize is empty")
    Xty[idx_new] = (Xty[idx_old]*sqrt(csize[idx_old]) + Xty[idx_new]*sqrt(csize[idx_new]))/sqrt(csize[idx_old]+csize[idx_new])
    Xty = Xty[-idx_old]
  }
  return(as.numeric(Xty))
}


# membership: cluster_id -> vid
# cluster: vid -> cluster_id
# have a constraint that new clster size is always less than old cluster size
splitCluster.loho<- function(mstgraph, k, membership, csize) {
  clust_split = sample.int(k, 1, prob = csize - 1)
  mst_subgraph = induced_subgraph(mstgraph, membership[[clust_split]])
  
  edge_cutted = E(mst_subgraph)[sample.int(csize[clust_split]-1, 1)]
  eid_cutted = edge_cutted$eid
  mst_subgraph = delete.edges(mst_subgraph, edge_cutted)
  connect_comp = components(mst_subgraph)
  clust_vid_new = split(V(mst_subgraph)$vid, connect_comp$membership)
  # vid_new = (V(mst_subgraph)$vid)[cluster_new == 2]  # vid for vertices belonging to new cluster
  # membership[vid_new] = k + 1
  
  #introduce the constraint; |new cluster| <= |old cluster| where |.| is size
  if(length(clust_vid_new[[2]])>length(clust_vid_new[[1]])){
    dummy <- clust_vid_new[[1]]
    clust_vid_new[[1]] <- clust_vid_new[[2]]
    clust_vid_new[[2]] <- dummy
  }
  return(list(vid_old = clust_vid_new[[1]], vid_new = clust_vid_new[[2]], eid_cutted = eid_cutted,
              clust_old = clust_split))
}

# membership: cluster_id -> vid
# cluster: vid -> cluster_id
# have a constraint that new clster size is always less than old cluster size
splitCluster.loho.waveletproposal <- function(mstgraph, k, membership, csize, edge_status, cluster) {
  
  n = length(V(mstgraph))
  indicator_within = (edge_status[E(mstgraph)$eid] == T) #length should be n-k, may be improved not using char vector
  proposal_weight = E(mstgraph)$proposalweight*indicator_within # between cluster edge weight will be 0
  sampled_idx = sample.int(n-1, 1, prob = proposal_weight) # although it is n-1, not n-k, some of the proposal_weights are zero
  # split proposal probability
  proposalprob_split = proposal_weight[sampled_idx]/sum(proposal_weight)
  # reverse move proposal probability
  indicator_between = (edge_status[E(mstgraph)$eid] == F)
  proposal_weight_rev = 1 - E(mstgraph)$proposalweight[indicator_between]
  proposalprob_merge = (1- proposal_weight[sampled_idx])/sum(proposal_weight_rev,1- proposal_weight[sampled_idx])
  
  edge_cutted = E(mstgraph)[sampled_idx]
  eid_cutted = edge_cutted$eid # see...
  #browser()
  edge_cutted_ends = ends(mstgraph, edge_cutted)
  cluster_idx = cluster[edge_cutted_ends]
  if(cluster_idx[1]!=cluster_idx[2]) stop("it cutted between-edges...")
  clust_split = cluster_idx[1]
  mst_subgraph = induced_subgraph(mstgraph, membership[[cluster_idx[1]]])
  
  mst_subgraph = delete.edges(mst_subgraph, which(E(mst_subgraph)$eid==eid_cutted))
  connect_comp = components(mst_subgraph)
  clust_vid_new = split(V(mst_subgraph)$vid, connect_comp$membership)
  
  
  # vid_new = (V(mst_subgraph)$vid)[cluster_new == 2]  # vid for vertices belonging to new cluster
  # membership[vid_new] = k + 1
  
  #introduce the constraint; |new cluster| <= |old cluster| where |.| is size
  if(length(clust_vid_new[[2]])>length(clust_vid_new[[1]])){
    dummy <- clust_vid_new[[1]]
    clust_vid_new[[1]] <- clust_vid_new[[2]]
    clust_vid_new[[2]] <- dummy
  }
  return(list(vid_old = clust_vid_new[[1]], vid_new = clust_vid_new[[2]], eid_cutted = eid_cutted,
              clust_old = clust_split, proposalprob=list(split= proposalprob_split, merge = proposalprob_merge)))
}

# function to merge two existing clusters
# always have a constraint that merged cluster inherits the larger cluster
mergeCluster.loho.waveletproposal <- function(mstgraph, eid_btw_mst, membership, csize, cluster, edge_list, edge_status,
                              change = F) {
  n = length(V(mstgraph))
  # edge for merging
  #edge_merge = sample.int(length(eid_btw_mst), 1)
  indicator_between = (edge_status[E(mstgraph)$eid] == T) #length should be n-k, may be improved not using char vector
  proposal_weight = E(mstgraph)$proposalweight*indicator_between
  sampled_idx = sample.int(n-1, 1, prob = proposal_weight)
  proposalprob_merge = proposal_weight[sampled_idx]/sum(proposal_weight)
  # reverse move proposal probability
  indicator_within = (edge_status[E(mstgraph)$eid] == F)
  proposal_weight_rev = 1 - E(mstgraph)$proposalweight[indicator_within]
  proposalprob_split = (1- proposal_weight[sampled_idx])/sum(proposal_weight_rev,1- proposal_weight[sampled_idx])
  
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = E(mstgraph)[sampled_idx]$eid
  edge_merge = which(eid_btw_mst == eid_merge)
  clusters_merge = cluster[edge_list[eid_merge, ]]
  #bigcluster = which.max(csize[clusters_merge])
  #c1 = clusters_merge[bigcluster]; c2 = clusters_merge[-bigcluster] # note that size of c1 > size of c2
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  if(c1 == c2) stop("merging cluster selection error")
  # merge c2 to c1
  
  # vid of vertices in c2
  vid_old = membership[[c2]]
  
  #added : vid of verticies in c1
  vid_new = membership[[c1]]
  ## #added
  if(length(vid_old) > length(vid_new)){
    oldinherit <- T
  }else{
    #cat("OLDINHERIT FALSE\n")
    oldinherit <- F
  }
  
  csize_new = NULL; clust_vid = NULL
  if(change) {
    clust_vid = membership
    clust_vid[[c1]] = c(membership[[c1]], vid_old)
    clust_vid[[c2]] = NULL
    
    csize_new = csize
    csize_new[c1] = length(clust_vid[[c1]])
    csize_new = csize_new[-c2]
  }
  # now drop c2
  return(list(vid_old = vid_old, clust_old = c2, clust_new = c1,
              edge_merge = edge_merge, clust_vid = clust_vid, csize = csize_new, oldinherit = oldinherit, vid_new = vid_new, 
              proposalprob=list(split= proposalprob_split, merge = proposalprob_merge)))
}


# function to merge two existing clusters
# always have a constraint that merged cluster inherits the larger cluster
mergeCluster.loho <- function(mstgraph, eid_btw_mst, membership, csize, cluster, edge_list,
                                  change = F) {
  # edge for merging
  edge_merge = sample.int(length(eid_btw_mst), 1)
  # update cluster information
  # clusters of endpoints of edge_merge
  eid_merge = eid_btw_mst[edge_merge]
  clusters_merge = cluster[edge_list[eid_merge, ]]
  
  #bigcluster = which.max(csize[clusters_merge])
  #c1 = clusters_merge[bigcluster]; c2 = clusters_merge[-bigcluster] # note that size of c1 > size of c2
  clusters_merge = sort(clusters_merge)
  c1 = clusters_merge[1]; c2 = clusters_merge[2] # note c1 < c2
  if(c1 == c2) stop("merging cluster selection error")
  # merge c2 to c1
  
  # vid of vertices in c2
  vid_old = membership[[c2]]
  
  #added : vid of verticies in c1
  vid_new = membership[[c1]]
  ## #added
  if(length(vid_old) > length(vid_new)){
    oldinherit <- T
  }else{
    #cat("OLDINHERIT FALSE\n")
    oldinherit <- F
  }
  
  csize_new = NULL; clust_vid = NULL
  if(change) {
    clust_vid = membership
    clust_vid[[c1]] = c(membership[[c1]], vid_old)
    clust_vid[[c2]] = NULL
    
    csize_new = csize
    csize_new[c1] = length(clust_vid[[c1]])
    csize_new = csize_new[-c2]
  }
  # now drop c2
  return(list(vid_old = vid_old, clust_old = c2, clust_new = c1,
              edge_merge = edge_merge, clust_vid = clust_vid, csize = csize_new, oldinherit = oldinherit, vid_new = vid_new))
}


# Auxiliary functions if the proposal is accepted --------------------

# function to update if a split move is accepted
updateSplit <- function(split_res, membership, k, csize, eid_btw_mst, cluster, edge_status,
                        adj_list=NULL, adj_edge_list=NULL) {
  clust_split = split_res$clust_old
  vid_old = split_res$vid_old; vid_new = split_res$vid_new
  
  membership[[clust_split]] = vid_old  # vid left in old cluster
  membership[[k+1]] = vid_new # vid in new cluster
  
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  cluster[vid_new] = k + 1
  eid_btw_mst = c(eid_btw_mst, split_res$eid_cutted)
  
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_btw = which(clust_adj_old != clust_split)
  eid_btw = adj_eid_old[idx_btw]
  edge_status[eid_btw] = T
  
  return(list(clust_vid = membership, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}


# function to update if a merge move is accepted
updateMerge <- function(res_merge, membership, csize, eid_btw_mst, cluster, edge_status,
                        adj_list=NULL, adj_edge_list=NULL) {
  clust_old = res_merge$clust_old; clust_new = res_merge$clust_new
  vid_old = membership[[clust_old]]
  membership[[clust_new]] = c(membership[[clust_new]], vid_old)
  membership[[clust_old]] = NULL
  
  csize[clust_new] = length(membership[[clust_new]])
  csize = csize[-clust_old]
  
  cluster[vid_old] = clust_new
  idx = which(cluster > clust_old)
  cluster[idx] = cluster[idx] - 1
  
  eid_btw_mst = eid_btw_mst[-res_merge$edge_merge]
  
  # update edge status
  adj_vid_old = unlist(adj_list[vid_old])
  adj_eid_old = unlist(adj_edge_list[vid_old])
  clust_adj_old = cluster[adj_vid_old]
  idx_within = which(clust_adj_old == clust_new)
  eid_within = adj_eid_old[idx_within]
  edge_status[eid_within] = F
  
  return(list(clust_vid = membership, csize = csize, cluster = cluster,
              eid_btw_mst = eid_btw_mst, estatus = edge_status))
}


# function to update Xtilde if move is accepted
updateXtilde.loho <- function(move, Xtilde, idx_old, idx_new, oldcol = NULL, newcol = NULL, csize = NULL){
  K = ncol(Xtilde)
  if(move == 'split'){
    Xtilde[,idx_old] <- oldcol
    if(idx_new <= K){
      Xtilde[,idx_new] <- newcol
    }else{
      Xtilde = cbind(Xtilde, newcol)
    }
  }
  if(move == 'merge') {
    Xtilde[, idx_new] = (Xtilde[, idx_new]*sqrt(csize[idx_new]) + Xtilde[, idx_old]*sqrt(csize[idx_old]))/sqrt(csize[idx_new]+csize[idx_old])
    Xtilde = Xtilde[, -idx_old, drop = F]
  }
  return(Xtilde)
}


# Auxiliary functions which are not used frequently but needed ---------------

# get Phi matrix
getPhi <- function(clust_vid, p){
  K = length(clust_vid)
  Phi = matrix(0, nrow = K, ncol = p)
  for(k in 1:K){
    Phi[k, clust_vid[[k]]] = 1/sqrt(length(clust_vid[[k]]))
  }
  #Phi = Matrix(Phi)
  return(Phi)
}

# get Phi*y, where Phi is k by p and y is p by 1 
# avoid construction of Phi and multiplication, which takes O(kp)
getPhi_y <- function(clust_vid, p, y){
  K = length(clust_vid)
  Phi_y = numeric(k)
  for(k in 1:K){
    Phi_y[k] = sum(y[clust_vid[[k]]]/sqrt(length(clust_vid[[k]])))
  }
  return(Phi_y)
}




# function to get whether an edge is within a cluster or bewteen two clusters
# 1 is cross-cluster(between) edge
# 0 is within-cluster edge
getEdgeStatus <- function(membership, inc_mat) {
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  iscrossing = rep(FALSE, nrow(inc_mat))
  iscrossing[membership_head != membership_tail] = TRUE # crossing
  return(iscrossing)
}



# David Dahl's method to get cluster point estimate from MCMC samples
Dahl <- function(s_save) {
  # Dahl’s method is equivalent to minimizing Binder’s loss function
  pihat <- salso::psm(s_save)
  estimate <- s_save[which.min(salso::binder(s_save, pihat)), ]
  estimate
  #return(aricode::NMI(true.cluster, estimate))
}

## functions for integrating out lambda





logHSmarginal <- function(y, sigma2, tau2){
  n = length(y)
  ybar = mean(y)
  logC = -(n-1)/2*log(2*pi*sigma2) - log(n)/2 - 1/(2*sigma2)*sum((y-ybar)^2)  # log normalizing constant
  integral = stats::integrate(integrandft, 0, Inf, ybar = ybar, t = sqrt(tau2), data.var = data.var, n = n, rel.tol = 1e-8)$value
  log(integral) + logC
}
integrandft <- function(u, ybar, t, data.var, n){
  dnorm(ybar, sd = sqrt(data.var*(1+u^2*t^2)/n))*2/(pi*(1+u^2))
}






