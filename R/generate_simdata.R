# generating data for the simulation. Assuming 30 x 30 lattice graph. 
generate_simdata <- function(n = 100, rhox = 0, SNR = 4, option = "twoclusters", ntest = 1000, seed = NULL){
  
  set.seed(seed)
  
  m = 30
  m2 = m^2
  beta <- numeric(m2)
  
  if(option == "twoclusters"){
    for(j in 0:5) beta[(460:465)+(j*30)] <- 1
    for(j in 0:5) beta[(286:291)+(j*30)] <- -1
    for(j in 0:5) beta[(280:285)+(j*30)] <- 1
    for(j in 0:5) beta[(466:471)+(j*30)] <- -1
    for(j in 0:5) beta[(283:285)+(j*30)] <- -1
    for(j in 0:5) beta[(466:468)+(j*30)] <- 1
  }else if(option == "fourclusters"){
    for(j in 0:5) beta[(458:465)+(j*30)] <- 1
    for(j in 0:5) beta[(286:293)+(j*30)] <- 2
    for(j in (-2):5) beta[(280:285)+(j*30)] <- -1
    for(j in 0:7) beta[(466:471)+(j*30)] <- -2
  }else if(option == "eightclusters"){
    for(j in 0:5) beta[(458:465)+(j*30)] <- 1
    for(j in 0:5) beta[(286:293)+(j*30)] <- 2
    for(j in (-2):5) beta[(280:285)+(j*30)] <- -1
    for(j in 0:7) beta[(466:471)+(j*30)] <- -2
    for(j in 0:5) beta[(1:5)+(j*30)] <- 3
    for(j in 0:5) beta[(26:30)+(j*30)] <- -3
    for(j in (-4):0) beta[(871:875)+(j*30)] <- -3
    for(j in (-4):0) beta[(896:900)+(j*30)] <- 3
  }else{
    stop("provide one of the options: 1. 'twoclusters', 2. 'fourclusters', 3. 'eightclusters' ")
  }
  
  if(rhox > 0){  
    D   <- as.matrix(dist(1:m))
    S   <- exp(-(D/rhox)^1)
    P   <- t(chol(S))
    P   <- kronecker(P,P)
  }
  
  # THIS PART IS NEEDED FOR FITTING Soft-thresholded GP(STGP)
  # comment the following if you are not fitting STGP of Kang et al.(2018)
  
  # Set up the basis functions when rhox > 0:
  knots <- expand.grid(seq(0,m+1,length=round(m/2)),
                      seq(0,m+1,length=round(m/2)))
  bw    <- min(dist(knots))
  # Define the spatial adjacency matrix: 1D
  ADJ1D <- function(m){
    A <- 0*diag(m)
    for(j in 2:m){
      A[j,j-1]<-A[j-1,j]<-1
    }
    A
  }

  A     <- ADJ1D(round(m/2))
  ADJ   <- kronecker(A,diag(round(m/2)))+kronecker(diag(round(m/2)),A)
  
  s     <- expand.grid(1:m,1:m)
  D     <- fields::rdist(s,knots)
  B     <- exp(-0.5*(D/bw)^2)
  B     <- ifelse(D>3*bw,0,B)

  S     <- sqrt(diag(B%*%solve(diag(colSums(ADJ))-0.99*ADJ)%*%t(B)))
  B     <- diag(1/S)%*%B
  S     <- sqrt(diag(B%*%solve(diag(colSums(ADJ))-0.99*ADJ)%*%t(B)))
  
  if(rhox > 0){
    X <- mvnfast::rmvn(n, mu = rep(0,m2), sigma = t(P), isChol = T)
    Xtest <- mvnfast::rmvn(ntest, mu = rep(0,m2), sigma = t(P), isChol = T)
  }else{
    X <- matrix(rnorm(n*m2), n, m2) # no dependency
    Xtest <- matrix(rnorm(ntest*m2), ntest, m2)
  }
  
  X = scale(X)/sqrt(n-1) # column standardization
  
  sigma2 = var(X%*%beta)/SNR
  
  Y <- as.numeric(X%*%beta+rnorm(n,0,sqrt(sigma2)))
  Ytest <- as.numeric(Xtest%*%beta+rnorm(ntest,0,sqrt(sigma2)))
  
  list(Y=Y, X=X, beta.true = beta, sigma2.true = sigma2, Xtest = Xtest, Ytest = Ytest, rhox = rhox, SNR = SNR)
}
