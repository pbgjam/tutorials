


bayesReg <- function(formula, data, ng = 5000, TOBIT=NULL){
  
  fc <- as.character(formula)
  
  yy <- unlist( strsplit( fc, '~' ) )
  yy <- yy[ nchar(yy) > 0]
  y  <- data[,yy[1]]
  
  ypos  <- which(y > 0)
  yzero <- which(y == 0)
  nzero <- length(yzero)
  npos  <- length(ypos)
  
  if(is.null(TOBIT)){
    TOBIT <- F
    if(nzero > 0)TOBIT <- T
  } 
  
  if(TOBIT)message('fitted as Tobit model')
  
  tmp <- model.frame(formula, data, na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames    <- colnames(x)
  snames    <- colnames(y)
  Q         <- ncol(x)
  n         <- nrow(x)
  predXcols <- 2:Q
  
  ymiss <- which(is.na(y))

  yy <- y
  xx <- x
  if(length(ymiss) > 0){
    yy <- yy[-ymiss]
    xx <- xx[-ymiss,]
  }

  XX  <- crossprod(xx)
  IXX <- solve(XX)
  bg  <- IXX%*%crossprod(xx,yy)
  py  <- x%*%bg
  w   <- y
  w[ymiss] <- mean(y,na.rm=T)
  if(TOBIT)w[w == 0] <- py[w == 0]
  wi  <- c(which(y == 0),ymiss)   
  XX  <- crossprod(x)
  
  
  priorB   <- matrix(0,Q,1)         #prior mean regression parameters
  priorIVB <- solve(diag(1000,Q))   #inverse prior covariance matrix
  s1       <- .1                    #variance prior values
  s2       <- .1
    
  bchains <- matrix(NA,ng,Q)
  schains <- rep(0,ng)  #store parameters
  colnames(bchains) <- xnames
  ychains <- matrix(0,ng,n)
  
  for(g in 1:ng){
    
    sg <- updateSigma(x, w, bg, s1, s2)
    bg <- updateBeta(x, w, sg, priorIVB, priorB, XX)
    mu <- x%*%bg
    py <- rnorm(n,mu,sqrt(sg))
    
    if(TOBIT){
      w[yzero] <- .tnorm(nzero,-Inf,0,mu[yzero],sqrt(sg))
      py[py <= 0] <- 0
    }
      
    w[ymiss] <- py[ymiss]
    
    bchains[g,] <- bg   #store estimates
    schains[g]  <- sg
    ychains[g,] <- py
  }
  
  beta <- signif( t( apply(bchains,2,quantile,c(.5,.025,.975)) ), 4)
  zero <- which(beta[,2] < 0 & beta[,3] > 0)
  notZero <- rep('*',Q)
  notZero[ zero ] <- ' '
  
  colnames(beta) <- c('median','0.025','0.975')
  
  py <- signif( t( apply( ychains, 2, quantile,c(.5,.025,.975) ) ), 3)
  py[,1] <- colMeans(ychains)
  py <- cbind(y,py)
  
  rmspe <- sqrt( mean((y - py[,2])^2, na.rm=T) )
  
  list(beta = beta, predict = py, sigma = median(schains), 
       rmspe = rmspe)
}
  
updateSigma <- function(x, y, beta, s1, s2){ # random value for residual variance
  n  <- length(y)
  u1 <- s1 + n/2
  u2 <- s2 + 0.5*crossprod(y - x%*%beta)
  1/rgamma(1,u1,u2)                          
}

updateBeta <- function(x, y, sigma, priorIVB, priorB, 
                       XX=NULL){  # random vector of coefficients
  
  if(is.null(XX))XX <- crossprod(x)
  
  V  <- solve( XX/sigma + priorIVB ) 
  v  <- crossprod(x,y)/sigma + priorIVB%*%priorB
  t( myrmvnorm(1,t(V%*%v),V) )                     
}

  myrmvnorm <- function (n, mu, sigma){
    
    # n - no. samples from one mu vector or nrow(mu) for matrix
    # mu - mean vector
    # sigma - variance
    
    if(!is.matrix(mu))mu <- matrix(mu,1)
    if(ncol(mu) == 1) mu <- t(mu)
    
    m <- ncol(sigma)
    
    if(n > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
    
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    
    mu + matrix(rnorm(n * m), nrow = n) %*% retval
  }
  
.tnorm <- function(n,lo,hi,mu,sig){   
    
    #normal truncated lo and hi
    
    tiny <- 10e-6
    
    if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
    if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
    
    q1 <- pnorm(lo,mu,sig)
    q2 <- pnorm(hi,mu,sig) 
    
    z <- runif(n,q1,q2)
    z <- qnorm(z,mu,sig)
    
    z[z == Inf]  <- lo[z == Inf] + tiny
    z[z == -Inf] <- hi[z == -Inf] - tiny
    z
  }
  