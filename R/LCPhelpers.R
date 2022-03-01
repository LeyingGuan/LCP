#'@export
dist <- function(x, x1){
  if(is.null(dim(x))){
    tmp = abs(x-x1)
  }else if(dim(x)[2]==1){
    tmp = abs(x-x1)
  }else{
    tmp = rep(0, nrow(x))
    for(j in 1:ncol(x)){
      tmp = tmp+(x[,j] - x1[j])^2
    }
    tmp = tmp/ncol(x)
    tmp = sqrt(tmp)
  }
  tmp
}

#'@export
train_test_split <- function(x, y, test_ratio, random_state = NULL){
  if (!is.null(random_state)) {
    # reinstate system seed after simulation
    sysSeed <- .GlobalEnv$.Random.seed
    on.exit({
      if (!is.null(sysSeed)) {
        .GlobalEnv$.Random.seed <- sysSeed
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  }
  n = length(y)
  p = ncol(x)
  lltest = ceiling(n * test_ratio)
  lltest = sample(1:n, lltest)
  lltrain = setdiff(1:n, ll_test)
  xtrain = x[lltrain,]
  ytrain = y[lltrain]
  xtest = x[lltest,]
  ytest = y[lltest]
  if(p == 1){
    xtest = matrix(xtest, ncol = p)
    xtrain = matrix(xtrain, ncol = p)
  }
  return(list(xtrain = xtrian, ytrain= ytrain,
              xtest = xtest, ytest = ytest))
}


#'@export
LCPdefault_distance <- function(utangent, uorthogonal, estimated_sds, xtrain, xcalibration, xtest){
  n0 =nrow(xcalibration); m = nrow(xtest)
  z =  estimated_sds[[1]]
  z0 =  estimated_sds[[2]]
  z1 =  estimated_sds[[3]]
  n = length(z)
  Hdist1 = array(0, dim = c(n0+m,n0+m))
  Hdist2 = array(0, dim = c(n0+m,n0+m))
  Hdist3 = array(0, dim = c(n0+m,n0+m))
  Hdist1cv = array(0, dim = c(n,n))
  Hdist2cv = array(0, dim = c(n,n))
  Hdist3cv = array(0, dim = c(n,n))
  

  for(j in 1:n0){
    Hdist1[1:n0,j] =dist(z0, z0[j])
    Hdist1[(n0+1):(n0+m),j] =dist(z1, z0[j])
  }
  for(j in 1:m){
    Hdist1[1:n0,j+n0] =dist(z0, z1[j])
    Hdist1[(n0+1):(n0+m),j+n0] =dist(z1, z1[j])
  }
  
  for(j in 1:n){
    Hdist1cv[1:n,j] =dist(z, z[j])
  }

  if(!is.null(utangent)){
    z0parallel = xcalibration%*%utangent
    z1parallel = xtest%*%utangent
    zparallel = xtrain%*%utangent
    if(ncol(z0parallel)==1){
      for(j in 1:n0){
        Hdist2[1:n0,j] =dist(z0parallel, z0parallel[j])
        Hdist2[(n0+1):(n0+m),j] =dist(z1parallel , z0parallel[j])
        
      }
      for(j in 1:m){
        Hdist2[1:n0,j+n0] =dist(z0parallel, z1parallel[j])
        Hdist2[(n0+1):(n0+m),j+n0] =dist(z1parallel , z1parallel[j])
      }
    }else{
      for(j in 1:n0){
        Hdist2[1:n0,j] =dist(z0parallel, z0parallel[j,])
        Hdist2[(n0+1):(n0+m),j] =dist(z1parallel , z0parallel[j,])
        
      }
      for(j in 1:m){
        Hdist2[1:n0,j+n0] =dist(z0parallel, z1parallel[j,])
        Hdist2[(n0+1):(n0+m),j+n0] =dist(z1parallel , z1parallel[j,])
      }
    }
    if(ncol(zparallel)==1){
      for(j in 1:n){
        Hdist2cv[1:n,j] =dist(zparallel , zparallel[j])
      }
    }else{
      for(j in 1:n){
        Hdist2cv[1:n,j] =dist(zparallel , zparallel[j,1])
      }
    }
    
  }
  if(!is.null(uorthogonal)){
    z0orthogonal = xcalibration%*%uorthogonal
    z1orthogonal = xtest%*%uorthogonal
    zorthogonal = xtrain%*%uorthogonal 
    if(ncol(z0orthogonal)==1){
      for(j in 1:n0){
        Hdist3[1:n0,j] =dist(z0orthogonal, z0orthogonal[j])
        Hdist3[(n0+1):(n0+m),j] =dist(z1orthogonal , z0orthogonal[j])
        
      }
      for(j in 1:m){
        Hdist3[1:n0,j+n0] =dist(z0orthogonal, z1orthogonal[j])
        Hdist3[(n0+1):(n0+m),j+n0] =dist(z1orthogonal , z1orthogonal[j])
      }
    }else{
      for(j in 1:n0){
        Hdist3[1:n0,j] =dist(z0orthogonal, z0orthogonal[j,])
        Hdist3[(n0+1):(n0+m),j] =dist(z1orthogonal , z0orthogonal[j,])
      }
      for(j in 1:m){
        Hdist3[1:n0,j+n0] =dist(z0orthogonal, z1orthogonal[j,])
        Hdist3[(n0+1):(n0+m),j+n0] =dist(z1orthogonal , z1orthogonal[j,])
      }
    }
    
    if(ncol(zorthogonal)==1){
      for(j in 1:n){
        Hdist3cv[1:n,j] =dist(zorthogonal , zorthogonal[j])
      }
    }else{
      for(j in 1:n){
        Hdist3cv[1:n,j] =dist(zorthogonal , zorthogonal[j,])
      }
    }
  }



  tmp = (row(Hdist1cv)!=col(Hdist1cv))
  s01 = mean( Hdist1cv[tmp]); 
  if(!is.null(utangent)){
    s02 = mean(Hdist2cv[tmp]);
  }else{
    s02 = 0
  }
  if(!is.null(uorthogonal)){
    s03 = mean(Hdist3cv[tmp])
  }else{
    s03 = 0
  }
  w = s02/(s02+s03)
  if(is.null(uorthogonal)){
    w  = 0
  }

  Hdist23 = (1-w)*Hdist2 + w*Hdist3
  Hdist23cv = (1-w) * Hdist2cv + w * Hdist3cv
  s023 = mean(Hdist23cv[tmp])
  a = 1
  Hcv =   a*Hdist1cv/s01+Hdist23cv[1:n,1:n]/s023; 
  tmp = a*Hdist1/s01+Hdist23/s023
  H = tmp[1:n0, 1:n0]
  Hnew1 =  tmp[-(1:n0),1:n0]; HnewT1 =  tmp[(1:n0),-(1:n0)]
  return(list(Hcv = Hcv, H = H, Hnew = Hnew1, HnewT = HnewT1))
}



#' sim_data_generator_1D_example2 - generate 1D simulated data for example2
#'@export
sim_data_generator_1D_example2 <- function(sim_name, n = 1000, n0 = 1000, m = 1000, alpha = 0.05){
  if(sim_name == "1D_setA"){
    #sin
    noise_generating = function(x){
      abs(sin(x))
    }
  }else if(sim_name == "1D_setB"){
    #cos
    noise_generating = function(x){
      abs(cos(x))
    }
    
  }else if(sim_name == "1D_setC"){
    #linear
    noise_generating = function(x){
      abs(x)
    }

  }else if(sim_name == "1D_setD"){
    #constant
    noise_generating = function(x){
      rep(1,length(x))
    }
  }
  x <- runif(n,-2,2);
  s <- noise_generating(x);
  y <- s*rnorm(n);
  x1 <- runif(m,-2,2); x0 <-runif(n0,-2,2)
  x1 <- x1[order(x1)];
  s1 <-noise_generating(x1); s0 = noise_generating(x0)
  y1 = s1*rnorm(m); y0 = s0 *rnorm(n0)
  xtrain = matrix(x, ncol = 1)
  xcalibration = matrix(x0, ncol = 1)
  xtest = matrix(x1, ncol = 1)
  ytrain = matrix(y, ncol = 1)
  ycalibration = matrix(y0, ncol = 1)
  ytest = matrix(y1, ncol = 1)
  truePI = matrix(0, ncol = 2, nrow = m)
  truePI[,2] = qnorm(1-alpha/2)*s1
  truePI[,1] = qnorm(alpha/2)*s1
  return(list(xtrain = xtrain, ytrain = ytrain,
              xcalibration = xcalibration, ycalibration = ycalibration,
              xtest = xtest, ytest = ytest, truePI = truePI))
}

#' sim_data_generator_1D_example1 - generate 1D simulated data for example1
#'@export
sim_data_generator_1D_example1 <- function(sim_name, n = 1000, n0 = 1000, m = 1000, alpha = 0.05){
  if(sim_name == "1D_setA"){
    #sin
    noise_generating = function(x){
      abs(sin(x))
    }
  }else if(sim_name == "1D_setB"){
    #cos
    noise_generating = function(x){
      abs(cos(x))
    }
    
  }else if(sim_name == "1D_setC"){
    #linear
    noise_generating = function(x){
      abs(x)
    }
    
  }else if(sim_name == "1D_setD"){
    #constant
    noise_generating = function(x){
      rep(1,length(x))
    }
  }
  x <- rnorm(n = n, mean = 0, sd = 1)
  s <- noise_generating(x);
  y <- s*rnorm(n);
  x1 <- rnorm(n = m, mean = 0, sd = 1) ; x0 <- rnorm(n = n0, mean = 0, sd = 1)
  x1 <- x1[order(x1)];
  s1 <-noise_generating(x1); s0 = noise_generating(x0)
  y1 = s1*rnorm(m); y0 = s0 *rnorm(n0)
  xtrain = matrix(x, ncol = 1)
  xcalibration = matrix(x0, ncol = 1)
  xtest = matrix(x1, ncol = 1)
  ytrain = matrix(y, ncol = 1)
  ycalibration = matrix(y0, ncol = 1)
  ytest = matrix(y1, ncol = 1)
  truePI = matrix(0, ncol = 2, nrow = m)
  truePI[,2] = qnorm(1-alpha/2)*s1
  truePI[,1] = qnorm(alpha/2)*s1
  return(list(xtrain = xtrain, ytrain = ytrain,
              xcalibration = xcalibration, ycalibration = ycalibration,
              xtest = xtest, ytest = ytest, truePI = truePI))
}
