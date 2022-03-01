library(R6)
LCP_alpha <- function(vs, alphas, alpha){
  if(is.null(dim(alphas))){
    idxes =which(alphas< alpha)
    l1 = length(idxes)
    if(l1>0){
      return(vs[idxes[l1]])
    }else{
      return(Inf)
    }
  }else{
    idxes = apply(alphas<alpha, 1, which)
    v_bounds = rep(Inf, nrow(alphas))
    tmp = sum(class(idxes)=="matrix")
    if(tmp == 0){
      for(i in 1:length(idxes)){
        l1 =length(idxes[[i]])
        if( l1> 0){
          v_bounds[i] =vs[idxes[[i]][l1]]
        }
      }
    }else{
      for(i in 1:ncol(idxes)){
        l1 =length(idxes[,i])
        if( l1> 0){
          v_bounds[i] =vs[idxes[,i][l1]]
        }
      }
    }

    return(v_bounds)
  }
}


autoTune_distance <- function(V, n, hs, D, alpha = 0.05,  delta = 0.05, B = 10, trace = T, lambda = 1){
  n0 = length(V)
  J = length(hs)
  n1 = min(n0, n+1)
  if(n1>= 2/3 * n0){
    B1 = 1
  }else{
    B1 = B
  }
  C1 = array(0, dim = c(B1, 2, J))
  C2 = array(0, dim = c(n0, B, J))
  mu = array(0, dim = c(n0, J))
  varphi = array(Inf, dim = c(n0, B, J))
  Ts = array(0, dim = c(J, 3))
  orders = order(V)
  V = V[orders]
  D = D[orders, orders]
  thr = max(V)
  print("###estimate means######")
  for(b in 1:B1){
    if(trace){
      print(b)
    }
    I = sort(sample(1:n0, n1, replace = F))
    Vb = V[I]
    Db = D[I, I]
    for(j in 1:J){
      h = hs[j]
      H = exp(-Db/h)
      Qcumsumb = t(apply(H,1,cumsum))
      ret1 = LCP_construction_distance_loop(V = Vb,  Qcumsum = Qcumsumb, H = H)
      tmp = rep(0, n0)
      for(i in 1:n0){
        Vb0 = Vb[-i]
        tmp[i] = LCP_alpha( vs = c(Vb0, Inf), alphas = ret1$alphas[i,], alpha = 1-alpha)
      }
      C1[b,1,j] = mean(tmp > thr)
      C1[b,2,j] = mean(tmp[tmp<=thr])
    }
  }
  print("###standard deviation estimation#####")
  if(lambda > 0){
    for(b in 1:B){
      if(trace){
        print(b)
      }
      I = sort(sample(1:n0, n, replace = T))
      Vb = V[I]
      Db = D[I,I]
      Dnew = D[, I]
      DnewT = D[I,]
      Vb = c(Vb, Inf)
      id_low = id_low_search(Vb)
      for(j in 1:J){
        h = hs[j]
        H = exp(-Db/h)
        Hnew = exp(-Dnew/h)
        HnewT = exp(-DnewT/h)
        Qcumsumb = t(apply(H,1,cumsum))
        q_low = q_low_compute(id_low[-length(id_low)],  Qcumsumb)
        qn =  Qcumsumb[,n]
        ret = LCP_construction_path(alpha = alpha, V = Vb,  id_low = id_low, q_low = q_low,  qn = qn, 
                                    Hnew = Hnew, HnewT = HnewT, type = "distance")
        tmp = LCP_alpha( vs = Vb, alphas = ret$Smat, alpha = 1-alpha)
        C2[,b,j] = tmp
      }
    }
    for(j in 1:J){
      for(i in 1:n0){
        tmp = C2[i,,j]
        tmp1 = mean(tmp <= thr)
        if(tmp1 > 0){
          mu[i,j] = mean(tmp[tmp<=thr])
        }
        for(b in 1:B){
          if(C2[i,b,j] < Inf){
            varphi[i,b,j] = (C2[i,b,j] - mu[i,j])^2
          }
        }
      }
    }
  }
  Ts[,1] = hs
  Ts[,2] = apply(C1[,1,, drop = FALSE],3,mean)
  Ts[,3] = apply(C1[,2,, drop = FALSE],3,mean)
  intermediants = list(probs = Ts[,2], means = Ts[,3], sds = rep(0, J))
  if(lambda > 0){
    for(j in 1:J){
      tmp = varphi[,,j]
      intermediants$sds[j] = sqrt(mean(tmp[tmp < Inf]))
      Ts[j,3] = Ts[j,3]+lambda*intermediants$sds[j]
    }
  }
  idx = which(Ts[,2] <= delta)
  if(length(idx) == 0){
    h = NA
  }else{
    hs0 = hs[idx]
    c0 = Ts[idx,3]
    h = hs0[which.min(c0)]
  }
  return(list(Ts = Ts, h = h, intermediants = intermediants))
  
}

LCP_construction_path = function(alpha, V,  id_low, q_low,  qn,  Hnew, HnewT,type = "distance",
                                 neighbor_size = 100, idx_boundary = NULL, size_boundary= NULL,distance_boundary = NULL){
  if(type == "neighbor"){
    if(is.null(idx_boundary) | is.null(distance_boundary)){
      stop("boundary conditions for nearest neighbor localizers missing.")
    }
    m = nrow(Hnew)
    n = length(V)
    n0 = n-1
    Smat = matrix(0, nrow = m, ncol = n)
    deltaLCP = rep(0, m)
    for(i in 1:m){
      hnew = rank(Hnew[i,], ties.method = "random")
      hnewT = HnewT[,i]
      hnew = ifelse(hnew<neighbor_size, 1, 0)
      qni = qn
      q_lowi = q_low 
      for(j in 1:n0){
        if(hnewT[j] < distance_boundary[j]){
          hnewT[j] = 1
          qni[j] = qni[j]-1
          if(idx_boundary[j] <= id_low[i]+1){
            q_lowi[j] = q_lowi[j]-1
          }
        }else if(hnewT[j] > distance_boundary[j]){
          hnewT[j] = 0
        }else{
          #hnewT[j] = distance_boundary[j]
          indicator = rbinom(n = 1, size = 1, prob = 1/(size_boundary[j]+1))
          if(indicator == 1){
            hnewT[j] = 1
            qni[j] = qni[j]-1
            if(idx_boundary[j] <= id_low[i]+1){
              q_lowi[j] = q_lowi[j]-1
            }
          }else{
            hnewT[j] = 0
          }
        }
      }
      hnew = matrix(hnew, nrow = 1)
      hnewT  = matrix(hnewT, ncol = 1)
      ret = LCP_construction_path_distance(V = V, id_low = id_low,
                                           q_low = q_lowi,qn = qni, Hnew = hnew, HnewT = hnewT)
      deltaLCP[i] = LCP_alpha( vs = V, alphas = ret$alphas[1,], alpha = 1-alpha)
      Smat[i,] = ret$alphas[1,]
    }

  }else{
    ret= LCP_construction_path_distance(V = V, id_low = id_low,
                                        q_low = q_low,qn = qn, Hnew = Hnew,
                                        HnewT = HnewT)
    deltaLCP = LCP_alpha( vs = V, alphas = ret$alphas, alpha = 1-alpha)
    Smat =  ret$alphas
    
  }
  return(list(deltaLCP = deltaLCP, Smat = Smat))
  
}


#'@import R6
#'@import Rcpp
#'@import XRPython
#'@useDynLib LCPcpp, .registration=TRUE
#'@export
#'\examples{
#' man/examples/example1.R
#'}
LCPmodule <- R6Class(classname = "LCP",
                     list(
                       #quantities and functions to feed in
                       ##quantities
                       H = NULL,
                       h = 1,
                       Hnew = NULL,
                       HnewT = NULL,
                       Hdistance = NULL,
                       Hrank = NULL,
                       idx_boundary = NULL,
                       distance_boundary = NULL,
                       size_boundary = NULL,
                       V = NULL,
                       n = NULL,
                       Qcumsum = NULL,
                       qlow0 = NULL,
                       qn0 = NULL,
                       alpha = NULL,
                       type = "distance",
                       ##functions
                       invert_func = NULL,
                       id_low = NULL,
                       band_V = NULL,
                       band_Y = NULL,
                       Smat = NULL,
                       
                       ##initialization
                       initialize = function(H, V, h = 1,alpha = 0.05, type = "distance", invert_func = NULL){
                         self$H = H
                         self$h = h
                         self$V = c(V, Inf)
                         self$n = length(self$V)-1
                         self$alpha = alpha
                         self$type = type
                         self$invert_func = invert_func
                         if(is.null(self$invert_func)){
                           #if no invert func is provided, use the identify function
                           self$invert_func = function(x){
                             return(x)
                           }
                         }
                       },
                       ##LCP inference functions
                       lower_idx = function(){
                         self$id_low = id_low_search(self$V)
                       },
                       cumsum_unnormalized = function(){
                         if(self$type == "distance"){
                           self$Hdistance = exp(-self$H/self$h)
                           self$Qcumsum = t(apply(self$Hdistance,1,cumsum))
                         }else if(self$type == "neighbor"){
                           self$Hrank = t(apply(self$H,1,rank, ties_method = "random"))
                           self$idx_boundary = apply(self$Hrank,1,function(z) which(z == self$h))
                           self$distance_boundary = apply(self$H, 1, function(z) sort(z)[self$h])
                           n = nrow(self$Hrank)
                           for(i in 1:n){
                             tmp = self$Hrank[i,]
                             self$Hrank[i,tmp  <= self$h] = 1
                             self$Hrank[i,tmp  >self$h] = 0
                           }
                           self$Qcumsum = t(apply(self$Hrank,1,cumsum))
                         }else{
                           stop("unsupported localizer type.")
                         }
                         self$qlow0 = q_low_compute(self$id_low[-length(self$id_low)],  self$Qcumsum)
                         self$qn0 = self$Qcumsum[,self$n]
                       },

                       LCP_construction = function(Hnew, HnewT){
                         if(self$type == "distance"){
                           Hnew = exp(-Hnew/self$h)
                           HnewT = exp(-HnewT/self$h)
                         }
                         if(is.null(dim(Hnew))){
                           Hnew = matrix(Hnew, nrow = 1)
                           HnewT  = matrix(HnewT, ncol = 1)
                         }
                         n = nrow(self$H)
                         m = nrow(Hnew)
                         if(self$type == "neighbor"){
                           ret = LCP_construction_path(alpha = self$alpha,V = self$V, id_low = self$id_low,
                                                       q_low = self$qlow0, qn = self$qn0, 
                                                       Hnew = Hnew, HnewT = HnewT, type = self$type,
                                                       neighbor_size = self$h, idx_boundary = self$idx_boundary,
                                                       size_boundary = self$size_boundary ,
                                                       distance_boundary = self$distance_boundary)
             
                         }else{

                           ret = LCP_construction_path(alpha = self$alpha, V = self$V, id_low = self$id_low,
                                                       q_low = self$qlow0, qn = self$qn0, 
                                                       Hnew = Hnew, HnewT = HnewT, type = self$type)
                         }
                         self$band_V = ret$deltaLCP
                         self$Smat = ret$Smat

                       },
                       
                       LCP_auto_tune = function(V0, H0, hs, B = 5, delta =0.05, lambda = 0, trace = TRUE){
                         if(self$type == "distance"){
                           ret = autoTune_distance(V = V0, n = self$n, hs = hs, D = H0, alpha = self$alpha,  
                                                   delta =delta, B = B, trace = trace, lambda = lambda)
                         }else{
                           stop("unsupported localizeer for auto-tune.")
                         }
                         return(ret)
                       }
                       
                     ) )










