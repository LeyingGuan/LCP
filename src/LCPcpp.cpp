#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



/*
 * Note here the V is extended with (n+1) elements and the last one is Infinity.
 */
//[[Rcpp::export]]
arma::vec id_low_search(arma::vec& V){
  int n = V.n_elem;
  arma::vec id_low = arma::ones<arma::vec>(n) * (-1);
  int it = -1;
  for(int i = 0; i < n; i++){
    while(V[it+1] < V[i]){
      it += 1;
    }
    id_low[i] = it;
  }
  return id_low;
}

//[[Rcpp::export]]
arma::vec q_low_compute(arma::vec& id_low, arma::mat& Qcumsum){
  int n = id_low.n_elem;
  arma::vec q_low = arma::zeros<arma::vec>(n);
  for(int i = 0; i< n;  i++){
    if(id_low(i) >= 0){
      q_low(i) = Qcumsum(i, id_low(i));
    }
  }
  return q_low;
}

/*
 * V is the extended (n+1) vector.
 * id_low: (n+1) vector
 * q_low, qn: n-vector. 
 */

//[[Rcpp::export]]
List LCP_construction_path_distance(arma::vec& V, arma::vec& id_low, arma::vec& q_low, arma::vec& qn, arma::mat& Hnew, arma::mat& HnewT){
  int m = Hnew.n_rows;
  int n = V.n_elem;
  int n0 = n-1;
  //note n: <-n+1, with infinite V at the end
  arma::mat theta = arma::zeros<arma::mat>(m, n);
  arma::mat theta0 = arma::zeros<arma::mat>(m, n);
  arma::vec qn_new =  arma::zeros<arma::vec>(m);
  arma::mat alphas = arma::zeros<arma::mat>(m, n);
  //calculate unnormalized ttheta in the paper
  //for j test sample, the unnormalized ttheta is ttheta(j,.)
  for(int i = 1; i < n; i++){
    theta.col(i) = theta.col(i-1) * 1.0;
    if((id_low(i) - id_low(i-1)) > 0){
      int i1 = (id_low(i-1)+1);
      while(i1 <=id_low(i)){
        theta.col(i) += Hnew.col(i1);
        i1 += 1;
      }
    }
  }
  qn_new = theta.col(n0) * 1.0;
  theta0 = theta;
  for(int i = 0; i < m; i++){
    arma::vec normalizers = qn + HnewT.col(i);
    //for each training sample->normalized theta_i^1
    arma::vec p1 = (q_low+HnewT.col(i))/normalizers;
    //for each training sample->normalized theta_i^2
    arma::vec p2 = q_low/normalizers;
    arma::vec ttheta_i = arma::zeros<arma::vec>(n);
    arma::vec ttheta_i0 = arma::zeros<arma::vec>(n0);
    for(int i1 = 0; i1 < n0; i1++){
      ttheta_i(i1) =  theta(i,i1);
      ttheta_i0(i1) = theta(i,i1);
    }
    ttheta_i(n0) =theta(i, n0);
    //get the normailzed ttheta for a test sample
    ttheta_i = ttheta_i/(qn_new(i)+1.0);
    ttheta_i0 = ttheta_i0/(qn_new(i)+1.0);
    arma::uvec idA1 = find(p1 < ttheta_i0);
    arma::uvec idA2 = find(p2 >= ttheta_i0);
    arma::uvec idA3 = find((p2 < ttheta_i0) && (p1 >= ttheta_i0));
    arma::vec pA1 = p1(idA1);
    arma::vec pA2 = p2(idA2);
    //arma::vec pA3 = idxs(idA3);
    arma::vec pA3 = id_low(idA3);
    int nA1 = pA1.n_elem;
    int nA2 = pA2.n_elem;
    int nA3 = pA3.n_elem;
    arma::vec pA1sort;
    arma::vec pA2sort;
    arma::vec pA3sort;
    if(nA1 > 0){
      pA1sort = arma::sort(pA1);
    }else{
      pA1sort = pA1;
    }
    if(nA2 > 0){
      pA2sort = arma::sort(pA2);
    }else{
      pA2sort = pA2;
    }
    if(nA3 > 0){
      pA3sort = arma::sort(pA3);
    }else{
      pA3sort = pA3;
    }
    int itA1=0;
    int itA2 = 0;
    int itA3 = 0;
    int k = -1;
    //Rcout << "step1" << "\n";
    while(k < (n-1)){
      k += 1;
      //check if someone in A list satisfies the condition
      if(nA1 > 0 & itA1 < nA1){
        while(itA1< (nA1-1) & pA1sort(itA1) < ttheta_i(k)){
          itA1 += 1;
        }
        if(itA1 == (nA1-1)){
          if(pA1sort(itA1)<  ttheta_i(k)){
            itA1 += 1;
          }
        }
      }
      //check if someone in B list satisfies the condition
      if(nA2 > 0 & itA2 < nA2){
        while(itA2< (nA2-1) & pA2sort(itA2) < ttheta_i(k)){
          itA2 += 1;
        }
        if(itA2 == (nA2-1)){
          if(pA2sort(itA2)< ttheta_i(k)){
            itA2 += 1;
          }
        }
      }
      //check if someone in B list satisfies the condition
      if(nA3 > 0 & itA3 < nA3){
        int id_low_k = id_low(k);
        while(itA3< (nA3-1) & pA3sort(itA3) < id_low_k){
          itA3 += 1;
        }
        if(itA3 == (nA3-1)){
          if(pA3sort(itA3) < id_low_k){
            itA3 += 1;
          }
        }
      }      
      // Calculate condition value
      alphas(i, k) = (itA1+itA2+itA3)/(n0+1.0);
    }
  }
  List ret = List::create(Named("alphas") = alphas , _["vs"] = V);
  return ret;
}


//[[Rcpp::export]]
List LCP_construction_distance_loop(arma::vec& V, arma::mat& Qcumsum, arma::mat& H){
  int n = V.n_elem;
  arma::mat alphas = arma::zeros<arma::mat>(n, n);
  arma::vec id_low = id_low_search(V);
  arma::vec q_low = q_low_compute(id_low, Qcumsum);
  for(int i = 0; i < n; i++){
    arma::vec Vi = arma::zeros<arma::vec>(n);
    arma::vec id_low_i = arma::zeros<arma::vec>(n);
    arma::vec q_low_i = arma::zeros<arma::vec>(n-1);
    arma::vec qn_i = arma::zeros<arma::vec>(n-1);
    arma::mat Hnew = arma::zeros<arma::mat>(1,n-1);
    arma::mat HnewT = arma::zeros<arma::mat>(n-1,1);
    for(int i1 = 0; i1 < n-1; i1++){
      int i2 = i1;
      if(i1 >= i){
        i2 += 1;
      }
      q_low_i(i1) = q_low(i2);
      id_low_i(i1) = id_low(i2);
      Vi(i1) = V(i2);
      if(i <= id_low(i2)){
        q_low_i(i1) -= H(i2,i);
        id_low_i(i1) -= 1;
      }
      qn_i(i1) = Qcumsum(i2,n-1) - H(i2,i);
      Hnew(0,i1) = H(i,i2);
      HnewT(i1,0) = H(i2,i);
    }
    Vi(n-1) = R_PosInf;
    id_low_i(n-1) = n-2;
    List ret_tmp = LCP_construction_path_distance(Vi, id_low_i, q_low_i, qn_i,Hnew, HnewT);
    arma::vec tmp1 = ret_tmp(0);
    for(int i1 = 0; i1 < n; i1++){
      alphas(i,i1) = tmp1(i1) * 1.0;
    }
  }
  List ret = List::create(Named("alphas") = alphas , _["vs"] = V);
  return ret;
}


