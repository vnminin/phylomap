
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <list>
#include <vector>
#include <sstream>

using namespace Rcpp;


NumericMatrix maketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N);
NumericMatrix SPARSEmaketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N);
NumericMatrix maketreelistEXP(List& x,NumericMatrix& Q,NumericVector& pid,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericMatrix& lefts,NumericMatrix& rights,NumericMatrix& d);


