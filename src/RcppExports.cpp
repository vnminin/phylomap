#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// maketreelistMCMC
NumericMatrix maketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N);
RcppExport SEXP phylomap_maketreelistMCMC(SEXP xSEXP,SEXP QSEXP,SEXP pidSEXP,SEXP BSEXP,SEXP OmegaSEXP,SEXP nenSEXP,SEXP nodelistSEXP,SEXP rootSEXP,SEXP NSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List& >::type x(xSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type Q(QSEXP);
        Rcpp::traits::input_parameter< NumericVector& >::type pid(pidSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type B(BSEXP);
        Rcpp::traits::input_parameter< double >::type Omega(OmegaSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nen(nenSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nodelist(nodelistSEXP);
        Rcpp::traits::input_parameter< int >::type root(rootSEXP);
        Rcpp::traits::input_parameter< int >::type N(NSEXP);

        NumericMatrix __result = maketreelistMCMC(x,Q,pid,B,Omega,nen,nodelist,root,N);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


// SPARSEmaketreelistMCMC
NumericMatrix SPARSEmaketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N);
RcppExport SEXP phylomap_SPARSEmaketreelistMCMC(SEXP xSEXP,SEXP QSEXP,SEXP pidSEXP,SEXP BSEXP,SEXP OmegaSEXP,SEXP nenSEXP,SEXP nodelistSEXP,SEXP rootSEXP,SEXP NSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List& >::type x(xSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type Q(QSEXP);
        Rcpp::traits::input_parameter< NumericVector& >::type pid(pidSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type B(BSEXP);
        Rcpp::traits::input_parameter< double >::type Omega(OmegaSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nen(nenSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nodelist(nodelistSEXP);
        Rcpp::traits::input_parameter< int >::type root(rootSEXP);
        Rcpp::traits::input_parameter< int >::type N(NSEXP);

        NumericMatrix __result = maketreelistMCMC(x,Q,pid,B,Omega,nen,nodelist,root,N);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

// maketreelistEXP
NumericMatrix maketreelistEXP(List& x,NumericMatrix& Q,NumericVector& pid,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericMatrix& lefts,NumericMatrix& rights,NumericMatrix& d);
RcppExport SEXP phylomap_maketreelistEXP(SEXP xSEXP,SEXP QSEXP,SEXP pidSEXP,SEXP nenSEXP,SEXP nodelistSEXP,SEXP rootSEXP,SEXP NSEXP,SEXP leftsSEXP,SEXP rightsSEXP,SEXP dSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List& >::type x(xSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type Q(QSEXP);
        Rcpp::traits::input_parameter< NumericVector& >::type pid(pidSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nen(nenSEXP);
        Rcpp::traits::input_parameter< IntegerVector& >::type nodelist(nodelistSEXP);
        Rcpp::traits::input_parameter< int >::type root(rootSEXP);
        Rcpp::traits::input_parameter< int >::type N(NSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type lefts(leftsSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type rights(rightsSEXP);
        Rcpp::traits::input_parameter< NumericMatrix& >::type d(dSEXP);

        NumericMatrix __result = maketreelistEXP(x,Q,pid,nen,nodelist,root,N,lefts,rights,d);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
