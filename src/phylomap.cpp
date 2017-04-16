#include "phylomap.h"


//
//
// This file contains the C code required to sample full state histories of a phylogenetic tree
// conditional on tip states in two different ways. The first way is a conditional MCMC sampler, the second
// uses a direct sampler that uses matrix exponentiation. 
//
//

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Branch functions /////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// the Branch structure contains a list of doubles corresponding to dwell times and a list of integers corresponding to states
struct Branch{
  std::list<double> branch;
  std::list<int> names; 
};

// makeabranch returns a Branch structure given a NumericVector of dwell times (branch), and an integer vector of states (names)
Branch makeabranch(NumericVector branch, IntegerVector names) {
  int i;
  std::list<double> branchlist;
  std::list<int> nameslist;
  for(i=1; i<=branch.size(); ++i) branchlist.push_back(branch(i-1));
  for(i=1; i<=names.size(); ++i) nameslist.push_back(names(i-1)-1);
  Branch lam; 
  lam.branch = branchlist;
  lam.names = nameslist;
  return lam;
}

// deletebranchjump takes a pointer to a Branch and removes the elements pointed to by the iterators, bit and nit
void deletebranchjump(Branch *lam, std::list<double>::iterator bit,std::list<int>::iterator nit) {
 (*lam).branch.erase(bit);
 (*lam).names.erase(nit);
 return;
}

// shortener removes self transitions from a Branch
void shortener(Branch *branch,arma::mat* dwelltimes,int ss,int iteration) {
  std::list<double>::iterator bit, bit2;
  bit = ((*branch).branch).begin();
  std::list<int>::iterator nit,nit2;
  nit = ((*branch).names).begin();
  int i;
  int n=((*branch).branch).size();
  for(i=0;i<(n-1);i++) {
    nit2=nit;bit2=bit;++nit2;++bit2;
    if((*nit)!=(*nit2)) {++nit;++bit;} else {
      (*bit)=(*bit)+(*bit2);
      deletebranchjump(branch,bit2,nit2);
    }
  }

  //
  nit = ((*branch).names).begin();
  nit2 = ((*branch).names).begin();
  ++nit2;
  n=((*branch).branch).size();
  for(i=1;i<n;i++) {
    if((*nit)<(*nit2)) (*dwelltimes)(iteration,ss+(*nit)*(ss-1)+(*nit2)-1)+=1;
    if((*nit)>(*nit2)) (*dwelltimes)(iteration,ss+(*nit)*(ss-1)+(*nit2))+=1;
    ++nit;
    ++nit2;
  }
  //
 
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Code for sampling branch states        //////////////////////////////////////////////////
////////////////////////////////////////// used in matrix exponentiation workflow //////////////////////////////////////////////////
//////////////////////////////////////////        2 functions                     //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int sampleOnce(arma::colvec weights, double rUnif) {
  double total = sum(weights);
  double cumProb=0;
  int i;
  for (i = 0; i < weights.n_elem; i++) {
    cumProb += weights(i) / total;
    if (rUnif < cumProb) break;
  }
  return(i);
}

// sample states along a branch conditional on starting and ending states
void newunifSample(int startState,int endState,double elapsedTime,double transProb,std::vector<Branch> *brancharray,int counter,arma::mat* dwelltimes,int iteration,int numStates,double poissonRate,arma::mat* B2) {

  RNGScope scope;

  // to avoid matrix*matrix calculations:
  arma::cube bpws = arma::cube(numStates,1,1);
  bpws.zeros();
  bpws(endState,0,0)=1;

  // simulate the number of jumps, including fictitious ones
  double rU = as<double>(runif(1));
  double cum = 0;  
 
  if (startState==endState){
    cum = ::Rf_dpois(0,poissonRate*elapsedTime,0)/transProb;
  }
 
  bool notExceed = true;
 
  if (cum > rU){ 
    notExceed = false; 
  }
 
  int numJumps = 0;
  double nextProb;  
  while(notExceed){
    numJumps++;
    if(numJumps>300) {
      Rcout << "newunifSample problem" << std::endl;
      List broken;
      broken["states"]=NumericVector::create(-1,-1);
      return;
    }
    bpws.insert_slices(numJumps,1,false);
    bpws.slice(numJumps) = (*B2)*bpws.slice(numJumps-1);
    nextProb = ::Rf_dpois(numJumps,poissonRate*elapsedTime,0)*bpws(startState,0,numJumps)/transProb;
    cum += nextProb;
    if (cum > rU){
      notExceed = false;
    }
  }
 
  List simPath = List::create(Named("states") = NULL, Named("times") = NULL);
 
  // if numJumps = 0: done
  if (numJumps == 0 || ((numJumps == 1) && (startState == endState))){
 
    simPath["states"] = NumericVector::create(startState, endState);
    simPath["times"] = NumericVector::create(0, elapsedTime);
 
  }else{ // if one true jump: done
    if ((numJumps == 1) && (startState != endState)){
 
      simPath["states"] = NumericVector::create(startState,endState,endState);
      simPath["times"] = NumericVector::create(0,elapsedTime*(as<double>(runif(1))),elapsedTime);
   }else{ // Case (nJmp >= 2) 
 
      // Simulate jumping times
      NumericVector dominJumpTimes = elapsedTime*runif(numJumps);
      std::sort(dominJumpTimes.begin(),dominJumpTimes.end());
      NumericVector dominStates(numJumps+1);
      dominStates[0] = startState;
      dominStates[numJumps] = endState;
 
      // Simulate states of the dominating Markov chain
      for (int i = 1; i < numJumps; i++){
        dominStates[i] = sampleOnce((*B2).row(dominStates[i-1]).t()%bpws.slice(numJumps-i), as<double>(runif(1)));
      }
 
      // Remove virtual substitutions
      std::vector<int> trueStates(1);
      std::vector<double> trueTimes(1);
      trueStates[0] = startState;
      trueTimes[0] = 0.0;
 
      for (int i = 1; i < dominStates.size(); i++){
        if (dominStates[i-1] != dominStates[i]){
          trueStates.push_back(dominStates[i]);
          trueTimes.push_back(dominJumpTimes[i-1]);
        }
      }
 
      trueStates.push_back(endState);
      trueTimes.push_back(elapsedTime);
 
      simPath["states"] = NumericVector(trueStates.begin(), trueStates.end());
      simPath["times"] = NumericVector(trueTimes.begin(), trueTimes.end());
    }
  }
  int i;
  IntegerVector states = as<IntegerVector>(simPath["states"]);
  arma::rowvec times = as<arma::rowvec>(simPath["times"]);

  IntegerVector fixstates(states.size()-1);
  for(i=0;i<fixstates.size();i++) fixstates(i)=states(i)+1;
  NumericVector fixtimes(times.size()-1);
  for(i=0;i<fixtimes.size();i++) fixtimes(i)=times(i+1)-times(i);
  Branch rm=makeabranch(fixtimes,fixstates);
  (*brancharray)[counter]=rm;

  // update dwell times
  std::list<int>::iterator nit,nit2;
  nit = ((rm).names).begin();
  nit2 = ((rm).names).begin();
  ++nit2;
  int n=((rm).branch).size();
  for(i=1;i<n;i++) {
    if((*nit)<(*nit2)) (*dwelltimes)(iteration,numStates+(*nit)*(numStates-1)+(*nit2)-1)+=1;
    if((*nit)>(*nit2)) (*dwelltimes)(iteration,numStates+(*nit)*(numStates-1)+(*nit2))+=1;
    ++nit;
    ++nit2;
  }
  //
 
  return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// code for resampling branches /////////////////////////////////////////////////////////
/////////////////////////// used in MCMC workflow        /////////////////////////////////////////////////////////
///////////////////////////       4 functions            /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// resamples states on a branch conditional on end points
void SPARSEresamplebranchstates(Branch *branch,arma::mat* B2,arma::sp_mat* B3) {

  // ss is the number of segments on the branch
  int ss=((*branch).branch).size();
  // if a branch only has 1 or 2 states there is nothing to be done
  if(ss==1){return;}
  if(ss==2){return;}
  // n is the total number of states
  int n = (*B2).n_rows;
  
  // vec starts off as a vector of zeros with a 1 in the position corresponding to the state at the end of the branch
  arma::colvec vec(n);
  vec.zeros();
  std::list<int>::iterator nit;
  nit = (*branch).names.end();
  --nit;
  vec(*nit)=1;

  // bpws is a matrix, column i will represent B^i %*% vec
  arma::mat bpws(n,ss);
  bpws.zeros();
  bpws.col(0)=vec;

  int i;
  int j;

  for(j=1;j<(ss-1);j++) bpws.col(j)=(*B3)*bpws.col(j-1); //SPARSE

  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  NumericVector prob(n);
  arma::mat top = arma::ones<arma::mat>(1,n);
  NumericVector tops(n);
  nit = (*branch).names.begin();
  for(i=1;i<(ss-1);i++) {
    top.row(0) = ((*B2).row(*nit)%bpws.col(ss-i-1).t()); //NOT SPARSE
    ++nit;
    for(j=0;j<n;j++) tops(j)=top(0,j);  
    *nit = as<int>(RcppArmadillo::sample(sts,1,1,tops));
  }

  return;
}

// resamplebranchstates resamples states given exisiting transitions and the states at either end
void resamplebranchstates(Branch *branch,arma::mat* B2) {

  // ss is the number of segments on the branch
  int ss=((*branch).branch).size();
  // if a branch only has 1 or 2 states there is nothing to be done
  if(ss==1){return;}
  if(ss==2){return;}
  // n is the total number of states
  int n = (*B2).n_rows;
  
  // vec starts off as a vector of zeros with a 1 in the position corresponding to the state at the end of the branch
  arma::colvec vec(n);
  vec.zeros();
  std::list<int>::iterator nit;
  nit = (*branch).names.end();
  --nit;
  vec(*nit)=1;

  // bpws is a matrix, column i will represent B^i %*% vec
  arma::mat bpws(n,ss);
  bpws.zeros();
  bpws.col(0)=vec;

  int i;
  int j;

  for(j=1;j<(ss-1);j++) bpws.col(j)=(*B2)*bpws.col(j-1);

 
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  NumericVector prob(n);
  arma::mat top = arma::ones<arma::mat>(1,n);
  NumericVector tops(n);
  nit = (*branch).names.begin();
  for(i=1;i<(ss-1);i++) {
    top.row(0) = ((*B2).row(*nit)%bpws.col(ss-i-1).t());
    ++nit;
    for(j=0;j<n;j++) tops(j)=top(0,j);  
    *nit = as<int>(RcppArmadillo::sample(sts,1,1,tops));
  }

  return;
}





// sampleabranch resamples states on a branch conditional on end point states and the locations of jumps
// it removes all the self transitions this creates
// and then resamples new self transitions

void SPARSEsampleabranch(Branch *branch,arma::mat* B2,arma::sp_mat* B3,double Omega,NumericMatrix* Q,arma::mat* dwelltimes,int statesize,int iteration) {
 RNGScope scope;

  std::list<double>::iterator bit;
  std::list<int>::iterator nit;

  // resample states given exisiting transitions
  SPARSEresamplebranchstates(branch,B2,B3); //SPARSE
  // remove self jumps
  shortener(branch,dwelltimes,statesize,iteration);

  int i;
  int k;
  int s;

  // n is the number of segments on the branch
  int n = ((*branch).names).size();
  bit = (*branch).branch.begin();
  nit = (*branch).names.begin();
  double segmentlength;
  double totallengthinserted;
  double rl;
  double r;
  for(i=0;i<n;i++) {
    segmentlength=(*bit);
    totallengthinserted=0;
    s=*nit;
    r=Omega+(*Q)(s,s);
    rl=0;
    while(totallengthinserted<segmentlength) {
     rl=as<double>(rexp(1,r));
     if((totallengthinserted+rl)<segmentlength) {
      (*branch).branch.insert(bit,rl);
      totallengthinserted+=rl;
      (*branch).names.insert(nit,s);
     } else {
      (*branch).branch.insert(bit,segmentlength-totallengthinserted);
      bit=(*branch).branch.erase(bit);
      ++nit;
      totallengthinserted=segmentlength;
     }
    }
  }


  return;
}

// sampleabranch first  resamples states on a branch given the exisiting transitions
//               second removes self transitions
//               third  samples new self transitions

void sampleabranch(Branch *branch,arma::mat* B2,double Omega,NumericMatrix* Q,arma::mat* dwelltimes,int statesize,int iteration) {
 RNGScope scope;

  std::list<double>::iterator bit;
  std::list<int>::iterator nit;

  resamplebranchstates(branch,B2);
  shortener(branch,dwelltimes,statesize,iteration);

  int i;
  int k;
  int s;
   
  // n is the number of segments on the branch
  int n = ((*branch).names).size();
  bit = (*branch).branch.begin();
  nit = (*branch).names.begin();
  double segmentlength;
  double totallengthinserted;
  double rl;
  double r;
  for(i=0;i<n;i++) {
    segmentlength=(*bit);
    totallengthinserted=0;
    s=*nit;
    r=Omega+(*Q)(s,s);
    rl=0;
    while(totallengthinserted<segmentlength) {
     rl=as<double>(rexp(1,r));
     if((totallengthinserted+rl)<segmentlength) {
      (*branch).branch.insert(bit,rl);
      totallengthinserted+=rl;
      (*branch).names.insert(nit,s);
     } else {
      (*branch).branch.insert(bit,segmentlength-totallengthinserted);
      bit=(*branch).branch.erase(bit);
      ++nit;
      totallengthinserted=segmentlength;
     }
    }
  }
  
  return;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// B is the original NumericMatrix provided
// B2 is the armadillo version
// B3 is the sparse armadillo version
// B4 is the transpose of B2


// Tvmmp returns vec %*% B^n with timing pointers
arma::vec Tvmmp(arma::mat* B4,arma::vec vec,int n) {
  if(n==0) return(vec);
  int i;
  for(i=0;i<n;i++) vec=(*B4)*vec;
  return(vec); 
}

// and the sparse version
arma::rowvec spvmmmm(arma::sp_mat* B,arma::rowvec vec,int n) {
  if(n==0) return(vec);
  int i;
  for(i=0;i<n;i++) vec=vec*(*B); //SPARSE 
  return(vec); 
}
// mmmmv returns B^n %*% vec
arma::vec mmmmvFORpl(arma::mat* B2,arma::vec* vec,int n) {
  if(n==0) return((*vec));
  for(int i=0;i<n;i++) (*vec)=(*B2)*(*vec);
  return((*vec)); 
}

// and the sparse version...
arma::vec spmmmmvFORpl(arma::sp_mat* B,arma::vec* vec,int n) {
  if(n==0) return((*vec));
  for(int i=0;i<n;i++) (*vec)=(*B)*(*vec); //SPARSE
  return((*vec)); 
}


void updatenodestates(std::vector<Branch> *brancharray,arma::imat* edge,arma::imat* nodestatesmatrix,arma::irowvec nodestates) {
  std::list<int>::iterator nit;
  int numedge = (*edge).n_rows;
  for(int i=0;i<numedge;i++) {
    //nodestatesmatrix starts the state space at 1
    (*nodestatesmatrix)(i,0)=nodestates((*edge)(i,0)-1);
    (*nodestatesmatrix)(i,1)=nodestates((*edge)(i,1)-1);
    //brancharray starts the state space at 0
    nit = ((*brancharray)[i]).names.begin();
    (*nit)=(nodestates((*edge)(i,0)-1)-1);
    nit = ((*brancharray)[i]).names.end();
    --nit;
    (*nit)=(nodestates((*edge)(i,1)-1)-1);
  }
  return;
}












// makePLrcpp creates the partial likelihood matrix

void SPARSEmakePLrcpp(arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::mat* PL,arma::irowvec* ne,arma::sp_mat* B3,arma::irowvec branchlengths){

 arma::vec first;
 arma::vec second;
  for(int i=0;i<Nnode;i++) {
    first=trans((*PL).row((*edge2)((*ne)(2*i+1)-1)-1));
    second=trans((*PL).row((*edge2)((*ne)(2*i)-1)-1));
    (*PL).row((*edge1)((*ne)(2*i)-1)-1)=trans(spmmmmvFORpl(B3,&first,branchlengths((*ne)(2*i+1)-1)-1)%spmmmmvFORpl(B3,&second,branchlengths((*ne)(2*i)-1)-1)); //SPARSE
  } 

 return;
}

void makePLrcpp(arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::mat* PL,arma::irowvec* ne,arma::mat* B2,arma::irowvec branchlengths){

 arma::vec first;
 arma::vec second;
 for(int i=0;i<Nnode;i++) {
   first=trans((*PL).row((*edge2)((*ne)(2*i+1)-1)-1));
   second=trans((*PL).row((*edge2)((*ne)(2*i)-1)-1));
   (*PL).row((*edge1)((*ne)(2*i)-1)-1)=trans(mmmmvFORpl(B2,&first,branchlengths((*ne)(2*i+1)-1)-1)%mmmmvFORpl(B2,&second,branchlengths((*ne)(2*i)-1)-1));
 }

 return;
}



// sampleinternalnodesSPARSE samples internal nodes conditional on tip states

arma::irowvec sampleinternalnodesSPARSE(std::vector<Branch> *brancharray,int branchnumber,arma::mat* PL,arma::rowvec* pid,arma::mat* B2,arma::sp_mat* B3,int root,IntegerVector nodelist,arma::irowvec* ne,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::irowvec* states) {

  RNGScope scope;
  int i;
  int j;

  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchnumber);
  for(i=0;i<branchnumber;i++) branchlengths(i)=((*brancharray)[i]).names.size();
  int n = (*B2).n_rows;
  int sss=(*states).size(); 
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  for(i=0;i<sss;i++) rm(i)=(*states)(i)-1;
  SPARSEmakePLrcpp(edge1,edge2,Nnode,PL,ne,B3,branchlengths); //SPARSE
  arma::rowvec covec = (*pid)%(*PL).row(root-1);
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));
  arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  NumericVector veccer(n);
  
  int cn;
  int pn;
  int ps;
  int ae;
  int nll = nodelist.size();
  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      cn=nodelist(i)-1; //child node
      j=0;
      while((*edge)(j,1)!=nodelist(i)) j++;
      pn = (*edge)(j,0)-1; //parent node
      ps = rm(pn); //parent state
      ae = j; //appropriate edge

      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = arma::zeros<arma::rowvec>(n);
      vecc(ps) = 1;
      vecc = spvmmmm(B3,vecc,branchlengths(ae)-1)%(*PL).row(cn); //SPARSE
      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }

  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;
  
  return rm;

}

// sampleinternalnodesMCMC samples internal nodes conditional on tip states

arma::irowvec sampleinternalnodesMCMC(std::vector<Branch> *brancharray,int branchnumber,arma::mat* PL,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,int root,IntegerVector nodelist,arma::irowvec* ne,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::irowvec* states) {

  RNGScope scope;
 
  int i;
  int j;

  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchnumber);
  for(i=0;i<branchnumber;i++) branchlengths(i)=((*brancharray)[i]).names.size();

  int n = (*B2).n_rows;

  // states: vector containing tip node states
  // sss: number of tips
  int sss=(*states).size(); 
  // edge: matrix containing parent labels (column 1) and child labels (column 2) for each edge (row)
  // arma::imat edge = as<arma::imat>(xtree["edge"]);

  // rm: vector to be returned, containing the states of all the nodes (tip nodes included)
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  // tip states remain unchanged
  for(i=0;i<sss;i++) rm(i)=(*states)(i)-1;
  // pid: the probability vector for states at the root

  // PL is the partial likelihood matrix, it the probability of the observed tip states under each node (row) assuming that node is in: state 1 (column 1), state 2 (column 2), ...
  makePLrcpp(edge1,edge2,Nnode,PL,ne,B2,branchlengths);
  // covec is proportional to the vector of probabilities, Pr(root|tip states)
  arma::rowvec covec = (*pid)%(*PL).row(root-1);
  // tops is a NumericVector version of covec
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  // sts is a vector of state names from 0 to n-1
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  // sample the root state, armadillo's sample function requires a NumericVector of probabilities
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));

  //arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  arma::vec vecc = arma::zeros<arma::vec>(n);
  NumericVector veccer(n);

  int cn;
  int pn;
  int ps;
  int ae;
  int nll = nodelist.size();
  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      cn=nodelist(i)-1; //child node
      j=0;
      while((*edge)(j,1)!=nodelist(i)) j++;
      pn = (*edge)(j,0)-1; //parent node
      ps = rm(pn); //parent state
      ae = j; //appropriate edge

      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans((*PL).row(cn));

      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }

  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;

  return rm;

}







void updatedwelltimes(int iteration,Branch* branch,arma::mat* dwelltimes) {
  std::list<int>::iterator nit;
  std::list<double>::iterator bit;
  nit=(*branch).names.begin();
  bit=(*branch).branch.begin();
  // i tracks which segment we are in
  for(int i=0;i<(*branch).names.size();i++) {
    (*dwelltimes)(iteration,*nit)+=*bit;
    ++nit;
    ++bit;
  }
  return;
}


// treesample samples a new tree conditional on the tip states and the number of jumps on each branch of the tree provided
void SPARSEtreesample(NumericMatrix* Q,arma::rowvec* pid,arma::mat* B2,arma::sp_mat* B3,double Omega,arma::irowvec* ne,std::vector<Branch> *brancharray,int branchnumber,arma::imat* nodestatesmatrix,arma::mat* dwelltimes,int root,arma::mat* PL,IntegerVector nodelist,int counter,int N,arma::irowvec* states,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode) {
RNGScope scope;

  int i;
  updatenodestates(brancharray,edge,nodestatesmatrix,sampleinternalnodesSPARSE(brancharray,branchnumber,PL,pid,B2,B3,root,nodelist,ne,edge,edge1,edge2,Nnode,states)); //SPARSE
  for(i=0;i<branchnumber;i++) SPARSEsampleabranch(&((*brancharray)[i]),B2,B3,Omega,Q,dwelltimes,(*Q).nrow(),counter);
  for(i=0;i<branchnumber;i++) updatedwelltimes(counter,&((*brancharray)[i]),dwelltimes);
 
  return;
}


// treesample samples a new tree conditional on the tip states and the number of jumps on each branch of the tree provided

void treesample(NumericMatrix* Q,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,double Omega,arma::irowvec* ne,std::vector<Branch> *brancharray,int branchnumber,arma::imat* nodestatesmatrix,arma::mat* dwelltimes,int root,arma::mat* PL,IntegerVector nodelist,int counter,int N,arma::irowvec* states,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode) {
RNGScope scope;
 
  int i;
  updatenodestates(brancharray,edge,nodestatesmatrix,sampleinternalnodesMCMC(brancharray,branchnumber,PL,pid,B2,B4,root,nodelist,ne,edge,edge1,edge2,Nnode,states));
  // for each branch, resample segmental states
  for(i=0;i<branchnumber;i++) sampleabranch(&((*brancharray)[i]),B2,Omega,Q,dwelltimes,(*Q).nrow(),counter);
  for(i=0;i<branchnumber;i++) updatedwelltimes(counter,&((*brancharray)[i]),dwelltimes);

  return;
}




// convert arma::mat to arma::sp_mat
arma::sp_mat matTospmat(arma::mat B) {
 
 int nrow=B.n_rows;
 int ncol=B.n_cols;
 arma::sp_mat A(nrow,ncol);
 int i;
 int j;
 for(i=0;i<nrow;i++) {
  for(j=0;j<ncol;j++) {
    //if(B(i,j)!=0) A(i,j)=B(i,j);
    if(B(i,j)>1e-7) A(i,j)=B(i,j);
  }
 }

 return A;
}


// SPARSEmaketreelistMCMC returns a list of (correlated) trees sampled iteratively via treesample

// [[Rcpp::export]]
NumericMatrix SPARSEmaketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N) {
 
  int i;
  RNGScope scope;


  List maps=x["maps"]; 
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchcount);
  for(i=0;i<maps.size();i++) branchlengths(i)=(as<arma::rowvec>(as<NumericVector>(as<SEXP>(maps[i])))).size();
  //Branch brancharray[maps.size()];
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  


  int Nnode = as<int>(x["Nnode"]); 
  arma::irowvec states=as<arma::irowvec>(x["states"]);
  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;

  // n is the size of the state space
  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::sp_mat B3=matTospmat(B2); //SPARSE
 
  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));
  
  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*(n-1));

  arma::rowvec rootdist =as<arma::rowvec>(pid);
 
  arma::irowvec ne=as<arma::irowvec>(nen);


  List w = x;
  List ret;
  for(i=0;i<N;i++){
    SPARSEtreesample(&Q,&rootdist,&B2,&B3,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);
    printf("%i \r",i);
  }
 
  return wrap(dwelltimes);
}




// maketreelistMCMC returns a list of (correlated) trees sampled iteratively via treesample
// x is a tree, in a tweaked version of ape's format
// Q is the rate matrix
// pid is the probability distribution over the states at the root node
// Omega is the rate of the dominating Poisson process    Omega>max(diag(Q))
// B is a probability transition matrix,   B=I+Q/Omega
// nen is a vector containing the order 
//    you should act on the branches if you want to start at the tips and 
//    work back to the root of the tree
// root is an integer corresponding to the root node
// nodelist is a vector containing the order
//    you should act on the nodes if you want to start below the root
//    and work down towards the tips of the tree
// N is an integer, it determines the number of state history samples returned

// [[Rcpp::export]]
NumericMatrix maketreelistMCMC(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N) {

  RNGScope scope;

  int i;
  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  //Branch brancharray[branchcount];
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

 arma::imat edge =  as<arma::imat>(x["edge"]);
 arma::irowvec edge1 = trans(edge.col(0));
 arma::irowvec edge2 = trans(edge.col(1));
 int Nnode = as<int>(x["Nnode"]);
 
  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);

  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;
 
  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::irowvec ne=as<arma::irowvec>(nen);

  List w = x;
  List ret;

  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*(n-1));

  for(i=0;i<N;i++){
    treesample(&Q,&rootdist,&B2,&B4,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);
    printf("%i \r",i);
  }

  return wrap(dwelltimes);

}





//////////////////////////////////////// Update Q matrix approach for 2 states ////////////////////////////////////////////////////////////





// shortener removes self transitions from a Branch
void shortenerbf(Branch *branch,arma::mat* dwelltimes,int ss,int iteration) {
  std::list<double>::iterator bit, bit2;
  bit = ((*branch).branch).begin();
  std::list<int>::iterator nit,nit2;
  nit = ((*branch).names).begin();
  int i;
  int n=((*branch).branch).size();

  //
  nit = ((*branch).names).begin();
  nit2 = ((*branch).names).begin();
  ++nit2;
  n=((*branch).branch).size();
  for(i=1;i<n;i++) {
    (*dwelltimes)(iteration,ss+(*nit)*ss+(*nit2))+=1;
    ++nit;
    ++nit2;
  }
  //

  nit = ((*branch).names).begin();
  bit = ((*branch).branch).begin();
  for(i=0;i<(n-1);i++) {
    nit2=nit;bit2=bit;++nit2;++bit2;
    if((*nit)!=(*nit2)) {++nit;++bit;} else {
      (*bit)=(*bit)+(*bit2);
      deletebranchjump(branch,bit2,nit2);
    }
  }

  return;
}


void sampleabranchbf(Branch *branch,arma::mat* B2,double Omega,NumericMatrix* Q,arma::mat* dwelltimes,int statesize,int iteration) {
 RNGScope scope;

  std::list<double>::iterator bit;
  std::list<int>::iterator nit;

  resamplebranchstates(branch,B2);
  shortenerbf(branch,dwelltimes,statesize,iteration);

  int i;
  int k;
  int s;
   
  // n is the number of segments on the branch
  int n = ((*branch).names).size();
  bit = (*branch).branch.begin();
  nit = (*branch).names.begin();
  double segmentlength;
  double totallengthinserted;
  double rl;
  double r;
  for(i=0;i<n;i++) {
    segmentlength=(*bit);
    totallengthinserted=0;
    s=*nit;
    r=Omega+(*Q)(s,s);
    rl=0;
    while(totallengthinserted<segmentlength) {
     rl=as<double>(rexp(1,r));
     if((totallengthinserted+rl)<segmentlength) {
      (*branch).branch.insert(bit,rl);
      totallengthinserted+=rl;
      (*branch).names.insert(nit,s);
     } else {
      (*branch).branch.insert(bit,segmentlength-totallengthinserted);
      bit=(*branch).branch.erase(bit);
      ++nit;
      totallengthinserted=segmentlength;
     }
    }
  }
  
  return;
}


void makePLnormalized(arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::mat* PL,arma::irowvec* ne,arma::mat* B2,arma::irowvec branchlengths){
 arma::vec first;
 arma::vec second;
 for(int i=0;i<Nnode;i++) {
   first=trans((*PL).row((*edge2)((*ne)(2*i+1)-1)-1));
   second=trans((*PL).row((*edge2)((*ne)(2*i)-1)-1));
   (*PL).row((*edge1)((*ne)(2*i)-1)-1)=trans(mmmmvFORpl(B2,&first,branchlengths((*ne)(2*i+1)-1)-1)%mmmmvFORpl(B2,&second,branchlengths((*ne)(2*i)-1)-1));
   (*PL).row((*edge1)((*ne)(2*i)-1)-1)=(*PL).row((*edge1)((*ne)(2*i)-1)-1)/sum((*PL).row((*edge1)((*ne)(2*i)-1)-1));
 }

 return;
}


arma::irowvec sampleinternalnodesMCMCbf(std::vector<Branch> *brancharray,int branchnumber,arma::mat* PL,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,int root,IntegerVector nodelist,arma::irowvec* ne,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::irowvec* states,arma::mat* dwelltimes,int iteration) {

  RNGScope scope;
 
  int i;
  int j;

  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchnumber);
  for(i=0;i<branchnumber;i++) branchlengths(i)=((*brancharray)[i]).names.size();

  int n = (*B2).n_rows;

  // states: vector containing tip node states
  // sss: number of tips
  int sss=(*states).size(); 
  // edge: matrix containing parent labels (column 1) and child labels (column 2) for each edge (row)
  // arma::imat edge = as<arma::imat>(xtree["edge"]);

  // rm: vector to be returned, containing the states of all the nodes (tip nodes included)
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  // tip states remain unchanged
  for(i=0;i<sss;i++) rm(i)=(*states)(i)-1;
  // pid: the probability vector for states at the root
  // PL is the partial likelihood matrix, it the probability of the observed tip states under each node (row) assuming that node is in: state 1 (column 1), state 2 (column 2), ...
  //makePLrcpp(edge1,edge2,Nnode,PL,ne,B2,branchlengths); 

  makePLnormalized(edge1,edge2,Nnode,PL,ne,B2,branchlengths); 

  // covec is proportional to the vector of probabilities, Pr(root|tip states)
  arma::rowvec covec = (*pid)%(*PL).row(root-1);
  // tops is a NumericVector version of covec
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  // sts is a vector of state names from 0 to n-1
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  // sample the root state, armadillo's sample function requires a NumericVector of probabilities
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));
  (*dwelltimes)(iteration,8)=rm(root-1);
  //arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  arma::vec vecc = arma::zeros<arma::vec>(n);
  NumericVector veccer(n);

  int cn;
  int pn;
  int ps;
  int ae;
  int nll = nodelist.size();
  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      cn=nodelist(i)-1; //child node
      j=0;
      while((*edge)(j,1)!=nodelist(i)) j++;
      pn = (*edge)(j,0)-1; //parent node
      ps = rm(pn); //parent state
      ae = j; //appropriate edge

      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans((*PL).row(cn));

      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }

  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;
  return rm;

}




void treesamplebf(NumericMatrix* Q,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,double Omega,arma::irowvec* ne,std::vector<Branch> *brancharray,int branchnumber,arma::imat* nodestatesmatrix,arma::mat* dwelltimes,int root,arma::mat* PL,IntegerVector nodelist,int iteration,int N,arma::irowvec* states,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode) {
RNGScope scope;

  int i;
  updatenodestates(brancharray,edge,nodestatesmatrix,sampleinternalnodesMCMCbf(brancharray,branchnumber,PL,pid,B2,B4,root,nodelist,ne,edge,edge1,edge2,Nnode,states,dwelltimes,iteration));
  // for each branch, resample segmental states
  for(i=0;i<branchnumber;i++) sampleabranchbf(&((*brancharray)[i]),B2,Omega,Q,dwelltimes,(*Q).nrow(),iteration);
  for(i=0;i<branchnumber;i++) updatedwelltimes(iteration,&((*brancharray)[i]),dwelltimes);

  return;
}

void recordQ(NumericMatrix* Q,arma::mat* dwelltimes,int iteration) {
  (*dwelltimes)(iteration,6)=(*Q)(0,1);
  (*dwelltimes)(iteration,7)=(*Q)(1,0);
  return;
}




void updatel01(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {

    int n00=(*dwelltimes)(iteration,2);
    int n01=(*dwelltimes)(iteration,3);
    double t0=(*dwelltimes)(iteration,0);
    double l01=(*Q)(0,1);
   
    int n10=(*dwelltimes)(iteration,4);
    int n11=(*dwelltimes)(iteration,5);
    double t1=(*dwelltimes)(iteration,1);
    double l10=(*Q)(1,0);

    double newl01=::Rf_rgamma(prior(0)+n01, 1/(prior(1)+t0));
    double newl10=l10;

    if(newl01>Omega) return;

    double accept=pow((Omega-newl01)/(Omega-l01),n00)*exp(t0*(newl01-l01));

    if(accept>1) accept=1;
    double compare= as<double>(runif(1));

    (*Q)(0,0)=-newl01;
    (*Q)(0,1)=newl01;
    (*B2)(0,0)=1-newl01/Omega;
    (*B2)(0,1)=newl01/Omega;
    (*B4)(0,0)=1-newl01/Omega;
    (*B4)(1,0)=newl01/Omega;

  return;
}

void updatel10(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {
  
    int n00=(*dwelltimes)(iteration,2);
    int n01=(*dwelltimes)(iteration,3);
    double t0=(*dwelltimes)(iteration,0);
    double l01=(*Q)(0,1);
   
    int n10=(*dwelltimes)(iteration,4);
    int n11=(*dwelltimes)(iteration,5);
    double t1=(*dwelltimes)(iteration,1);
    double l10=(*Q)(1,0);

    double newl01=l01;
    double newl10=::Rf_rgamma(prior(2)+n10, 1/(prior(3)+t1));

    if(newl10>Omega) return;

    double accept=pow((Omega-newl10)/(Omega-l10),n11)*exp(t1*(newl10-l10));

    if(accept>1) accept=1;
    double compare= as<double>(runif(1));


    (*Q)(1,0)=newl10;
    (*Q)(1,1)=-newl10;
    (*B2)(1,0)=newl10/Omega;
    (*B2)(1,1)=1-newl10/Omega;
    (*B4)(0,1)=newl10/Omega;
    (*B4)(1,1)=1-newl10/Omega;

  return;
}



// [[Rcpp::export]]
NumericMatrix maketreelistMCMCbf(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericVector& prior) {

  RNGScope scope;

  int i;
  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  //Branch brancharray[branchcount];
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

 arma::imat edge =  as<arma::imat>(x["edge"]);
 arma::irowvec edge1 = trans(edge.col(0));
 arma::irowvec edge2 = trans(edge.col(1));
 int Nnode = as<int>(x["Nnode"]);
 
  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);

  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;
 
  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::irowvec ne=as<arma::irowvec>(nen);

  List w = x;
  List ret;

  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*n+2+1);

  for(i=0;i<N;i++){
    recordQ(&Q,&dwelltimes,i);
    treesamplebf(&Q,&rootdist,&B2,&B4,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);
    updatel01(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    updatel10(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    printf("%i \r",i);
  }

  return wrap(dwelltimes);

}




//////////////////////////////////////// Update Q matrix approach for k states ////////////////////////////////////////////////////////////



arma::irowvec sampleinternalnodesMCMCks(std::vector<Branch> *brancharray,int branchnumber,arma::mat* PL,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,int root,IntegerVector nodelist,arma::irowvec* ne,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::irowvec* states,arma::mat* dwelltimes,int iteration) {
  RNGScope scope;
 
  int i;
  int j;

  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchnumber);
  for(i=0;i<branchnumber;i++) branchlengths(i)=((*brancharray)[i]).names.size();

  int n = (*B2).n_rows;

  // states: vector containing tip node states
  // sss: number of tips
  int sss=(*states).size(); 
  // edge: matrix containing parent labels (column 1) and child labels (column 2) for each edge (row)
  // arma::imat edge = as<arma::imat>(xtree["edge"]);

  // rm: vector to be returned, containing the states of all the nodes (tip nodes included)
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  // tip states remain unchanged
  for(i=0;i<sss;i++) rm(i)=(*states)(i)-1;
  // pid: the probability vector for states at the root
  // PL is the partial likelihood matrix, it the probability of the observed tip states under each node (row) assuming that node is in: state 1 (column 1), state 2 (column 2), ...
  makePLnormalized(edge1,edge2,Nnode,PL,ne,B2,branchlengths); 

  // covec is proportional to the vector of probabilities, Pr(root|tip states)
  arma::rowvec covec = (*pid)%(*PL).row(root-1);
  // tops is a NumericVector version of covec
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  // sts is a vector of state names from 0 to n-1
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  // sample the root state, armadillo's sample function requires a NumericVector of probabilities
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));
  int k=n/2-1;
  int um=n+n*n+2+3*k;
  (*dwelltimes)(iteration,um)=rm(root-1);
  //arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  arma::vec vecc = arma::zeros<arma::vec>(n);
  NumericVector veccer(n);

  int cn;
  int pn;
  int ps;
  int ae;
  int nll = nodelist.size();

  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      cn=nodelist(i)-1; //child node
      j=0;
      while((*edge)(j,1)!=nodelist(i)) j++;
      pn = (*edge)(j,0)-1; //parent node
      ps = rm(pn); //parent state
      ae = j; //appropriate edge

      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans((*PL).row(cn));

      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }

  // update tip node states - they don't actually remain the same anymore
  for(i=0;i<(*edge2).size();i++){
    if((*edge)(i,1)<=sss) {
      cn=(*edge)(i,1)-1;
      pn=(*edge)(i,0)-1;
      ps = rm(pn);
      ae = i;
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans((*PL).row(cn));
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }
  //

  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;
  return rm;

}

void updatenodestates4s(std::vector<Branch> *brancharray,arma::imat* edge,arma::imat* nodestatesmatrix,arma::irowvec nodestates) {
  std::list<int>::iterator nit;
  int numedge = (*edge).n_rows;
  for(int i=0;i<numedge;i++) {
    //nodestatesmatrix starts the state space at 1
    (*nodestatesmatrix)(i,0)=nodestates((*edge)(i,0)-1);
    (*nodestatesmatrix)(i,1)=nodestates((*edge)(i,1)-1);
    //brancharray starts the state space at 0
    nit = ((*brancharray)[i]).names.begin();
    (*nit)=(nodestates((*edge)(i,0)-1)-1);
    nit = ((*brancharray)[i]).names.end();
    --nit;
    (*nit)=(nodestates((*edge)(i,1)-1)-1);
  }
  return;
}

void treesampleks(NumericMatrix* Q,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,double Omega,arma::irowvec* ne,std::vector<Branch> *brancharray,int branchnumber,arma::imat* nodestatesmatrix,arma::mat* dwelltimes,int root,arma::mat* PL,IntegerVector nodelist,int iteration,int N,arma::irowvec* states,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode) {
RNGScope scope;

  int i;
  updatenodestates4s(brancharray,edge,nodestatesmatrix,sampleinternalnodesMCMCks(brancharray,branchnumber,PL,pid,B2,B4,root,nodelist,ne,edge,edge1,edge2,Nnode,states,dwelltimes,iteration));
  // for each branch, resample segmental states
  for(i=0;i<branchnumber;i++) sampleabranchbf(&((*brancharray)[i]),B2,Omega,Q,dwelltimes,(*Q).nrow(),iteration);
  for(i=0;i<branchnumber;i++) updatedwelltimes(iteration,&((*brancharray)[i]),dwelltimes);

  return;
}


void updateksl01(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {

  int i;
  int j;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*dwelltimes)(iteration,i);
  for(i=0;i<n;i++) for(j=0;j<n;j++) transitioncounts(i*n+j)=(*dwelltimes)(iteration,n+i*n+j);

  double alphaprime = prior(0);
  for(i=0;i<=k;i++) alphaprime=alphaprime+transitioncounts(2*i*n+2*i+1);
  double betaprime = prior(1);
  for(i=0;i<=k;i++) betaprime=betaprime+gammas(i)*sojourntimes(2*i);

    double newl01=::Rf_rgamma(alphaprime, 1/betaprime);

    double gammatimes=betaprime-prior(1);
    double logaccept = (newl01-lambdas(0))*gammatimes;
    logaccept=logaccept+transitioncounts(0)*log((Omega-rkappas(0)-gammas(0)*newl01)/(Omega-rkappas(0)-gammas(0)*lambdas(0)));
    for(i=1;i<k;i++) logaccept=logaccept+transitioncounts(2*i*n+2*i)*log((Omega-rkappas(i)-lkappas(i-1)-gammas(i)*newl01)/(Omega-rkappas(i)-lkappas(i-1)-gammas(i)*lambdas(0)));
    logaccept=logaccept+transitioncounts(2*k*n+2*k)*log((Omega-lkappas(k-1)-gammas(k)*newl01)/(Omega-lkappas(k-1)-gammas(k)*lambdas(0)));

    double compare= as<double>(runif(1));

    if(newl01+rkappas(0)>Omega) return;
    for(i=1;i<k;i++) if(gammas(i)*newl01+rkappas(i)+lkappas(i-1)>Omega) return;
    if(gammas(k)*newl01+lkappas(k-1)>Omega) return;
    // the following line is wrong but we run up against machine precision otherwise...
    if(newl01<1e-300) return;

    if(logaccept<log(compare)) return;

    (*Q)(0,0)=-rkappas(0)-gammas(0)*newl01;
    (*Q)(0,1)=gammas(0)*newl01;
    (*B2)(0,0)=1+(*Q)(0,0)/Omega;
    (*B2)(0,1)=(*Q)(0,1)/Omega;
    (*B4)(0,0)=1+(*Q)(0,0)/Omega;
    (*B4)(1,0)=(*Q)(0,1)/Omega;

    for(i=1;i<k;i++) {
      (*Q)(2*i,2*i)=-lkappas(i-1)-rkappas(i)-gammas(i)*newl01;
      (*Q)(2*i,2*i+1)=gammas(i)*newl01;
      (*B2)(2*i,2*i)=1+(*Q)(2*i,2*i)/Omega;
      (*B2)(2*i,2*i+1)=(*Q)(2*i,2*i+1)/Omega;
      (*B4)(2*i,2*i)=1+(*Q)(2*i,2*i)/Omega;
      (*B4)(2*i+1,2*i)=(*Q)(2*i,2*i+1)/Omega;
    }

    (*Q)(2*k,2*k)=-lkappas(k-1)-gammas(k)*newl01;
    (*Q)(2*k,2*k+1)=gammas(k)*newl01;
    (*B2)(2*k,2*k)=1+(*Q)(2*k,2*k)/Omega;
    (*B2)(2*k,2*k+1)=(*Q)(2*k,2*k+1)/Omega;
    (*B4)(2*k,2*k)=1+(*Q)(2*k,2*k)/Omega;
    (*B4)(2*k+1,2*k)=(*Q)(2*k,2*k+1)/Omega;

  return;
}



void updateksl10(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {

  int i;
  int j;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);


  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*dwelltimes)(iteration,i);
  for(i=0;i<n;i++) for(j=0;j<n;j++) transitioncounts(i*n+j)=(*dwelltimes)(iteration,n+i*n+j);

  double alphaprime = prior(0);
  for(i=0;i<=k;i++) alphaprime=alphaprime+transitioncounts((2*i+1)*n+2*i);
  double betaprime = prior(1);
  for(i=0;i<=k;i++) betaprime=betaprime+gammas(i)*sojourntimes(2*i+1);

    double newl10=::Rf_rgamma(alphaprime, 1/betaprime);

    double gammatimes=betaprime-prior(1);
    double logaccept = (newl10-lambdas(1))*gammatimes;
    logaccept=logaccept+transitioncounts(n+1)*log((Omega-rkappas(0)-gammas(0)*newl10)/(Omega-rkappas(0)-gammas(0)*lambdas(1)));
    for(i=1;i<k;i++) logaccept=logaccept+transitioncounts((2*i+1)*n+2*i+1)*log((Omega-rkappas(i)-lkappas(i-1)-gammas(i)*newl10)/(Omega-rkappas(i)-lkappas(i-1)-gammas(i)*lambdas(1)));
    logaccept=logaccept+transitioncounts((2*k+1)*n+2*k+1)*log((Omega-lkappas(k-1)-gammas(k)*newl10)/(Omega-lkappas(k-1)-gammas(k)*lambdas(1)));

    double compare= as<double>(runif(1));

    if(newl10+rkappas(0)>Omega) return;
    for(i=1;i<k;i++) if(gammas(i)*newl10+rkappas(i)+lkappas(i-1)>Omega) return;
    if(gammas(k)*newl10+lkappas(k-1)>Omega) return;
    // the following line is wrong but we run up against machine precision otherwise...
    if(newl10<1e-300) return;

    if(logaccept<log(compare)) return;

    (*Q)(1,1)=-rkappas(0)-gammas(0)*newl10;
    (*Q)(1,0)=gammas(0)*newl10;
    (*B2)(1,1)=1+(*Q)(1,1)/Omega;
    (*B2)(1,0)=(*Q)(1,0)/Omega;
    (*B4)(1,1)=1+(*Q)(1,1)/Omega;
    (*B4)(0,1)=(*Q)(1,0)/Omega;
    


    for(i=1;i<k;i++) {
      (*Q)(2*i+1,2*i+1)=-lkappas(i-1)-rkappas(i)-gammas(i)*newl10;
      (*Q)(2*i+1,2*i)=gammas(i)*newl10;
      (*B2)(2*i+1,2*i+1)=1+(*Q)(2*i+1,2*i+1)/Omega;
      (*B2)(2*i+1,2*i)=(*Q)(2*i+1,2*i)/Omega;
      (*B4)(2*i+1,2*i+1)=1+(*Q)(2*i+1,2*i+1)/Omega;
      (*B4)(2*i,2*i+1)=(*Q)(2*i+1,2*i)/Omega;
    }

    (*Q)(2*k+1,2*k+1)=-lkappas(k-1)-gammas(k)*newl10;
    (*Q)(2*k+1,2*k)=gammas(k)*newl10;
    (*B2)(2*k+1,2*k+1)=1+(*Q)(2*k+1,2*k+1)/Omega;
    (*B2)(2*k+1,2*k)=(*Q)(2*k+1,2*k)/Omega;
    (*B4)(2*k+1,2*k+1)=1+(*Q)(2*k+1,2*k+1)/Omega;
    (*B4)(2*k,2*k+1)=(*Q)(2*k+1,2*k)/Omega;



  return;
}

void updaterkappas(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*dwelltimes)(iteration,i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*dwelltimes)(iteration,n+i*n+j2);

  double alphaprime = prior(2)+transitioncounts((2*j)*n+2*j+2)+transitioncounts((2*j+1)*n+2*j+3);
  double betaprime = prior(3)+sojourntimes(2*j)+sojourntimes(2*j+1);

    double newrkappaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newrkappaj-rkappas(j))*(sojourntimes(2*j)+sojourntimes(2*j+1));
    if(j==0) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-newrkappaj-gammas(j)*lambdas(0))/(Omega-rkappas(j)-gammas(j)*lambdas(0)));
    if(j==0) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-newrkappaj-gammas(j)*lambdas(1))/(Omega-rkappas(j)-gammas(j)*lambdas(1)));
    if(j>0)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-newrkappaj-gammas(j)*lambdas(0))/(Omega-lkappas(j-1)-rkappas(j)-gammas(j)*lambdas(0)));
    if(j>0)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-newrkappaj-gammas(j)*lambdas(1))/(Omega-lkappas(j-1)-rkappas(j)-gammas(j)*lambdas(1)));

    double compare= as<double>(runif(1));

    if(j==0) if(newrkappaj+gammas(j)*lambdas(0)>Omega) return;
    if(j==0) if(newrkappaj+gammas(j)*lambdas(1)>Omega) return;
    if(j>0)  if(newrkappaj+gammas(j)*lambdas(0)+lkappas(j-1)>Omega) return;
    if(j>0)  if(newrkappaj+gammas(j)*lambdas(1)+lkappas(j-1)>Omega) return;
    // the following line is wrong but we run up against machine precision otherwise...
    if(newrkappaj<1e-300) return;

    if(logaccept<log(compare)) return;


    (*Q)(2*j,2*j+2)=newrkappaj;
    (*Q)(2*j+1,2*j+3)=newrkappaj;
    if(j==0) (*Q)(0,0)=-newrkappaj-gammas(j)*lambdas(0);
    if(j==0) (*Q)(1,1)=-newrkappaj-gammas(j)*lambdas(1);
    if(j>0)  (*Q)(2*j,2*j)=-newrkappaj-lkappas(j-1)-gammas(j)*lambdas(0);
    if(j>0)  (*Q)(2*j+1,2*j+1)=-newrkappaj-lkappas(j-1)-gammas(j)*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j+2)=(*Q)(2*j,2*j+2)/Omega;
    (*B2)(2*j+1,2*j+3)=(*Q)(2*j+1,2*j+3)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j+2,2*j)=(*Q)(2*j,2*j+2)/Omega;
    (*B4)(2*j+3,2*j+1)=(*Q)(2*j+1,2*j+3)/Omega;
  

  return;
}




void updatelkappas(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*dwelltimes)(iteration,i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*dwelltimes)(iteration,n+i*n+j2);

  double alphaprime = prior(2)+transitioncounts((2*j)*n+2*j-2)+transitioncounts((2*j+1)*n+2*j-1);
  double betaprime = prior(3)+sojourntimes(2*j)+sojourntimes(2*j+1);

    double newlkappaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newlkappaj-lkappas(j-1))*(sojourntimes(2*j)+sojourntimes(2*j+1));
    if(j==k) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-newlkappaj-gammas(j)*lambdas(0))/(Omega-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j==k) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-newlkappaj-gammas(j)*lambdas(1))/(Omega-lkappas(j-1)-gammas(j)*lambdas(1)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-rkappas(j)-newlkappaj-gammas(j)*lambdas(0))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-rkappas(j)-newlkappaj-gammas(j)*lambdas(1))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(1)));
 
    double compare= as<double>(runif(1));

    if(j==k) if(newlkappaj+gammas(j)*lambdas(0)>Omega) return;
    if(j==k) if(newlkappaj+gammas(j)*lambdas(1)>Omega) return;
    if(j<k)  if(newlkappaj+gammas(j)*lambdas(0)+rkappas(j)>Omega) return;
    if(j<k)  if(newlkappaj+gammas(j)*lambdas(1)+rkappas(j)>Omega) return;
    // the following line is wrong but we run up against machine precision otherwise...
    if(newlkappaj<1e-300) return;

    if(logaccept<log(compare)) return;

    (*Q)(2*j,2*j-2)=newlkappaj;
    (*Q)(2*j+1,2*j-1)=newlkappaj;
    if(j==k) (*Q)(2*j,2*j)=-newlkappaj-gammas(j)*lambdas(0);
    if(j==k) (*Q)(2*j+1,2*j+1)=-newlkappaj-gammas(j)*lambdas(1);
    if(j<k)  (*Q)(2*j,2*j)=-newlkappaj-rkappas(j)-gammas(j)*lambdas(0);
    if(j<k)  (*Q)(2*j+1,2*j+1)=-newlkappaj-rkappas(j)-gammas(j)*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j-2)=(*Q)(2*j,2*j-2)/Omega;
    (*B2)(2*j+1,2*j-1)=(*Q)(2*j+1,2*j-1)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j-2,2*j)=(*Q)(2*j,2*j-2)/Omega;
    (*B4)(2*j-1,2*j+1)=(*Q)(2*j+1,2*j-1)/Omega;
  

  return;
}



void updategammas(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*dwelltimes)(iteration,i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*dwelltimes)(iteration,n+i*n+j2);

  double alphaprime = prior(4)+transitioncounts((2*j)*n+2*j+1)+transitioncounts((2*j+1)*n+2*j);
  double betaprime = prior(5)+sojourntimes(2*j)*lambdas(0)+sojourntimes(2*j+1)*lambdas(1);

    double newgammaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newgammaj-gammas(j))*(sojourntimes(2*j)*lambdas(0)+sojourntimes(2*j+1)*lambdas(1));
    if(j==k) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-newgammaj*lambdas(0))/(Omega-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j==k) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-newgammaj*lambdas(1))/(Omega-lkappas(j-1)-gammas(j)*lambdas(1)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(0))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(1))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(1)));
 
    double compare= as<double>(runif(1));

    if(j==k) if(lkappas(j-1)+newgammaj*lambdas(0)>Omega) return;
    if(j==k) if(lkappas(j-1)+newgammaj*lambdas(1)>Omega) return;
    if(j<k)  if(lkappas(j-1)+newgammaj*lambdas(0)+rkappas(j)>Omega) return;
    if(j<k)  if(lkappas(j-1)+newgammaj*lambdas(1)+rkappas(j)>Omega) return;
    // the following line is wrong but we run up against machine precision otherwise...
    if(newgammaj<1e-300) return;

    if(logaccept<log(compare)) return;


    (*Q)(2*j,2*j+1)=newgammaj*lambdas(0);
    (*Q)(2*j+1,2*j)=newgammaj*lambdas(1);
    if(j==k) (*Q)(2*j,2*j)=-lkappas(j-1)-newgammaj*lambdas(0);
    if(j==k) (*Q)(2*j+1,2*j+1)=-lkappas(j-1)-newgammaj*lambdas(1);
    if(j<k)  (*Q)(2*j,2*j)=-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(0);
    if(j<k)  (*Q)(2*j+1,2*j+1)=-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j+1)=(*Q)(2*j,2*j+1)/Omega;
    (*B2)(2*j+1,2*j)=(*Q)(2*j+1,2*j)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j+1,2*j)=(*Q)(2*j,2*j+1)/Omega;
    (*B4)(2*j,2*j+1)=(*Q)(2*j+1,2*j)/Omega;


  return;
}


void recordQks(NumericMatrix* Q,arma::mat* dwelltimes,int iteration,int n,int k) {
  (*dwelltimes)(iteration,n+n*n)=(*Q)(0,1);
  (*dwelltimes)(iteration,n+n*n+1)=(*Q)(1,0);
  int i;
  for(i=0;i<k;i++) (*dwelltimes)(iteration,n+n*n+2+i)=(*Q)(2*i,2*i+2);
  for(i=0;i<k;i++) (*dwelltimes)(iteration,n+n*n+2+k+i)=(*Q)(2*i+2,2*i);
  for(i=0;i<k;i++)(*dwelltimes)(iteration,n+n*n+2+2*k+i)=(*Q)(2*(i+1),2*(i+1)+1)/(*Q)(0,1);

  return;
}


// [[Rcpp::export]]
NumericMatrix maketreelistMCMCks(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericVector& prior) {

  RNGScope scope;

  int i;
  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  //Branch brancharray[branchcount];
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

 arma::imat edge =  as<arma::imat>(x["edge"]);
 arma::irowvec edge1 = trans(edge.col(0));
 arma::irowvec edge2 = trans(edge.col(1));
 int Nnode = as<int>(x["Nnode"]);
 
 int k=Q.nrow()/2-1;
 arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
 lambdas(0)=Q(0,1);
 lambdas(1)=Q(1,0); 
 arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) rkappas(i)=Q(2*i,2*i+2);
 arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) lkappas(i)=Q(2*i+2,2*i);
 arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
 gammas(0)=1;
 for(i=1;i<=k;i++) gammas(i)=Q(2*i,2*i+1)/lambdas(0);

  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);

  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int j;
   for(i=0;i<states.size();i++) {
     if((states(i)-2*floor(states(i)/2))==0) {
       for(j=1;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
     if((states(i)-2*floor(states(i)/2))==1) {
       for(j=0;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
    }

  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::irowvec ne=as<arma::irowvec>(nen);

  List w = x;
  List ret;

  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*n+2+3*k+1);

  for(i=0;i<N;i++){
    recordQks(&Q,&dwelltimes,i,n,k);
    treesampleks(&Q,&rootdist,&B2,&B4,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);
    updateksl01(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    updateksl10(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    for(j=0;j<k;j++)  updaterkappas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    for(j=1;j<=k;j++) updatelkappas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    for(j=1;j<=k;j++) updategammas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    printf("%i \r",i);
  }

  return wrap(dwelltimes);

}


//////////////////////////////////////// Rate Matrix Updated ////////////////////////////////////////////////////////////
////////////////////////////////////////   Multiple Trees    ////////////////////////////////////////////////////////////
////////////////////////////////////////     2 state         ////////////////////////////////////////////////////////////


// shortener removes self transitions from a Branch
void shortenermtNS(Branch *branch,arma::rowvec* jodt,int ss) {

  std::list<double>::iterator bit, bit2;
  bit = ((*branch).branch).begin();
  std::list<int>::iterator nit,nit2;
  nit = ((*branch).names).begin();
  int i;
  int n=((*branch).branch).size();

  //
  nit = ((*branch).names).begin();
  nit2 = ((*branch).names).begin();
  ++nit2;
  n=((*branch).branch).size();
  for(i=1;i<n;i++) {
    (*jodt)(ss+(*nit)*ss+(*nit2))+=1;
    ++nit;
    ++nit2;
  }
  //

  nit = ((*branch).names).begin();
  bit = ((*branch).branch).begin();
  for(i=0;i<(n-1);i++) {
    nit2=nit;bit2=bit;++nit2;++bit2;
    if((*nit)!=(*nit2)) {++nit;++bit;} else {
      (*bit)=(*bit)+(*bit2);
      deletebranchjump(branch,bit2,nit2);
    }
  }

  return;
}




void updatenodestatesmt(arma::irowvec nodestates,int wt,std::vector< std::vector<Branch> >* BranchArray,std::vector< arma::imat >* Edge,std::vector< arma::imat >* NodeStatesMatrix) {
  std::list<int>::iterator nit;
  int numedge = ((*Edge)[wt]).n_rows;
  for(int i=0;i<numedge;i++) {
    //nodestatesmatrix starts the state space at 1
    ((*NodeStatesMatrix)[wt])(i,0)=nodestates(((*Edge)[wt])(i,0)-1);
    ((*NodeStatesMatrix)[wt])(i,1)=nodestates(((*Edge)[wt])(i,1)-1);
    //brancharray starts the state space at 0
    nit = (((*BranchArray)[wt])[i]).names.begin();
    (*nit)=(nodestates(((*Edge)[wt])(i,0)-1)-1);
    nit = (((*BranchArray)[wt])[i]).names.end();
    --nit;
    (*nit)=(nodestates(((*Edge)[wt])(i,1)-1)-1);
  }


  return;
}


void makePLrcppmt(int Nnode,arma::mat* B2,arma::irowvec branchlengths,int wt,std::vector< arma::irowvec >* Edge1,std::vector< arma::irowvec >* Edge2,std::vector< arma::mat >* PartialLikelihood,arma::imat* ne_m){
 arma::vec first;
 arma::vec second;
 for(int i=0;i<Nnode;i++) {
   first=trans(((*PartialLikelihood)[wt]).row(((*Edge2)[wt])(((*ne_m).row(wt))(2*i+1)-1)-1));
   second=trans(((*PartialLikelihood)[wt]).row(((*Edge2)[wt])(((*ne_m).row(wt))(2*i)-1)-1));
   ((*PartialLikelihood)[wt]).row(((*Edge1)[wt])(((*ne_m).row(wt))(2*i)-1)-1)=trans(mmmmvFORpl(B2,&first,branchlengths(((*ne_m).row(wt))(2*i+1)-1)-1)%mmmmvFORpl(B2,&second,branchlengths(((*ne_m).row(wt))(2*i)-1)-1));

 }

 return;
}


arma::irowvec sampleinternalnodesMCMCmt(int branchnumber,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,int root,int Nnode,int wt,std::vector< std::vector<Branch> >* BranchArray,std::vector< arma::mat >* PartialLikelihood,IntegerMatrix* nodelist_m,arma::imat* ne_m,std::vector< arma::imat >* Edge,std::vector< arma::irowvec >* Edge1,std::vector< arma::irowvec >* Edge2,std::vector< arma::irowvec >* States) {

  RNGScope scope;
 
  int i;
  int j;

  arma::irowvec branchlengths=arma::zeros<arma::irowvec>(branchnumber);
  for(i=0;i<branchnumber;i++) branchlengths(i)=(((*BranchArray)[wt])[i]).names.size();

  int n = (*B2).n_rows;

  // states: vector containing tip node states
  // sss: number of tips
  //int sss=(*states).size(); 
  int sss=((*States)[wt]).size(); 

  // edge: matrix containing parent labels (column 1) and child labels (column 2) for each edge (row)
  // arma::imat edge = as<arma::imat>(xtree["edge"]);

  // rm: vector to be returned, containing the states of all the nodes (tip nodes included)
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  // tip states remain unchanged
  for(i=0;i<sss;i++) rm(i)=((*States)[wt])(i)-1;
  // pid: the probability vector for states at the root
  // PL is the partial likelihood matrix, it the probability of the observed tip states under each node (row) assuming that node is in: state 1 (column 1), state 2 (column 2), ...
  makePLrcppmt(Nnode,B2,branchlengths,wt,Edge1,Edge2,PartialLikelihood,ne_m); 
  // covec is proportional to the vector of probabilities, Pr(root|tip states)
  arma::rowvec covec = (*pid)%((*PartialLikelihood)[wt]).row(root-1);
  // tops is a NumericVector version of covec
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  // sts is a vector of state names from 0 to n-1
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  // sample the root state, armadillo's sample function requires a NumericVector of probabilities
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));
  //arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  arma::vec vecc = arma::zeros<arma::vec>(n);
  NumericVector veccer(n);

  int cn;
  int pn;
  int ps;
  int ae;

  int nll = ((*nodelist_m).row(wt)).size();
  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      //cn=nodelist(i)-1; //child node
      cn=(*nodelist_m)(wt,i)-1; //child node

      j=0;
      //while((*edge)(j,1)!=nodelist(i)) j++;
      while(((*Edge)[wt])(j,1)!=(*nodelist_m)(wt,i)) j++;

      //pn = (*edge)(j,0)-1; //parent node
      pn = ((*Edge)[wt])(j,0)-1; //parent node

      ps = rm(pn); //parent state
      ae = j; //appropriate edge

      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans(((*PartialLikelihood)[wt]).row(cn));

      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }

  // update tip node states - they don't actually remain the same anymore
  for(i=0;i<((*Edge2)[wt]).size();i++){
    if(((*Edge)[wt])(i,1)<=sss) {
      cn=((*Edge)[wt])(i,1)-1;
      pn=((*Edge)[wt])(i,0)-1;
      ps = rm(pn);
      ae = i;
      vecc = arma::zeros<arma::vec>(n);
      vecc(ps) = 1;
      vecc = Tvmmp(B4,vecc,branchlengths(ae)-1)%trans(((*PartialLikelihood)[wt]).row(cn));
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }
  //

  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;
  return rm;

}


// resamplebranchstates resamples states given exisiting transitions and the states at either end
void resamplebranchstatesmt(Branch *branch,arma::mat* B2) {

  // ss is the number of segments on the branch
  int ss=((*branch).branch).size();
  // if a branch only has 1 or 2 states there is nothing to be done
  if(ss==1){return;}
  if(ss==2){return;}
  // n is the total number of states
  int n = (*B2).n_rows;
  
  // vec starts off as a vector of zeros with a 1 in the position corresponding to the state at the end of the branch
  arma::colvec vec(n);
  vec.zeros();
  std::list<int>::iterator nit;
  nit = (*branch).names.end();
  --nit;
  vec(*nit)=1;

  // bpws is a matrix, column i will represent B^i %*% vec
  arma::mat bpws(n,ss);
  bpws.zeros();
  bpws.col(0)=vec;

  int i;
  int j;

  for(j=1;j<(ss-1);j++) bpws.col(j)=(*B2)*bpws.col(j-1);

 
  IntegerVector sts(n);
  sts(0)=0;
  for(i=1;i<n;i++) sts(i)=sts(i-1)+1;
  NumericVector prob(n);
  arma::mat top = arma::ones<arma::mat>(1,n);
  NumericVector tops(n);
  nit = (*branch).names.begin();
  for(i=1;i<(ss-1);i++) {
    top.row(0) = ((*B2).row(*nit)%bpws.col(ss-i-1).t());
    ++nit;
    for(j=0;j<n;j++) tops(j)=top(0,j);  
    *nit = as<int>(RcppArmadillo::sample(sts,1,1,tops));
  }

  return;
}


void sampleabranchmtNS(Branch *branch,arma::mat* B2,double Omega,NumericMatrix* Q,arma::rowvec* jodt,int statesize) {
 RNGScope scope;

  std::list<double>::iterator bit;
  std::list<int>::iterator nit;

  resamplebranchstatesmt(branch,B2);
  shortenermtNS(branch,jodt,statesize);

  int i;
  int k;
  int s;
   
  // n is the number of segments on the branch
  int n = ((*branch).names).size();
  bit = (*branch).branch.begin();
  nit = (*branch).names.begin();
  double segmentlength;
  double totallengthinserted;
  double rl;
  double r;
  for(i=0;i<n;i++) {
    segmentlength=(*bit);
    totallengthinserted=0;
    s=*nit;
    r=Omega+(*Q)(s,s);
    rl=0;
    while(totallengthinserted<segmentlength) {
     rl=as<double>(rexp(1,r));
     if((totallengthinserted+rl)<segmentlength) {
      (*branch).branch.insert(bit,rl);
      totallengthinserted+=rl;
      (*branch).names.insert(nit,s);
     } else {
      (*branch).branch.insert(bit,segmentlength-totallengthinserted);
      bit=(*branch).branch.erase(bit);
      ++nit;
      totallengthinserted=segmentlength;
     }
    }
  }
  
  return;
}

void updatedwelltimesmtNS(Branch* branch,arma::rowvec* jodt) {
  std::list<int>::iterator nit;
  std::list<double>::iterator bit;
  nit=(*branch).names.begin();
  bit=(*branch).branch.begin();
  // i tracks which segment we are in
  for(int i=0;i<(*branch).names.size();i++) {
    (*jodt)(*nit)+=*bit;
    ++nit;
    ++bit;
  }
  return;
}


void treesamplemtNS(NumericMatrix* Q,arma::rowvec* pid,arma::mat* B2,arma::mat* B4,double Omega,int branchnumber,int root,int N,int Nnode,int wt, arma::imat* ne_m,std::vector< std::vector<Branch> >* BranchArray,std::vector< arma::imat >* NodeStatesMatrix,std::vector< arma::rowvec >* DwellTimes,std::vector< arma::mat >* PartialLikelihood,IntegerMatrix* nodelist_m,std::vector< arma::irowvec >* States,std::vector< arma::imat >* Edge,std::vector< arma::irowvec >* Edge1,std::vector< arma::irowvec >* Edge2) {
RNGScope scope;
 
  int i;
  updatenodestatesmt(sampleinternalnodesMCMCmt(branchnumber,pid,B2,B4,root,Nnode,wt,BranchArray,PartialLikelihood,nodelist_m,ne_m,Edge,Edge1,Edge2,States),wt,BranchArray,Edge,NodeStatesMatrix);
  // for each branch, resample segmental states
  for(i=0;i<branchnumber;i++) sampleabranchmtNS(&(((*BranchArray)[wt])[i]),B2,Omega,Q,&((*DwellTimes)[wt]),(*Q).nrow());
  for(i=0;i<branchnumber;i++) updatedwelltimesmtNS(&(((*BranchArray)[wt])[i]),&((*DwellTimes)[wt]));

  return;
}

void recordQmtNS(NumericMatrix* Q,arma::rowvec* jodt) {
  (*jodt)(6)=(*Q)(0,1);
  (*jodt)(7)=(*Q)(1,0);
  return;
}


double metropolis(double l01,double l10,int n00,int n01, int n10,int n11,double Omega,double t0,double t1,NumericVector prior,double newl01,double newl10) {
  double rm=1;
  rm=rm*exp(-Omega*(t0+t1))*pow((Omega-newl01),n00)*pow((Omega-newl10),n11)*pow(newl01,n01)*pow(newl10,n10);
  rm=rm*exp(::Rf_dgamma(newl01,prior(0), 1/prior(1), 1))*exp(::Rf_dgamma(newl10,prior(2), 1/prior(3), 1));
  rm=rm/(exp(-Omega*(t0+t1))*pow((Omega-l01),n00)*pow((Omega-l10),n11)*pow(l01,n01)*pow(l10,n10));
  rm=rm/(exp(::Rf_dgamma(l01,prior(0), 1/prior(1), 1))*exp(::Rf_dgamma(l10,prior(2), 1/prior(3), 1)));
  return(rm);
}

double hastings(double l01,double l10,int n00,int n01, int n10,int n11,double Omega,double t0,double t1,NumericVector prior,double newl01,double newl10) {
  double rm=1;
  rm=rm*exp(::Rf_dgamma(l01,prior(0)+n01, 1/(prior(1)+t0), 1))/exp(::Rf_dgamma(newl01,prior(0)+n01, 1/(prior(1)+t0), 1));
  rm=rm*exp(::Rf_dgamma(l10,prior(2)+n10, 1/(prior(3)+t1), 1))/exp(::Rf_dgamma(newl10,prior(2)+n10, 1/(prior(3)+t1), 1));
  return(rm);
}

void updatel01mtNS(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior) {

    int n00=(*jodt)(2);
    int n01=(*jodt)(3);
    double t0=(*jodt)(0);
    double l01=(*Q)(0,1);
   
    int n10=(*jodt)(4);
    int n11=(*jodt)(5);
    double t1=(*jodt)(1);
    double l10=(*Q)(1,0);

    double newl01=::Rf_rgamma(prior(0)+n01, 1/(prior(1)+t0));
    double newl10=l10;

    if(newl01>Omega) return;

    double accept=pow((Omega-newl01)/(Omega-l01),n00)*exp(t0*(newl01-l01));
    double acceptcompare=metropolis(l01,l10,n00,n01,n10,n11,Omega,t0,t1,prior,newl01,newl10)*hastings(l01,l10,n00,n01,n10,n11,Omega,t0,t1,prior,newl01,newl10);

    if(accept>1) accept=1;
    double compare= as<double>(runif(1));
  
    if(accept<compare) return;
    

    (*Q)(0,0)=-newl01;
    (*Q)(0,1)=newl01;
    (*B2)(0,0)=1-newl01/Omega;
    (*B2)(0,1)=newl01/Omega;
    (*B4)(0,0)=1-newl01/Omega;
    (*B4)(1,0)=newl01/Omega;

  return;
}

void updatel10mtNS(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior) {
  
    int n00=(*jodt)(2);
    int n01=(*jodt)(3);
    double t0=(*jodt)(0);
    double l01=(*Q)(0,1);
   
    int n10=(*jodt)(4);
    int n11=(*jodt)(5);
    double t1=(*jodt)(1);
    double l10=(*Q)(1,0);

    double newl01=l01;
    double newl10=::Rf_rgamma(prior(2)+n10, 1/(prior(3)+t1));

    if(newl10>Omega) return;

    double accept=pow((Omega-newl10)/(Omega-l10),n11)*exp(t1*(newl10-l10));
    double acceptcompare=metropolis(l01,l10,n00,n01,n10,n11,Omega,t0,t1,prior,newl01,newl10)*hastings(l01,l10,n00,n01,n10,n11,Omega,t0,t1,prior,newl01,newl10);
  
    if(accept>1) accept=1;
    double compare= as<double>(runif(1));

    if(accept<compare) return;
    

    (*Q)(1,0)=newl10;
    (*Q)(1,1)=-newl10;
    (*B2)(1,0)=newl10/Omega;
    (*B2)(1,1)=1-newl10/Omega;
    (*B4)(0,1)=newl10/Omega;
    (*B4)(1,1)=1-newl10/Omega;

  return;
}



// [[Rcpp::export]]
NumericMatrix maketreelistMCMCmt(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerMatrix& nen,IntegerMatrix& nodelist_m,IntegerVector roots,int N,NumericVector& prior) {

  RNGScope scope;

  int i;
  int j;
  int k;
  int treecount=x.size();


  List test=x[0];
  int Nnode = as<int>(test["Nnode"]);
  List maps=test["maps"];
  List mapnames=test["mapnames"];
  arma::imat edge =  as<arma::imat>(test["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));
  arma::imat nodestatesmatrix=as<arma::imat>(test["node.states"]);
  arma::irowvec states = as<arma::irowvec>(test["states"]);
  arma::mat PL((2*Nnode+1),Q.nrow());

  const int branchcount = maps.size();
  std::vector<Branch> brancharray(branchcount);

  std::vector< std::vector<Branch> > BranchArray(treecount);
  std::vector< arma::imat > Edge(treecount);
  std::vector< arma::irowvec > Edge1(treecount);
  std::vector< arma::irowvec > Edge2(treecount);
  std::vector< arma::imat > NodeStatesMatrix(treecount);
  std::vector< arma::irowvec > States(treecount);
  std::vector< arma::mat > PartialLikelihood(treecount);

  for(i=0;i<treecount;i++) {
   test=x[i];
   maps=test["maps"];
   mapnames=test["mapnames"];
   for(j=0;j<branchcount;j++) brancharray[j] = makeabranch(maps(j),mapnames(j));
   BranchArray[i]=brancharray;
   Edge[i]=as<arma::imat>(test["edge"]);
   Edge1[i]=trans((Edge[i]).col(0));
   Edge2[i]=trans((Edge[i]).col(1));
   nodestatesmatrix=as<arma::imat>(test["node.states"]);
   NodeStatesMatrix[i]=nodestatesmatrix;
   states = as<arma::irowvec>(test["states"]);
   States[i]=states;
   PL.zeros();
   for(j=0;j<states.size();j++) PL(j,states(j)-1) = 1;
   PartialLikelihood[i]=PL;
  }
 

  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::imat ne_m=as<arma::imat>(nen);

  List ret;

  arma::rowvec jodt(n+n*n+2);
  std::vector< arma::rowvec > DwellTimes(treecount);
  for(i=0;i<treecount;i++) DwellTimes[i]=jodt;

  // which tree should we use next
  double wt = as<double>(runif(1));
  arma::colvec weights(treecount);
  weights.ones();

  // to be returned matrix
  arma::mat tbr=arma::zeros<arma::mat>(N,n+n*n+2+1);

  for(i=0;i<N;i++){

    for(j=0;j<treecount;j++) {
      (DwellTimes[j]).zeros();
      recordQmtNS(&Q,&(DwellTimes[j]));
      treesamplemtNS(&Q,&rootdist,&B2,&B4,Omega,branchcount,roots[j],N,Nnode,j,&ne_m,&BranchArray,&NodeStatesMatrix,&DwellTimes,&PartialLikelihood,&nodelist_m,&States,&Edge,&Edge1,&Edge2);
    }

    wt = as<double>(runif(1));
    j=sampleOnce(weights,wt);

    for(k=0;k<n+n*n+2;k++) tbr(i,k)=(DwellTimes[j])(k);
    tbr(i,n+n*n+2)=j;

    updatel01mtNS(&Q,&B2,&B4,&(DwellTimes[j]),Omega,prior);
    updatel10mtNS(&Q,&B2,&B4,&(DwellTimes[j]),Omega,prior);

    printf("%i \r",i);
  }
		       				       

  return wrap(tbr);

}

//////////////////////////////////////// Rate Matrix Updated ////////////////////////////////////////////////////////////
////////////////////////////////////////   Multiple Trees    ////////////////////////////////////////////////////////////
////////////////////////////////////////     k state         ////////////////////////////////////////////////////////////




void updateksl01mt(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior) {

  int i;
  int j;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);
  
  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
    
  for(i=0;i<n;i++) sojourntimes(i)=(*jodt)(i);
  for(i=0;i<n;i++) for(j=0;j<n;j++) transitioncounts(i*n+j)=(*jodt)(n+i*n+j);
  
  double alphaprime = prior(0);
  for(i=0;i<=k;i++) alphaprime=alphaprime+transitioncounts(2*i*n+2*i+1);
  double betaprime = prior(1);
  for(i=0;i<=k;i++) betaprime=betaprime+gammas(i)*sojourntimes(2*i);
  
    double newl01=::Rf_rgamma(alphaprime, 1/betaprime);
  
    double gammatimes=betaprime-prior(1);
    double logaccept = (newl01-lambdas(0))*gammatimes;
    logaccept=logaccept+transitioncounts(0)*log((Omega-rkappas(0)-gammas(0)*newl01)/(Omega-rkappas(0)-gammas(0)*lambdas(0)));
    for(i=1;i<k;i++) logaccept=logaccept+transitioncounts(2*i*n+2*i)*log((Omega-rkappas(i)-lkappas(i-1)-gammas(i)*newl01)/(Omega-rkappas(i)-lkappas(i-1)-gammas(i)*lambdas(0)));
    logaccept=logaccept+transitioncounts(2*k*n+2*k)*log((Omega-lkappas(k-1)-gammas(k)*newl01)/(Omega-lkappas(k-1)-gammas(k)*lambdas(0)));
  
    double compare= as<double>(runif(1));

    if(newl01+rkappas(0)>Omega) return;
    for(i=1;i<k;i++) if(gammas(i)*newl01+rkappas(i)+lkappas(i-1)>Omega) return;
    if(gammas(k)*newl01+lkappas(k-1)>Omega) return;
 
    if(logaccept<log(compare)) return;

    (*Q)(0,0)=-rkappas(0)-gammas(0)*newl01;
    (*Q)(0,1)=gammas(0)*newl01;
    (*B2)(0,0)=1+(*Q)(0,0)/Omega;
    (*B2)(0,1)=(*Q)(0,1)/Omega;
    (*B4)(0,0)=1+(*Q)(0,0)/Omega;
    (*B4)(1,0)=(*Q)(0,1)/Omega;

    for(i=1;i<k;i++) {
      (*Q)(2*i,2*i)=-lkappas(i-1)-rkappas(i)-gammas(i)*newl01;
      (*Q)(2*i,2*i+1)=gammas(i)*newl01;
      (*B2)(2*i,2*i)=1+(*Q)(2*i,2*i)/Omega;
      (*B2)(2*i,2*i+1)=(*Q)(2*i,2*i+1)/Omega;
      (*B4)(2*i,2*i)=1+(*Q)(2*i,2*i)/Omega;
      (*B4)(2*i+1,2*i)=(*Q)(2*i,2*i+1)/Omega;
    }

    (*Q)(2*k,2*k)=-lkappas(k-1)-gammas(k)*newl01;
    (*Q)(2*k,2*k+1)=gammas(k)*newl01;
    (*B2)(2*k,2*k)=1+(*Q)(2*k,2*k)/Omega;
    (*B2)(2*k,2*k+1)=(*Q)(2*k,2*k+1)/Omega;
    (*B4)(2*k,2*k)=1+(*Q)(2*k,2*k)/Omega;
    (*B4)(2*k+1,2*k)=(*Q)(2*k,2*k+1)/Omega;
  
  return;
}



void updateksl10mt(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior) {

  int i;
  int j;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*jodt)(i);
  for(i=0;i<n;i++) for(j=0;j<n;j++) transitioncounts(i*n+j)=(*jodt)(n+i*n+j);

  double alphaprime = prior(2);
  for(i=0;i<=k;i++) alphaprime=alphaprime+transitioncounts((2*i+1)*n+2*i);
  double betaprime = prior(3);
  for(i=0;i<=k;i++) betaprime=betaprime+gammas(i)*sojourntimes(2*i+1);

    double newl10=::Rf_rgamma(alphaprime, 1/betaprime);
    double gammatimes=betaprime-prior(1);
    double logaccept = (newl10-lambdas(1))*gammatimes;
    logaccept=logaccept+transitioncounts(n+1)*log((Omega-rkappas(0)-gammas(0)*newl10)/(Omega-rkappas(0)-gammas(0)*lambdas(1)));
    for(i=1;i<k;i++) logaccept=logaccept+transitioncounts((2*i+1)*n+2*i+1)*log((Omega-rkappas(i)-lkappas(i-1)-gammas(i)*newl10)/(Omega-rkappas(i)-lkappas(i-1)-gammas(i)*lambdas(1)));
    logaccept=logaccept+transitioncounts((2*k+1)*n+2*k+1)*log((Omega-lkappas(k-1)-gammas(k)*newl10)/(Omega-lkappas(k-1)-gammas(k)*lambdas(1)));

    double compare= as<double>(runif(1));

    if(newl10+rkappas(0)>Omega) return;
    for(i=1;i<k;i++) if(gammas(i)*newl10+rkappas(i)+lkappas(i-1)>Omega) return;
    if(gammas(k)*newl10+lkappas(k-1)>Omega) return;

    if(logaccept<log(compare)) return;

    (*Q)(1,1)=-rkappas(0)-gammas(0)*newl10;
    (*Q)(1,0)=gammas(0)*newl10;
    (*B2)(1,1)=1+(*Q)(1,1)/Omega;
    (*B2)(1,0)=(*Q)(1,0)/Omega;
    (*B4)(1,1)=1+(*Q)(1,1)/Omega;
    (*B4)(0,1)=(*Q)(1,0)/Omega;

    for(i=1;i<k;i++) {
      (*Q)(2*i+1,2*i+1)=-lkappas(i-1)-rkappas(i)-gammas(i)*newl10;
      (*Q)(2*i+1,2*i)=gammas(i)*newl10;
      (*B2)(2*i+1,2*i+1)=1+(*Q)(2*i+1,2*i+1)/Omega;
      (*B2)(2*i+1,2*i)=(*Q)(2*i+1,2*i)/Omega;
      (*B4)(2*i+1,2*i+1)=1+(*Q)(2*i+1,2*i+1)/Omega;
      (*B4)(2*i,2*i+1)=(*Q)(2*i+1,2*i)/Omega;
    }

    (*Q)(2*k+1,2*k+1)=-lkappas(k-1)-gammas(k)*newl10;
    (*Q)(2*k+1,2*k)=gammas(k)*newl10;
    (*B2)(2*k+1,2*k+1)=1+(*Q)(2*k+1,2*k+1)/Omega;
    (*B2)(2*k+1,2*k)=(*Q)(2*k+1,2*k)/Omega;
    (*B4)(2*k+1,2*k+1)=1+(*Q)(2*k+1,2*k+1)/Omega;
    (*B4)(2*k,2*k+1)=(*Q)(2*k+1,2*k)/Omega;

  return;
}



void updaterkappasmt(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);
  
  for(i=0;i<n;i++) sojourntimes(i)=(*jodt)(i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*jodt)(n+i*n+j2);

  double alphaprime = prior(4)+transitioncounts((2*j)*n+2*j+2)+transitioncounts((2*j+1)*n+2*j+3);
  double betaprime = prior(5)+sojourntimes(2*j)+sojourntimes(2*j+1);

    double newrkappaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newrkappaj-rkappas(j))*(sojourntimes(2*j)+sojourntimes(2*j+1));
    if(j==0) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-newrkappaj-gammas(j)*lambdas(0))/(Omega-rkappas(j)-gammas(j)*lambdas(0)));
    if(j==0) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-newrkappaj-gammas(j)*lambdas(1))/(Omega-rkappas(j)-gammas(j)*lambdas(1)));
    if(j>0)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-newrkappaj-gammas(j)*lambdas(0))/(Omega-lkappas(j-1)-rkappas(j)-gammas(j)*lambdas(0)));
    if(j>0)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-newrkappaj-gammas(j)*lambdas(1))/(Omega-lkappas(j-1)-rkappas(j)-gammas(j)*lambdas(1)));

    double compare= as<double>(runif(1));

    if(j==0) if(newrkappaj+gammas(j)*lambdas(0)>Omega) return;
    if(j==0) if(newrkappaj+gammas(j)*lambdas(1)>Omega) return;
    if(j>0)  if(newrkappaj+gammas(j)*lambdas(0)+lkappas(j-1)>Omega) return;
    if(j>0)  if(newrkappaj+gammas(j)*lambdas(1)+lkappas(j-1)>Omega) return;

    if(logaccept<log(compare)) return;


    (*Q)(2*j,2*j+2)=newrkappaj;
    (*Q)(2*j+1,2*j+3)=newrkappaj;
    if(j==0) (*Q)(0,0)=-newrkappaj-gammas(j)*lambdas(0);
    if(j==0) (*Q)(1,1)=-newrkappaj-gammas(j)*lambdas(1);
    if(j>0)  (*Q)(2*j,2*j)=-newrkappaj-lkappas(j-1)-gammas(j)*lambdas(0);
    if(j>0)  (*Q)(2*j+1,2*j+1)=-newrkappaj-lkappas(j-1)-gammas(j)*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j+2)=(*Q)(2*j,2*j+2)/Omega;
    (*B2)(2*j+1,2*j+3)=(*Q)(2*j+1,2*j+3)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j+2,2*j)=(*Q)(2*j,2*j+2)/Omega;
    (*B4)(2*j+3,2*j+1)=(*Q)(2*j+1,2*j+3)/Omega;

  return;
}




void updatelkappasmt(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);

  for(i=0;i<n;i++) sojourntimes(i)=(*jodt)(i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*jodt)(n+i*n+j2);

  double alphaprime = prior(4)+transitioncounts((2*j)*n+2*j-2)+transitioncounts((2*j+1)*n+2*j-1);
  double betaprime = prior(5)+sojourntimes(2*j)+sojourntimes(2*j+1);

    double newlkappaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newlkappaj-lkappas(j-1))*(sojourntimes(2*j)+sojourntimes(2*j+1));
    if(j==k) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-newlkappaj-gammas(j)*lambdas(0))/(Omega-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j==k) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-newlkappaj-gammas(j)*lambdas(1))/(Omega-lkappas(j-1)-gammas(j)*lambdas(1)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-rkappas(j)-newlkappaj-gammas(j)*lambdas(0))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-rkappas(j)-newlkappaj-gammas(j)*lambdas(1))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(1)));
 
    double compare= as<double>(runif(1));

    if(j==k) if(newlkappaj+gammas(j)*lambdas(0)>Omega) return;
    if(j==k) if(newlkappaj+gammas(j)*lambdas(1)>Omega) return;
    if(j<k)  if(newlkappaj+gammas(j)*lambdas(0)+rkappas(j)>Omega) return;
    if(j<k)  if(newlkappaj+gammas(j)*lambdas(1)+rkappas(j)>Omega) return;
 
    if(logaccept<log(compare)) return;

    (*Q)(2*j,2*j-2)=newlkappaj;
    (*Q)(2*j+1,2*j-1)=newlkappaj;
    if(j==k) (*Q)(2*j,2*j)=-newlkappaj-gammas(j)*lambdas(0);
    if(j==k) (*Q)(2*j+1,2*j+1)=-newlkappaj-gammas(j)*lambdas(1);
    if(j<k)  (*Q)(2*j,2*j)=-newlkappaj-rkappas(j)-gammas(j)*lambdas(0);
    if(j<k)  (*Q)(2*j+1,2*j+1)=-newlkappaj-rkappas(j)-gammas(j)*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j-2)=(*Q)(2*j,2*j-2)/Omega;
    (*B2)(2*j+1,2*j-1)=(*Q)(2*j+1,2*j-1)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j-2,2*j)=(*Q)(2*j,2*j-2)/Omega;
    (*B4)(2*j-1,2*j+1)=(*Q)(2*j+1,2*j-1)/Omega;

  return;
}


void updategammasmt(NumericMatrix* Q,arma::mat* B2,arma::mat* B4,arma::rowvec* jodt,double Omega,NumericVector prior,int j) {

  int i;
  int j2;
  int n = (*Q).nrow();
  int k=n/2-1;
  arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
  lambdas(0)=(*Q)(0,1);
  lambdas(1)=(*Q)(1,0); 
  arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) rkappas(i)=(*Q)(2*i,2*i+2);
  arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
  for(i=0;i<k;i++) lkappas(i)=(*Q)(2*i+2,2*i);
  arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
  gammas(0)=1;
  for(i=1;i<=k;i++) gammas(i)=(*Q)(2*i,2*i+1)/lambdas(0);

  arma::rowvec sojourntimes = arma::zeros<arma::rowvec>(n);
  arma::rowvec transitioncounts = arma::zeros<arma::rowvec>(n*n);

  for(i=0;i<n;i++) sojourntimes(i)=(*jodt)(i);
  for(i=0;i<n;i++) for(j2=0;j2<n;j2++) transitioncounts(i*n+j2)=(*jodt)(n+i*n+j2);

  double alphaprime = prior(6)+transitioncounts((2*j)*n+2*j+1)+transitioncounts((2*j+1)*n+2*j);
  double betaprime = prior(7)+sojourntimes(2*j)*lambdas(0)+sojourntimes(2*j+1)*lambdas(1);

    double newgammaj=::Rf_rgamma(alphaprime, 1/betaprime);
    double logaccept = (newgammaj-gammas(j))*(sojourntimes(2*j)*lambdas(0)+sojourntimes(2*j+1)*lambdas(1));
    if(j==k) logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-newgammaj*lambdas(0))/(Omega-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j==k) logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-newgammaj*lambdas(1))/(Omega-lkappas(j-1)-gammas(j)*lambdas(1)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j)*n+2*j)*log((Omega-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(0))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(0)));
    if(j<k)  logaccept=logaccept+transitioncounts((2*j+1)*n+2*j+1)*log((Omega-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(1))/(Omega-rkappas(j)-lkappas(j-1)-gammas(j)*lambdas(1)));
 
    double compare= as<double>(runif(1));

    if(j==k) if(lkappas(j-1)+newgammaj*lambdas(0)>Omega) return;
    if(j==k) if(lkappas(j-1)+newgammaj*lambdas(1)>Omega) return;
    if(j<k)  if(lkappas(j-1)+newgammaj*lambdas(0)+rkappas(j)>Omega) return;
    if(j<k)  if(lkappas(j-1)+newgammaj*lambdas(1)+rkappas(j)>Omega) return;

    if(logaccept<log(compare)) return;


    (*Q)(2*j,2*j+1)=newgammaj*lambdas(0);
    (*Q)(2*j+1,2*j)=newgammaj*lambdas(1);
    if(j==k) (*Q)(2*j,2*j)=-lkappas(j-1)-newgammaj*lambdas(0);
    if(j==k) (*Q)(2*j+1,2*j+1)=-lkappas(j-1)-newgammaj*lambdas(1);
    if(j<k)  (*Q)(2*j,2*j)=-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(0);
    if(j<k)  (*Q)(2*j+1,2*j+1)=-lkappas(j-1)-rkappas(j)-newgammaj*lambdas(1);

    (*B2)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B2)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B2)(2*j,2*j+1)=(*Q)(2*j,2*j+1)/Omega;
    (*B2)(2*j+1,2*j)=(*Q)(2*j+1,2*j)/Omega;

    (*B4)(2*j,2*j)=1+(*Q)(2*j,2*j)/Omega;
    (*B4)(2*j+1,2*j+1)=1+(*Q)(2*j+1,2*j+1)/Omega;
    (*B4)(2*j+1,2*j)=(*Q)(2*j,2*j+1)/Omega;
    (*B4)(2*j,2*j+1)=(*Q)(2*j+1,2*j)/Omega;

  return;
}




void recordQksmt(NumericMatrix* Q,arma::rowvec* jodt,int n,int k) {
  (*jodt)(n+n*n)=(*Q)(0,1);
  (*jodt)(n+n*n+1)=(*Q)(1,0);
  int i;
  for(i=0;i<k;i++) (*jodt)(n+n*n+2+i)=(*Q)(2*i,2*i+2);
  for(i=0;i<k;i++) (*jodt)(n+n*n+2+k+i)=(*Q)(2*i+2,2*i);
  for(i=0;i<k;i++)(*jodt)(n+n*n+2+2*k+i)=(*Q)(2*(i+1),2*(i+1)+1)/(*Q)(0,1);

  return;
}


// [[Rcpp::export]]
NumericMatrix maketreelistMCMCksmt(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerMatrix& nen,IntegerMatrix& nodelist_m,IntegerVector roots,int N,NumericVector& prior) {

 
  RNGScope scope;

  int i;
  int j;
  int kk;
  int h;
  int RandomlyChoosenTree;
  int treecount=x.size();

  // there are k+1 'rates' of evolution (2k+2 states total) 
 int k=Q.nrow()/2-1;
 arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
 lambdas(0)=Q(0,1);
 lambdas(1)=Q(1,0); 
 arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) rkappas(i)=Q(2*i,2*i+2);
 arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) lkappas(i)=Q(2*i+2,2*i);
 arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
 gammas(0)=1;
 for(i=1;i<=k;i++) gammas(i)=Q(2*i,2*i+1)/lambdas(0);

  List test=x[0];
  int Nnode = as<int>(test["Nnode"]);
  List maps=test["maps"];
  List mapnames=test["mapnames"];
  arma::imat edge =  as<arma::imat>(test["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));
  arma::imat nodestatesmatrix=as<arma::imat>(test["node.states"]);
  arma::irowvec states = as<arma::irowvec>(test["states"]);
  arma::mat PL((2*Nnode+1),Q.nrow());

  const int branchcount = maps.size();
  std::vector<Branch> brancharray(branchcount);

  std::vector< std::vector<Branch> > BranchArray(treecount);
  std::vector< arma::imat > Edge(treecount);
  std::vector< arma::irowvec > Edge1(treecount);
  std::vector< arma::irowvec > Edge2(treecount);
  std::vector< arma::imat > NodeStatesMatrix(treecount);
  std::vector< arma::irowvec > States(treecount);
  std::vector< arma::mat > PartialLikelihood(treecount);

  for(i=0;i<treecount;i++) {
   test=x[i];
   maps=test["maps"];
   mapnames=test["mapnames"];
   for(j=0;j<branchcount;j++) brancharray[j] = makeabranch(maps(j),mapnames(j));
   BranchArray[i]=brancharray;
   Edge[i]=as<arma::imat>(test["edge"]);
   Edge1[i]=trans((Edge[i]).col(0));
   Edge2[i]=trans((Edge[i]).col(1));
   nodestatesmatrix=as<arma::imat>(test["node.states"]);
   NodeStatesMatrix[i]=nodestatesmatrix;
   states = as<arma::irowvec>(test["states"]);
   States[i]=states;
   PL.zeros();
   for(kk=0;kk<states.size();kk++) {
     if((states(kk)-2*floor(states(kk)/2))==0) {
       for(j=1;j<Q.nrow();j=j+2) PL(kk,j) = 1;
     }
     if((states(kk)-2*floor(states(kk)/2))==1) {
       for(j=0;j<Q.nrow();j=j+2) PL(kk,j) = 1;
     }
    }
   PartialLikelihood[i]=PL;
  }
 

  // n is the number of states, n=2k+2
  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::imat ne_m=as<arma::imat>(nen);

  List ret;

  arma::rowvec jodt(n+n*n+2+3*k);
  std::vector< arma::rowvec > DwellTimes(treecount);
  for(i=0;i<treecount;i++) DwellTimes[i]=jodt;

  // which tree should we use next
  double wt = as<double>(runif(1));
  arma::colvec weights(treecount);
  weights.ones();

  // to be returned matrix
  arma::mat tbr=arma::zeros<arma::mat>(N,n+n*n+2+3*k+1);
  for(i=0;i<N;i++){

    for(j=0;j<treecount;j++) {
      (DwellTimes[j]).zeros();
      recordQksmt(&Q,&(DwellTimes[j]),n,k);
      treesamplemtNS(&Q,&rootdist,&B2,&B4,Omega,branchcount,roots[j],N,Nnode,j,&ne_m,&BranchArray,&NodeStatesMatrix,&DwellTimes,&PartialLikelihood,&nodelist_m,&States,&Edge,&Edge1,&Edge2);
    }


    wt = as<double>(runif(1));
    RandomlyChoosenTree=sampleOnce(weights,wt);
 
    for(kk=0;kk<(n+n*n+2+3*k);kk++) tbr(i,kk)=(DwellTimes[RandomlyChoosenTree])(kk);
    tbr(i,n+n*n+2+3*k)=RandomlyChoosenTree;

    updateksl01mt(&Q,&B2,&B4,&(DwellTimes[RandomlyChoosenTree]),Omega,prior);
    updateksl10mt(&Q,&B2,&B4,&(DwellTimes[RandomlyChoosenTree]),Omega,prior);
    for(h=0;h<k;h++)  updaterkappasmt(&Q,&B2,&B4,&(DwellTimes[RandomlyChoosenTree]),Omega,prior,h);
    for(h=1;h<=k;h++) updatelkappasmt(&Q,&B2,&B4,&(DwellTimes[RandomlyChoosenTree]),Omega,prior,h);
    for(h=1;h<=k;h++) updategammasmt(&Q,&B2,&B4,&(DwellTimes[RandomlyChoosenTree]),Omega,prior,h);
    

    printf("%i \r",i);
  }
       				       

  return wrap(tbr);

}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////














//////////////////////////////////////// Matrix Exponentiation approach below ////////////////////////////////////////////////////////////








// makePLold creates the partial likelihood matrix for independent tree samples
arma::mat makePLold(List x,arma::irowvec ne,NumericMatrix Q,arma::cube TransProb){
  int Nnode = as<int>(x["Nnode"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);
  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int i;
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;

  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));

  int n = as<int>(x["Nnode"]);
  // for each node (row) calculate the probability of the tip states below that node conditioned on the node's state (column)
  for(i=0;i<n;i++) {
    PL.row(edge1(ne(2*i)-1)-1)=((TransProb.slice(ne(2*i)-1)*(PL.row(edge2(ne(2*i)-1)-1)).t())%(TransProb.slice(ne(2*i+1)-1)*(PL.row(edge2(ne(2*i+1)-1)-1)).t())).t();
  }
  return PL;
}


// makePLexp creates the partial likelihood matrix for independent tree samples
void makePLexp(arma::irowvec* states,int Nnode,arma::irowvec* edge1,arma::irowvec* edge2,IntegerVector* ne,arma::cube* TransProb,arma::mat* PL){
  int i;
  // for each node (row) calculate the probability of the tip states below that node conditioned on the node's state (column)
  for(i=0;i<Nnode;i++) {
    (*PL).row((*edge1)((*ne)(2*i)-1)-1)=(((*TransProb).slice((*ne)(2*i)-1)*((*PL).row((*edge2)((*ne)(2*i)-1)-1)).t())%((*TransProb).slice((*ne)(2*i+1)-1)*((*PL).row((*edge2)((*ne)(2*i+1)-1)-1)).t())).t();
  }
  return;
}


// sampleinternalnodesEXP samples internal nodes conditional on tip states
arma::irowvec sampleinternalnodesEXP(arma::imat* edge,arma::irowvec* states,arma::mat* PL,arma::rowvec* pid,int root,IntegerVector* nodelist,arma::cube* TransProb) {

  int i;
  // n is the number of states
  int n = (*pid).n_elem;
  // states: a vector containing tip states
  // sss: number of tips
  int sss=(*states).size();
  // edge: matrix containing parent labels (column 1) and child labels (column 2) for each edge (row)
  // rm: vector to be returned, containing the states of all the nodes (tip nodes included)
  arma::irowvec rm=arma::zeros<arma::irowvec>(2*sss-1);
  int j;
  // tip node states should remain the same
  for(i=0;i<sss;i++) rm(i)=(*states)(i)-1;
  // pid: the probability vector for states at the root
  // covec is proportional to the vector of probabilities, Pr(root|tip states)
  arma::rowvec covec = (*pid)%(*PL).row(root-1);
  // tops is a NumericVector version of covec
  NumericVector tops(n);
  for(i=0;i<n;i++) tops(i)= covec(i);
  // sts is a vector of state names from 0 to n-1
  IntegerVector sts(n);
  for(i=0;i<n;i++) sts(i)=i;
  // sample the root state, armadillo's sample function requires a NumericVector of probabilities
  rm(root-1)=as<int>(RcppArmadillo::sample(sts,1,1,tops));
  arma::rowvec vecc = arma::zeros<arma::rowvec>(n);
  NumericVector veccer(n);
  int cn;
  int pn;
  int ps;
  int ae;
  int nll = (*nodelist).size();
  if(nll>0) {  
    // sample a state for each (non-root) internal node starting below the root and working down
    for(i=0;i<nll;i++) {
      cn=(*nodelist)(i)-1; //child node
      j=0;
      while((*edge)(j,1)!=(*nodelist)(i)) j++;
      pn = (*edge)(j,0)-1; //parent node
      ps = rm(pn); //parent state
      ae = j; //appropriate edge
      vecc = arma::zeros<arma::rowvec>(n);
      // (probability of transitioning from parent state to child state) * (probability of tips below child | child state)
      vecc = ((*TransProb).slice(ae)).row(ps)%(*PL).row(cn);
      // armadillo sample needs a NumericVector
      for(j=0;j<n;j++) veccer(j)=vecc(j);
      rm(cn)=as<int>(RcppArmadillo::sample(sts,1,1,veccer));
    }
  }
  for(i=0;i<(2*sss-1);i++) rm(i)=rm(i)+1;
  return rm;
}

// matexp returns e^(Q*edgelength) where Q=left*D*right
arma::mat matexp(arma::mat left,arma::mat right,arma::mat D,double edgelength) {
  int i;
  for(i=0;i<D.n_rows;i++) D(i,i)=exp(D(i,i)*edgelength);
  return left*D*right;
}







// treesampleEXP samples a tree conditional on tip states, re-exponentiates Q
void treesampleEXP(std::vector<Branch> *brancharray,int branchnumber,arma::imat* nodestatesmatrix,arma::mat* left,arma::mat* right,arma::mat* D,arma::mat* dwelltimes,arma::irowvec* states,arma::imat* edge,arma::irowvec* edge1,arma::irowvec* edge2,int Nnode,arma::rowvec* pid,IntegerVector* ne,IntegerVector* nodelist,int iteration,int root,NumericVector edgelengths,arma::cube* TransProb,arma::mat* PL,int n,double poissonRate,arma::mat* B2) {
  int i;

  for(i=0;i<branchnumber;i++) (*TransProb).slice(i) = abs(matexp((*left),(*right),(*D),edgelengths(i))); //abs is wrong, it is a workaround for elements near zero
  makePLexp(states,Nnode,edge1,edge2,ne,TransProb,PL);
 
  updatenodestates(brancharray,edge,nodestatesmatrix,sampleinternalnodesEXP(edge,states,PL,pid,root,nodelist,TransProb));

  int startState;
  int endState;

  for(i=0;i<branchnumber;i++) {
    startState = (*nodestatesmatrix)(i,0)-1;
    endState = (*nodestatesmatrix)(i,1)-1;
    newunifSample(startState,endState,edgelengths(i),((*TransProb).slice(i))(startState,endState),brancharray,i,dwelltimes,iteration,n,poissonRate,B2);
    updatedwelltimes(iteration,&((*brancharray)[i]),dwelltimes);
  }

  return;
}



// [[Rcpp::export]]
NumericMatrix maketreelistEXP(List& x,NumericMatrix& Q,NumericVector& pid,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericMatrix& lefts,NumericMatrix& rights,NumericMatrix& d) {

  RNGScope scope;

  int i;
  int n = Q.nrow();
  arma::mat Q2(Q.begin(),n,n,false);
  double poissonRate = -1.0*arma::min(Q2.diag());
  arma::mat B2(n,n);
  B2.eye();
  B2=B2+Q2/poissonRate;

  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  //Branch brancharray[branchcount];
  std::vector<Branch> brancharray(branchcount);

  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));
  int Nnode = as<int>(x["Nnode"]);

  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);
  

  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*(n-1));

  arma::rowvec rootdist =as<arma::rowvec>(pid);
  NumericVector edgelengths = x["edge.length"];
  arma::mat left = as<arma::mat>(lefts);  
  arma::mat right = as<arma::mat>(rights); 
  arma::mat D = as<arma::mat>(d); 

  arma::cube TransProb = arma::cube(Q.nrow(),Q.nrow(),branchnumber);
  TransProb.zeros();
  
  for(i=0;i<branchnumber;i++) TransProb.slice(i) = abs(matexp(left,right,D,edgelengths(i))); //abs is wrong, it is a janky workaround for elements near zero
  arma::mat PL=makePLold(x,nen,Q,TransProb);
 
  for(i=0;i<N;i++){
    treesampleEXP(&brancharray,branchnumber,&nodestatesmatrix,&left,&right,&D,&dwelltimes,&states,&edge,&edge1,&edge2,Nnode,&rootdist,&nen,&nodelist,i,root,edgelengths,&TransProb,&PL,n,poissonRate,&B2);
    printf("%i \r",i);
  }

  return wrap(dwelltimes);
}







///////////////////////////////////////////////////////////////
/////////////////////// DIC CODE //////////////////////////////
///////////////////////////////////////////////////////////////


void PPupdatel01(NumericMatrix* Q,arma::mat* Q2,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {

    int n00=(*dwelltimes)(iteration,2);
    int n01=(*dwelltimes)(iteration,3);
    double t0=(*dwelltimes)(iteration,0);
    double l01=(*Q)(0,1);
   
    int n10=(*dwelltimes)(iteration,4);
    int n11=(*dwelltimes)(iteration,5);
    double t1=(*dwelltimes)(iteration,1);
    double l10=(*Q)(1,0);

    double newl01=::Rf_rgamma(prior(0)+n01, 1/(prior(1)+t0));
    double newl10=l10;

    if(newl01>Omega) return;

    double accept=pow((Omega-newl01)/(Omega-l01),n00)*exp(t0*(newl01-l01));

    if(accept>1) accept=1;
    double compare= as<double>(runif(1));

    (*Q)(0,0)=-newl01;
    (*Q)(0,1)=newl01;
    (*Q2)(0,0)=-newl01;
    (*Q2)(0,1)=newl01;
    (*B2)(0,0)=1-newl01/Omega;
    (*B2)(0,1)=newl01/Omega;
    (*B4)(0,0)=1-newl01/Omega;
    (*B4)(1,0)=newl01/Omega;

  return;
}

void PPupdatel10(NumericMatrix* Q,arma::mat* Q2,arma::mat* B2,arma::mat* B4,arma::mat* dwelltimes,int iteration,double Omega,NumericVector prior) {
  
    int n00=(*dwelltimes)(iteration,2);
    int n01=(*dwelltimes)(iteration,3);
    double t0=(*dwelltimes)(iteration,0);
    double l01=(*Q)(0,1);
   
    int n10=(*dwelltimes)(iteration,4);
    int n11=(*dwelltimes)(iteration,5);
    double t1=(*dwelltimes)(iteration,1);
    double l10=(*Q)(1,0);

    double newl01=l01;
    double newl10=::Rf_rgamma(prior(2)+n10, 1/(prior(3)+t1));

    if(newl10>Omega) return;

    double accept=pow((Omega-newl10)/(Omega-l10),n11)*exp(t1*(newl10-l10));

    if(accept>1) accept=1;
    double compare= as<double>(runif(1));


    (*Q)(1,0)=newl10;
    (*Q)(1,1)=-newl10;
    (*Q2)(1,0)=newl10;
    (*Q2)(1,1)=-newl10;
    (*B2)(1,0)=newl10/Omega;
    (*B2)(1,1)=1-newl10/Omega;
    (*B4)(0,1)=newl10/Omega;
    (*B4)(1,1)=1-newl10/Omega;

  return;
}


// PPmakePL creates the partial likelihood matrix for independent tree samples
arma::mat PPmakePL(List x,arma::irowvec states,arma::irowvec ne,NumericMatrix Q,arma::cube TransProb){
  int Nnode = as<int>(x["Nnode"]);
  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int i;
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;

  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));

  int n = as<int>(x["Nnode"]);
  // for each node (row) calculate the probability of the tip states below that node conditioned on the node's state (column)
  for(i=0;i<n;i++) {
    PL.row(edge1(ne(2*i)-1)-1)=((TransProb.slice(ne(2*i)-1)*(PL.row(edge2(ne(2*i)-1)-1)).t())%(TransProb.slice(ne(2*i+1)-1)*(PL.row(edge2(ne(2*i+1)-1)-1)).t())).t();
    // normalized for big trees...
    PL.row(edge1(ne(2*i)-1)-1)=PL.row(edge1(ne(2*i)-1)-1)/sum(PL.row(edge1(ne(2*i)-1)-1));
  }
  return PL;
}


// PPmakePLD creates the partial likelihood matrix for independent tree samples and keeps track of the renormalizing factors
arma::mat PPmakePLD(List x,arma::irowvec states,arma::irowvec ne,NumericMatrix Q,arma::cube TransProb,double *S){
  int Nnode = as<int>(x["Nnode"]);
  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int i;
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;

  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));

  int n = as<int>(x["Nnode"]);
  // for each node (row) calculate the probability of the tip states below that node conditioned on the node's state (column)
  for(i=0;i<n;i++) {
    PL.row(edge1(ne(2*i)-1)-1)=((TransProb.slice(ne(2*i)-1)*(PL.row(edge2(ne(2*i)-1)-1)).t())%(TransProb.slice(ne(2*i+1)-1)*(PL.row(edge2(ne(2*i+1)-1)-1)).t())).t();
    // normalized for big trees...
    (*S)=(*S)+log(sum(PL.row(edge1(ne(2*i)-1)-1)));
    PL.row(edge1(ne(2*i)-1)-1)=PL.row(edge1(ne(2*i)-1)-1)/sum(PL.row(edge1(ne(2*i)-1)-1));
  }
  return PL;
}



// [[Rcpp::export]]
NumericMatrix maketreelistMCMC2sDICt(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericVector& prior) {

  RNGScope scope;

  int i;
  int j;
  int f;
  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  //Branch brancharray[branchcount];
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

 arma::imat edge =  as<arma::imat>(x["edge"]);
 arma::irowvec edge1 = trans(edge.col(0));
 arma::irowvec edge2 = trans(edge.col(1));
 int Nnode = as<int>(x["Nnode"]);
 
  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);

  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  for(i=0;i<states.size();i++) PL(i,states(i)-1) = 1;
 
  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::irowvec ne=as<arma::irowvec>(nen);

  List w = x;
  List ret;


  // DIC stuff
  arma::mat Q2(Q.begin(),n,n,false);
  arma::cube TransProb = arma::cube(Q.nrow(),Q.nrow(),branchnumber);
  TransProb.zeros();
  NumericVector edgelengths = x["edge.length"];
  for(i=0;i<branchnumber;i++) TransProb.slice(i) = expmat(Q2*edgelengths(i));
  arma::mat matexpPL=PPmakePL(x,states,ne,Q,TransProb);
  double S=0;
  double X=0;
  // end DIC stuff

  // dwell times, transition counts, l01 &l10, root state, log(p(y|Q))
  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*n+2+1+1);

  for(i=0;i<N;i++){
    recordQ(&Q,&dwelltimes,i);

    treesamplebf(&Q,&rootdist,&B2,&B4,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);

    // DIC stuff

    // matrix exponentiation log likelihood
    for(j=0;j<branchnumber;j++) TransProb.slice(j) = expmat(Q2*edgelengths(j));
    // update the Partial Likelihood matrix 
    S=0;
    X=0;
    matexpPL=PPmakePLD(x,states,ne,Q,TransProb,&S);
    // calculate and record log likelihood
    for(j=0;j<n;j++) X=X+matexpPL(root-1,j)*rootdist(j);
    dwelltimes(i,n+n*n+2+1)=log(X)+S; 

    // \end DIC Stuff


    PPupdatel01(&Q,&Q2,&B2,&B4,&dwelltimes,i,Omega,prior);
    PPupdatel10(&Q,&Q2,&B2,&B4,&dwelltimes,i,Omega,prior);

    printf("%i \r",i);
  }


  return wrap(dwelltimes);

}


// PPmakePLD creates the partial likelihood matrix for independent tree samples
arma::mat PPmakePLksD(List x,arma::irowvec states,arma::irowvec ne,NumericMatrix Q,arma::cube TransProb,double *S){
  int Nnode = as<int>(x["Nnode"]);
  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int i;
  int j;
   for(i=0;i<states.size();i++) {
     if((states(i)-2*floor(states(i)/2))==0) {
       for(j=1;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
     if((states(i)-2*floor(states(i)/2))==1) {
       for(j=0;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
    }


  arma::imat edge =  as<arma::imat>(x["edge"]);
  arma::irowvec edge1 = trans(edge.col(0));
  arma::irowvec edge2 = trans(edge.col(1));

  int n = as<int>(x["Nnode"]);
  // for each node (row) calculate the probability of the tip states below that node conditioned on the node's state (column)
  for(i=0;i<n;i++) {
    PL.row(edge1(ne(2*i)-1)-1)=((TransProb.slice(ne(2*i)-1)*(PL.row(edge2(ne(2*i)-1)-1)).t())%(TransProb.slice(ne(2*i+1)-1)*(PL.row(edge2(ne(2*i+1)-1)-1)).t())).t();
    // normalized for big trees...
    (*S)=(*S)+log(sum(PL.row(edge1(ne(2*i)-1)-1)));
    PL.row(edge1(ne(2*i)-1)-1)=PL.row(edge1(ne(2*i)-1)-1)/sum(PL.row(edge1(ne(2*i)-1)-1));
  }
  return PL;
}

// [[Rcpp::export]]
NumericMatrix maketreelistMCMCksDICt(List& x,NumericMatrix& Q,NumericVector& pid,NumericMatrix& B,double Omega,IntegerVector& nen,IntegerVector& nodelist,int root,int N,NumericVector& prior) {

  RNGScope scope;

  int i;
  int f;
  double g;
  List maps=x["maps"];
  List mapnames=x["mapnames"];
  const int branchcount = maps.size();
  std::vector<Branch> brancharray(branchcount);
  for(i=0;i<maps.size();i++) brancharray[i] = makeabranch(maps(i),mapnames(i));
  int branchnumber=mapnames.size();  

 arma::imat edge =  as<arma::imat>(x["edge"]);
 arma::irowvec edge1 = trans(edge.col(0));
 arma::irowvec edge2 = trans(edge.col(1));
 int Nnode = as<int>(x["Nnode"]);
 
 int k=Q.nrow()/2-1;
 arma::rowvec lambdas = arma::zeros<arma::rowvec>(2);
 lambdas(0)=Q(0,1);
 lambdas(1)=Q(1,0); 
 arma::rowvec rkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) rkappas(i)=Q(2*i,2*i+2);
 arma::rowvec lkappas = arma::zeros<arma::rowvec>(k);
 for(i=0;i<k;i++) lkappas(i)=Q(2*i+2,2*i);
 arma::rowvec gammas = arma::zeros<arma::rowvec>(k+1);
 gammas(0)=1;
 for(i=1;i<=k;i++) gammas(i)=Q(2*i,2*i+1)/lambdas(0);

  arma::imat nodestatesmatrix=as<arma::imat>(x["node.states"]);
  arma::irowvec states = as<arma::irowvec>(x["states"]);

  arma::mat PL((2*Nnode+1),Q.nrow());
  PL.zeros();
  int j;
   for(i=0;i<states.size();i++) {
     if((states(i)-2*floor(states(i)/2))==0) {
       for(j=1;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
     if((states(i)-2*floor(states(i)/2))==1) {
       for(j=0;j<Q.nrow();j=j+2) PL(i,j) = 1;
     }
    }

  int n = B.nrow();
  arma::mat B2(B.begin(),n,n,false);
  arma::mat B4=trans(B2);
  
  arma::rowvec rootdist =as<arma::rowvec>(pid);
  arma::irowvec ne=as<arma::irowvec>(nen);


  // DIC stuff
  arma::mat Q2(Q.begin(),n,n,false);
  arma::cube TransProb = arma::cube(Q.nrow(),Q.nrow(),branchnumber);
  TransProb.zeros();
  NumericVector edgelengths = x["edge.length"];
  for(i=0;i<branchnumber;i++) TransProb.slice(i) = expmat(Q2*edgelengths(i));
  arma::mat matexpPL=PPmakePL(x,states,ne,Q,TransProb);
  double S=0;
  double X=0;
  // end DIC stuff

  
  // time spent in state, n**, l01, l10, kappas & gammas, root state, log(p(y|Q))
  arma::mat dwelltimes=arma::zeros<arma::mat>(N,n+n*n+2+3*k+1+1);
  int startcolumn=n+n*n+2+3*k+1+1+1;

  for(i=0;i<N;i++){
    recordQks(&Q,&dwelltimes,i,n,k);
    treesampleks(&Q,&rootdist,&B2,&B4,Omega,&ne,&brancharray,branchnumber,&nodestatesmatrix,&dwelltimes,root,&PL,nodelist,i,N,&states,&edge,&edge1,&edge2,Nnode);

    // DIC Stuff

    // loglikelihood conditional on the number of jumps on each branch, m
    //X=0;
    //for(j=0;j<n;j++) X=X+PL(root-1,j)*rootdist(j);
    //dwelltimes(i,n+n*n+2+3*k+1+1)=log(X)+S; 

    // matrix exponentiation log likelihood
    for(j=0;j<n;j++) for(f=0;f<n;f++) Q2(j,f)=Q(j,f);
    for(j=0;j<branchnumber;j++) TransProb.slice(j) = expmat(Q2*edgelengths(j));
    // update the Partial Likelihood matrix 
    S=0;
    X=0;
    matexpPL=PPmakePLksD(x,states,ne,Q,TransProb,&S);
    // calculate and record log likelihood
    for(j=0;j<n;j++) X=X+matexpPL(root-1,j)*rootdist(j);
    dwelltimes(i,n+n*n+2+3*k+1)=log(X)+S; 
    // end DIC stuff

    updateksl01(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    updateksl10(&Q,&B2,&B4,&dwelltimes,i,Omega,prior);
    for(j=0;j<k;j++)  updaterkappas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    for(j=1;j<=k;j++) updatelkappas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    for(j=1;j<=k;j++) updategammas(&Q,&B2,&B4,&dwelltimes,i,Omega,prior,j);
    printf("%i \r",i);
  }

  return wrap(dwelltimes);

}



















