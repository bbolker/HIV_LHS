#include <Rcpp.h>
using namespace Rcpp;

// note both . and _ have potential semantic meaning in C++
//  (I think _ is safer, though ...)

//[[Rcpp::export]]
CharacterVector getListElement(List stuff) {
    return as<CharacterVector>(stuff["a"]);
}

// getListElement(list(a="abc"))

//[[Rcpp::export]]
List g(double t, NumericVector yini, List parameters,
       NumericMatrix cc_mat) {
	
	int i,j,k; // counters

  // http://gallery.rcpp.org/articles/subsetting/    

  // should we (1) pass 1 long vector and split it?
  //           (2) pass 1 long vector and use vector indexing (maybe 
  // .. #define n(i) N[i]
  // .. #define nn(j,k) N[nrisk+k*nrisk+j]

  // FIXME: is there a more compact but still reasonably clear
  //    way to do this?  macro def?
  int nrisk = as<int>(parameters["n.risk"]);
  NumericVector rrisk(nrisk);
  rrisk = parameters["r.risk"];
  
  double rho = as<double>(parameters["rho"]);
	
  // N <- yini[1:n.risk]
  NumericVector N(nrisk);
  
  for (i=0; i<nrisk; i++){
  	
  	N[i] = yini[i];
  }
  

  // copy last (nrisk^2) elements into NN matrix
  // NN <- matrix(yini[-(1:n.risk)], n.risk, n.risk)
  // NN_adj <- NN
  // should be able to replace this with a single incantation
  // something like std::fill(NN.begin(), NN.end(), N[nrisk,N.end()]);

  NumericMatrix NN(nrisk,nrisk);
  NumericMatrix NN_adj(nrisk,nrisk);
	
 	for (i=0; i<nrisk; i++){
 		k = i * nrisk;
 		
 		for (j=0; j<nrisk; j++)  {
 			NN_adj = NN(i,j) = N[nrisk+k];
 			k++;
 		}
 	}
		
  // double diagonal
  // diag(NN_adj) <- 2 * diag(NN_adj)
  for (i=0; i<nrisk; i++){
  	NN_adj(i,i) *= 2;
  }
	
  //           frate <- r.risk * rho * N
  //           f.sum <- sum(frate)
	
  NumericVector frate(nrisk);
 	frate = rrisk*rho*N;
  double fsum = sum(frate);

  NumericVector dN(nrisk);

  //           dN <- - frate + colSums(cc_mat * NN_adj)
  for (i=0; i<nrisk; i++) {
		dN(i) = -frate(i);
    for (int j=0; j<nrisk; j++) {
	    dN[i] += cc_mat(i,j)*NN_adj(i,j);
		}
  }

  //           f.NN <- outer(frate, frate, "*")/f.sum 
  //           diag(f.NN) <- diag(f.NN)/2
  NumericMatrix fNN(nrisk,nrisk);
  for (i=0; i<nrisk; i++)  {
		for (j=0; j<nrisk; j++) {
	    fNN(i,j) = frate(i)*frate(j)/fsum;
			
	    if (i==j) {
				fNN(i,j) /= 2.0;
	  	}
		}
  }
    
  // don't really need this, should just modify fNN in place ... 
  //           dNN <- f.NN - cc_mat * NN
  // might be vectorizable ... ?
  NumericMatrix dNN(nrisk,nrisk);
  for (i=0; i<nrisk; i++)  {
		for (j=0; j<nrisk; j++) {
	    dNN(i,j) = fNN(i,j) - cc_mat(i,j)*NN(i,j);
		}
  }
	
  //           list(c(dN, dNN))
  
  return Rcpp::List::create(Rcpp::Named("dN") = dN,
                           Rcpp::Named("dNN") = dNN);
  
}