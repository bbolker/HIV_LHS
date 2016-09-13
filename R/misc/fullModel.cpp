#include <Rcpp.h>
using namespace Rcpp;

// note both . and _ have potential semantic meaning in C++
//  (I think _ is safer, though ...)

//[[Rcpp::export]]
CharacterVector getListElement(List stuff) {
    return as<CharacterVector>(stuff["a"]);
}

// getListElement(list(a="abc"))

NumericMatrix outer2 (NumericVector x,
                  NumericVector y){
	unsigned int xn = x.length();
	unsigned int yn = y.length();
	
	NumericMatrix res(xn,yn);

	for(unsigned int i = 0; i < xn; i ++){
		for(unsigned int j = 0; j < yn; j ++){
			res(i, j) = x[i] * y[j];
		}
	}
	return res;
}

//[[Rcpp::export]]
List g(double t, NumericVector yini, List parameters) {
	
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
  NumericMatrix cc_mat(nrisk, nrisk);
  cc_mat = as<NumericMatrix>(parameters["cc_mat"]);
	
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
 			NN_adj(j,i) = NN(j,i) = yini[nrisk+k];
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

//[[Rcpp::export]]
List g2(double t, NumericVector yini, List parameters) {
	
	int i,j,k;
	
	int nrisk = as<int>(parameters["n.risk"]);
	int nalpha = as<int>(parameters["n.alpha"]);
	NumericVector rrisk(nrisk);
	rrisk = parameters["r.risk"];
	NumericMatrix cc_mat(nrisk, nrisk);
	cc_mat = as<NumericMatrix>(parameters["cc_mat"]);
	NumericVector rho2 = parameters["rho2"];
	NumericVector rho2I = parameters["rho2.I"];
	NumericMatrix cc_matI(nrisk, nrisk * nalpha);
	cc_matI = as<NumericMatrix>(parameters["cc_mat.I"]);
	NumericMatrix cc_matII(nrisk * nalpha, nrisk * nalpha);
	cc_matII = as<NumericMatrix>(parameters["cc_mat.II"]);
	
	//S <- yini[S.ind]
	NumericVector S(nrisk);
	for (i=0; i<nrisk; i++){
		S[i] = yini[i];
	}
	
	//I <- yini[I.ind]
	NumericVector I(nrisk*nalpha);
	for(i = 0; i<nrisk*nalpha; i++){
		j = i + nrisk;
		I[i] = yini[j];
	}
	
	//SS <- matrix(yini[SS.ind], n.risk, n.risk)
	//SS_adj <- SS
	//diag(SS_adj) <- 2 * diag(SS_adj)
	NumericMatrix SS(nrisk, nrisk);
	NumericMatrix SS_adj(nrisk, nrisk);
	for(i = 0; i< nrisk; i++){
		k = i * nrisk + nrisk+nrisk*nalpha;
		for(j=0;j <nrisk;j++){
			SS_adj(j,i) = SS(j, i) = yini[k];
			k++;
		}
		SS_adj(i,i) *= 2;
	}
	
	//SI <- matrix(yini[SI.ind], n.risk, n.risk * n.alpha)
	NumericMatrix SI(nrisk, nrisk*nalpha);
	for(i = 0; i < nrisk*nalpha; i++){
		k = i * nrisk + nrisk + nrisk * nalpha + nrisk * nrisk;
		for(j = 0; j < nrisk; j++){
			SI(j, i) = yini[k];
			k++;
		}
	}
	
	//II <- matrix(yini[II.ind], n.risk * n.alpha, n.risk * n.alpha)
	//II_adj <- II
	//diag(II_adj) <- 2 * diag(II_adj)
	NumericMatrix II(nrisk*nalpha, nrisk*nalpha);
	NumericMatrix II_adj(nrisk*nalpha, nrisk*nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		k = i * nrisk * nalpha + nrisk + nrisk * nalpha + nrisk * nrisk + nrisk * nrisk * nalpha;
		for(j = 0; j < nrisk * nalpha; j++){
			II_adj(j,i) =II(j, i) = yini[k];
			k++;
		}
		II_adj(i,i) *= 2;
	}
	
	//Partnership formation rate
	
	//N <- S + colSums2(I, n.alpha)
	NumericVector N(nrisk);
	for(i = 0; i < nrisk; i ++){
		N[i] = S[i];
	}
	for(i = 0; i < nrisk * nalpha; i ++){
		j = ceil(i/nalpha);
		N[j] += I[i];
	}
	
	//frate <- rho2 * N
	NumericVector frate(nrisk);
	frate = rho2*N;
	
	//fsum <- sum(frate)
	double fsum = std::accumulate(frate.begin(),
	                              frate.end(), 0.0);
	
	//f.SS <- outer(rho2 * S, rho2 * S)/f.sum
	//diag(f.SS) <- diag(f.SS)/2
	NumericMatrix fSS(nrisk, nrisk);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j <nrisk; j ++){
			fSS(i,j) = rho2[i] * S[i] * rho2[j] * S[j]/fsum;
		}
		fSS(i,i) /= 2;
	}
	
	//f.SI <- outer(rho2 * S, rho2.I * I)/f.sum
	NumericMatrix fSI(nrisk, nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			fSI(i, j) = rho2[i] * S[i] * rho2I[j] * I[j]/fsum;
		}
	}
		
	//f.II <- outer(rho2.I * I, rho2.I * I)/f.sum
	//diag(f.II) <- diag(f.II)/2
	NumericMatrix fII(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			fII(i, j) = rho2I[i] * I[i] * rho2[j] * I[j]/fsum;
		}
		fII(i,i) /= 2;
	}
	
	//NN <- colSums(SS_adj) + rowSums(SI) + colSums2(SI, n.alpha) + colSums2(II_adj, n.alpha)
	NumericVector NN(nrisk);
	for(i = 0; i < nrisk; i ++){
		NN[i] = sum(SS_adj(_,i));
		NN[i] += sum(SI(i,_));
	}
	
	for(i = 0; i < nrisk * nalpha; i ++){
		j = ceil(i/nalpha);
		NN[j] += sum(SI(_,i)) + sum(II_adj(_,i));
	}
	
	//d.SS <- cc_mat * SS
	NumericMatrix d_SS(nrisk, nrisk);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk; j ++){
			d_SS(i,j) = cc_mat(i,j) * SS(i,j);
		}
	}
	
	//d.SI <- cc_mat.I * SI
	NumericMatrix d_SI(nrisk, nrisk*nalpha);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			d_SI(i,j) = cc_matI(i,j) * SI(i,j);
		}
	}
	
	//dS <- - rho2 * S + colSums(cc_mat * SS_adj) + rowSums(d.SI) - inf.S.from + mort.I.sus +
	//	mort.SI.sus + mort.II.sus
	NumericVector dS(nrisk);
	dS = - rho2 * S;
	for(i = 0; i < nrisk; i ++){
		dS[i] += sum(cc_mat(_,i) * SS_adj(_,i));
		dS[i] += sum(d_SI(i,_));
	}
	
	//dI <- - rho2.I * I + colSums(d.SI) + colSums(cc_mat.II * II_adj) + inf.S.to - mort.I.out + mort.II.I
	NumericVector dI(nrisk * nalpha);
	dI = - rho2I * I;
	for(i = 0; i < nrisk * nalpha; i ++){
		dI[i] += sum(d_SI(_,i));
		dI[i] += sum(cc_matII(_,i) * II_adj(_,i));
	}
	
	//dSS <- f.SS - d.SS - inf.SS.from
	NumericMatrix dSS(nrisk, nrisk);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk; j ++){
			dSS(i,j) = fSS(i, j) - d_SS(i, j);
		}
	}
	
	
	//dSI <- f.SI - d.SI - inf.SI.from + inf.SS.to - inf.SI.from2 - mort.SI.out
	NumericMatrix dSI(nrisk, nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			dSI(i, j) = dSI(i, j) + d_SI(i,j);
		}
	}
	
	//dII <- f.II - cc_mat.II * II + inf.SI.to + inf.SI.to2 - mort.II.out
	NumericMatrix dII(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			dII(i,j) = fII(i,j) - cc_matII(i,j) * II(i,j);
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("dS") = dS,
                           Rcpp::Named("dI") = dI,
                           Rcpp::Named("dSS") = dSS,
                           Rcpp::Named("dSI") = dSI,
                           Rcpp::Named("dII") = dII);
	
}


