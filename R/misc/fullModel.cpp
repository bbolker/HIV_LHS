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
	
	int i,j,k,l;
	
	//partnership
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
	
	//infection
	NumericMatrix betaSI(nrisk, nrisk * nalpha);
	betaSI = as<NumericMatrix>(parameters["beta.SI"]);
	NumericMatrix p(nalpha, nalpha);
	p = as<NumericMatrix>(parameters["p"]);
	
	//epc infection
	NumericVector c_e(nalpha);
	c_e = parameters["c_e"];
	NumericVector c_u_ratio2(nrisk);
	c_u_ratio2 = parameters["c_u_ratio2"];
	NumericVector c_e_ratio2(nrisk);
	c_e_ratio2 = parameters["c_e_ratio2"];
	NumericVector c_u_ratio2I(nalpha * nrisk);
	c_u_ratio2I = parameters["c_u_ratio2.I"];
	NumericVector c_e_ratio2I(nalpha * nrisk);
	c_e_ratio2I = parameters["c_e_ratio2.I"];
	NumericMatrix c_u2(nrisk, nalpha);
	c_u2 = as<NumericMatrix>(parameters["c_u2"]);
	NumericMatrix rriskcouple(nrisk, nrisk);
	rriskcouple = as<NumericMatrix>(parameters["r.risk.couple"]);
	
	//mortality
	NumericVector lam2(nalpha * nrisk);
	lam2 = parameters["lam2"];
	NumericMatrix lammat_dis2(nalpha * nrisk, nalpha * nrisk);
	lammat_dis2 = as<NumericMatrix>(parameters["lammat_dis2"]);
	NumericMatrix lammat_adj2(nalpha * nrisk, nalpha * nrisk);
	lammat_adj2 = as<NumericMatrix>(parameters["lammat_adj2"]);
	
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
		j = floor(i/nalpha);
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
			fII(i, j) = rho2I[i] * I[i] * rho2I[j] * I[j]/fsum;
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
		j = floor(i/nalpha);
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
	
	//Infection
	
	//##within couple infection
	//inf.SI.from <- beta.SI * SI
	//inf.SI.tmp <- distribute.SIinf(inf.SI.from)
	//inf.SI.to <- inf.SI.tmp + t(inf.SI.tmp)
	//diag(inf.SI.to) <- diag(inf.SI.to)/2
	//distribute.SIinf <- function(SIinf){
	//	attach(pp)
	//	dist.SI <- matrix(0, n.risk * n.alpha, n.risk * n.alpha)
	//	for(i in 1:n.risk){
	//		vec.i <- c(1:n.alpha) + n.alpha * (i-1)
	//		for(j in 1:n.risk){
	//			vec.j <- c(1:n.alpha) + n.alpha * (j-1)
	//			dist.SI[vec.i,vec.j] <- SIinf[i,vec.j] * p
	//		}
	//	}
	//	detach(pp)
	//		return(dist.SI)
	//}
	NumericMatrix infSIfrom(nrisk, nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			infSIfrom(i, j) = betaSI(i,j) * SI(i, j);
		}
	}
	
	NumericMatrix infSItotmp(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		k = i % nalpha;
		for(j = 0; j < nrisk * nalpha; j ++){
			l = j % nalpha;
			infSItotmp(i,j) = infSIfrom(floor(i/nalpha), (k + floor(j/nalpha) * nalpha)) * p(k, l);
		}
	}
	
	NumericMatrix infSIto(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		k = i % nalpha;
		for(j = 0; j < nrisk * nalpha; j ++){
			l = j % nalpha;
			infSIto(j,i) = infSItotmp(i,j) + infSItotmp(j,i);
		}
		infSIto(i,i) /= 2;
	}
	
	//prop.mat <- matrix(c(c_u_ratio2.I * I + c_e_ratio2.I * (colSums(SI) + colSums(II_adj))), ncol = n.risk)
	//prop.num <- rowSums(prop.mat)
	//prop.denom <- sum(c_u_ratio2 * N + c_e_ratio2 * NN)
	//prop <- prop.num/prop.denom
	
	NumericMatrix propmat(nalpha, nrisk);
	for(i = 0; i < nrisk * nalpha; i ++){
		j = i % nalpha;
		k = floor(i/nalpha);
		propmat(j,k) = c_u_ratio2I[i] * I[i] + c_e_ratio2I[i] * (sum(SI(_,i)) + sum(II_adj(_,i)));
	}
	
	double propdenom = sum(c_u_ratio2 * N + c_e_ratio2 * NN);
	NumericVector prop(nalpha);
	for(i = 0; i < nalpha; i++){
		prop[i] = sum(propmat(i,_))/propdenom;
	}
	
	//inf.S.rate <- S * sweep(c_u2, 2, prop, "*")
	//inf.S.from <- rowSums(inf.S.rate)
	//inf.S.to <- c(apply(inf.S.rate, 1, function(x) x %*% p))
	NumericMatrix infSrate(nrisk, nalpha);
	NumericVector infSfrom(nrisk);
	for(i = 0; i < nrisk; i ++){
		infSrate(i,_) = S[i] * c_u2(i,_) * prop;
		infSfrom[i] = sum(infSrate(i,_));
	}	
	NumericVector infSto(nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		int k = i * nalpha;
		for(j = 0; j < nalpha; j ++){
			infSto[k] = sum(infSrate(i,_) * p(_,j));
			k++;
		}
	}
	
	//inf.ce.tmp <- c_e * prop
	//risk.SS <- r.risk.couple * SS
	//inf.SS.from <- risk.SS * sum(inf.ce.tmp)
	NumericVector infcetmp(nalpha);
	infcetmp = c_e * prop;
	NumericMatrix infSSfrom(nrisk, nrisk);
	for(i = 0; i < nrisk; i ++){
		infSSfrom(_,i) = rriskcouple(_,i) * SS(_,i) * sum(infcetmp);
	}
	
	//ce.mut <- inf.ce.tmp %*% p
	//inf.SS.tmp2 <- outer(c((r.risk * SS_adj)), c(ce.mut))
	//inf.SS.to <- matrix(c(t(inf.SS.tmp2)), nrow = n.risk, byrow = TRUE)
	NumericVector cemut(nalpha);
	for(i = 0; i < nalpha; i ++){
		cemut[i] = sum(infcetmp * p(_,i));
	}
	//FIXME:
	
	NumericMatrix infSStmp2(nrisk * nrisk, nalpha);
	for(i = 0; i < nrisk; i ++){
		k = nrisk * i;
		for(j = 0; j < nrisk; j ++){
			infSStmp2(k,_) = rrisk[j] * SS_adj(j,i) * cemut;
			k++;
		}
	}
	
	NumericMatrix infSSto(nrisk, nrisk * nalpha);
	for(i= 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			k = j % nalpha;
			l = i * nrisk + floor(j/nalpha);
			infSSto(i,j) =  infSStmp2(l,k);
		}
	}
	
	//SI2 <- r.risk * SI
	//inf.SI.from2 <- sum(inf.ce.tmp) * SI2
	//inf.SI.to2 <- matrix(0, n.risk * n.alpha, n.risk * n.alpha)
	//for(i in 1:n.risk){
  //	tmp.vec <- c(1:n.alpha) + (i-1) * n.alpha
	//	inf.SI.to2[,tmp.vec] <- outer(SI2[i,], ce.mut)
	//}
	//inf.SI.to2 <- inf.SI.to2 + t(inf.SI.to2)
	//diag(inf.SI.to2) <- diag(inf.SI.to2)/2
	NumericMatrix infSIfrom2(nrisk, nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		infSIfrom2(i,_) = rrisk[i] * sum(infcetmp) * SI(i,_);
	}
	NumericMatrix infSIto2tmp(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		j = i % nalpha;
		k = floor(i/nalpha);
		infSIto2tmp(_,i) = rrisk[k] * SI(k,_) * cemut[j];
	}
	NumericMatrix infSIto2(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			infSIto2(i,j) = infSIto2tmp(j,i) + infSIto2tmp(i,j);
		}
		infSIto2(i,i) /= 2;
	}
	
	//mortality
	//mort.I.out <- lam2 * I
	//mort.I.sus <- colSums2(mort.I.out, n.alpha)
	NumericVector mortIout(nrisk * nalpha);
	mortIout = lam2 * I;
	NumericVector mortIsus(nrisk);
	for(i = 0; i < nrisk * nalpha; i ++){
		j = floor(i/nalpha);
		mortIsus[j] += mortIout[i];
	}
	
	//mort.SI.out <- sweep(SI, 2, lam2, "*")
	//mort.SI.sus <- rowSums(mort.SI.out) + colSums2(mort.SI.out, n.alpha)
	NumericMatrix mortSIout(nrisk, nrisk * nalpha);
	NumericVector mortSIsus(nrisk);
	for(i = 0; i < nrisk; i ++){
  	mortSIout(i,_) = SI(i,_) * lam2;
		mortSIsus[i] = sum(mortSIout(i,_));
	}
	
	for(i = 0; i < nrisk * nalpha; i ++){
		j = floor(i/nalpha);
		mortSIsus[j] += sum(mortSIout(_,i));
	}
	
	//mort.II.out <- lammat_adj2 * II
	//mort.II.I <- colSums(lammat_dis2 * II)
	//mort.II.sus <- colSums2(t(lammat_dis2 * II), n.alpha)
	NumericMatrix mortIIout(nrisk * nalpha, nrisk * nalpha);
	NumericVector mortIII(nrisk * nalpha);
	NumericVector mortIIsus(nrisk);
	for(i = 0; i < nrisk * nalpha; i ++){
		mortIIout(i,_) = lammat_adj2(i,_) * II(i,_);
	}
	for(i = 0; i < nrisk * nalpha; i ++){
		j = floor(i/nalpha);
	 	mortIII[i] = sum(lammat_dis2(_,i) * II(_,i));
	 	mortIIsus[j] += sum(lammat_dis2(i,_) * II(i,_));
	}
	 
	//dS <- - rho2 * S + colSums(cc_mat * SS_adj) + rowSums(d.SI) - inf.S.from + mort.I.sus +
	//	mort.SI.sus + mort.II.sus
	NumericVector dS(nrisk);
	dS = - rho2 * S;
	for(i = 0; i < nrisk; i ++){
		dS[i] += sum(cc_mat(_,i) * SS_adj(_,i)) + sum(d_SI(i,_)) - infSfrom[i] + mortIsus[i] + mortSIsus[i] + mortIIsus[i];
	}
	
	//dI <- - rho2.I * I + colSums(d.SI) + colSums(cc_mat.II * II_adj) + inf.S.to - mort.I.out + mort.II.I
	NumericVector dI(nrisk * nalpha);
	dI = - rho2I * I;
	for(i = 0; i < nrisk * nalpha; i ++){
		dI[i] += sum(d_SI(_,i)) + sum(cc_matII(_,i) * II_adj(_,i)) + infSto[i] - mortIout[i] + mortIII[i];
	}
	
	//dSS <- f.SS - d.SS - inf.SS.from
	NumericMatrix dSS(nrisk, nrisk);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk; j ++){
			dSS(i,j) = fSS(i, j) - d_SS(i, j) - infSSfrom(i,j);
		}
	}
	
	
	//dSI <- f.SI - d.SI - inf.SI.from + inf.SS.to - inf.SI.from2 - mort.SI.out
	NumericMatrix dSI(nrisk, nrisk * nalpha);
	for(i = 0; i < nrisk; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			dSI(i, j) = fSI(i, j) - d_SI(i,j) - infSIfrom(i,j) + infSSto(i,j) - infSIfrom2(i,j) - mortSIout(i,j);
		}
	}
	
	//dII <- f.II - cc_mat.II * II + inf.SI.to + inf.SI.to2 - mort.II.out
	NumericMatrix dII(nrisk * nalpha, nrisk * nalpha);
	for(i = 0; i < nrisk * nalpha; i ++){
		for(j = 0; j < nrisk * nalpha; j ++){
			dII(i,j) = fII(i,j) - cc_matII(i,j) * II(i,j) + infSIto(i,j) + infSIto2(i,j) - mortIIout(i,j);
		}
	}
	
	return Rcpp::List::create(Rcpp::Named("dS") = dS,
                           Rcpp::Named("dI") = dI,
                           Rcpp::Named("dSS") = dSS,
                           Rcpp::Named("dSI") = dSI,
                           Rcpp::Named("dII") = dII,
                           Rcpp::Named("test") = mortSIout);
	
}


