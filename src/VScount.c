#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>

SEXP VScount(SEXP TH, SEXP cutoffs, SEXP m, SEXP B, SEXP c){

  int B_len = INTEGER(B)[0], m_len = INTEGER(m)[0], c_len = INTEGER(c)[0];
  int b, i, j;
  SEXP guessFT, THb, VS;

  PROTECT(guessFT=allocVector(INTSXP,1));
  PROTECT(THb=allocVector(REALSXP,m_len));
  PROTECT(VS=allocVector(INTSXP,c_len*B_len));

  for(b=0;b<B_len;b++){

    if((b%250==0.0) & (b>0.0))
      Rprintf("%d ",b);

    for(j=0;j<c_len;j++){
      INTEGER(guessFT)[0]=0;

      for(i=0;i<m_len;i++){
	REAL(THb)[i]=REAL(TH)[b*m_len+i];
	if(REAL(THb)[i]>REAL(cutoffs)[j]){
	  ++INTEGER(guessFT)[0];
	}
      }

      INTEGER(VS)[b*c_len+j]=INTEGER(guessFT)[0];
    }
  }

  Rprintf("%d\n",B_len);

  UNPROTECT(3);

  return(VS);
}
