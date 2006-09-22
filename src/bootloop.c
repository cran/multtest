#include <R.h>
#include <Rinternals.h>
#include <math.h>

void bootloop(double *X, double *W, int *p, int *n, int *B, double *muboot, int *samp, SEXP fbody, SEXP fenv){

  int B_len=*B, p_len=*p, num_samples=*n, b, i, j;
  SEXP Xb, Wb, Sb, Tb;

  PROTECT(Xb=allocVector(REALSXP,num_samples));
  PROTECT(Wb=allocVector(REALSXP,num_samples));
  PROTECT(Sb=allocVector(INTSXP,num_samples));
  PROTECT(Tb=allocVector(REALSXP,3));


  for(b=0;b<B_len;b++){
    if((b%100==0.0) & (b>0.0)) /* modulo 100 */
      Rprintf("%d ",b);
    for(j=0;j<p_len;j++){
      for(i=0;i<num_samples;i++){
	INTEGER(Sb)[i]=samp[num_samples*b+i];
	REAL(Xb)[i]=X[(samp[num_samples*b+i]-1)*p_len+j];
	REAL(Wb)[i]=W[(samp[num_samples*b+i]-1)*p_len+j];
      }

      defineVar(install("samp"),Sb,fenv);
      defineVar(install("x"),Xb,fenv);
      defineVar(install("w"),Wb,fenv);

      Tb=eval(fbody,fenv);
      muboot[p_len*b+j]=REAL(Tb)[0]/REAL(Tb)[1];
    }
  }

  Rprintf("%d\n",B_len);

  UNPROTECT(4);
}

