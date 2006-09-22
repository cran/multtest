/*****************************************************************/
/*           Header files                                        */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include "mt.h"
/************************************************************************************/
/*                      malloc_gene_data                                            */
/************************************************************************************/
/*Allocate the necessary space for the big data,
  see the comments about the structrue GENE_DATA
*/
int myDEBUG=0;
long int g_random_seed=3455660;
void print_b(int b,int B,char* prompt){
  static int p=0;
  if(b==0) p=0;
  if(!PROMPT_LEN){
    if((B<=100) ||(b%(B/100)==0))
      {
	   /*fprintf(stderr,"%s%d\t",prompt,b);*/
	Rprintf("%s%d\t",prompt,b);
	p++;
	if(PRINT_VAR_NUM && (p%PRINT_VAR_NUM==0))
	     /*fprintf(stderr,"\n");*/
	     Rprintf("\n");

      }
  }else if((b%(PROMPT_LEN+1))==0){ /*use PROMPT_LEN+1 to avoid the compiling warnings*/
    p++;
    /*fprintf(stderr,"%s%d",prompt,b);*/
    Rprintf("%s%d",prompt,b);
    if(PRINT_VAR_NUM && (p%PRINT_VAR_NUM==0))
	 /*fprintf(stderr,"\n");*/
	 Rprintf("\n");

  }
}
#ifdef WINDOWS
void win_print(FILE* fp, char* format,...)
{
  va_list ap;
  va_start(ap, format);
  REvprintf(format, ap);
  va_end(ap);
}
#endif
void malloc_gene_data(GENE_DATA* pdata)
{
  int i;
  int nrow=pdata->nrow;
  int ncol=pdata->ncol;

  assert(pdata->id=(char**)Calloc(nrow,char*));
  assert(pdata->d=(float**)Calloc(nrow,float*));
  assert(pdata->L=(int*)Calloc(ncol,int));

  /*initialization*/
  memset(pdata->L,0,sizeof(int)*ncol);
  for(i=0;i<ncol;i++) 
    pdata->L[i]=0;

  for (i=0; i<nrow; i++) {
    assert(pdata->id[i] = (char *) Calloc(MAX_ID,char));
    assert(pdata->d[i]=(float *) Calloc(ncol,float));
 }
}

/********************************************************************************/
/*                           free_gene_data                                      */
/********************************************************************************/
/*free the space allocated for pdata*/
void free_gene_data(GENE_DATA* pdata)
{
  int i;
  for (i=0; i<pdata->nrow; i++) {
    Free(pdata->d[i]);
    Free(pdata->id[i]);
  }
  Free(pdata->L);
  Free(pdata->d);
  Free(pdata->id);
}

/********************************************************************************/
/*                   compute_test_stat                                          */
/********************************************************************************/
void compute_test_stat(GENE_DATA* pdata, int* L,float* T,
		       FUNC_STAT func_stat,const void* extra)
     /*L is an array which contains 0,1,2 ... for specifying class label*/
     /*T is the test_test needs to return*/
     /*func_stat is a functin pointer with the following protocol
       float func_stat(float *Y, int* L,int n, float na,const void* extra)*/
{
  int i;
  for(i=0;i<pdata->nrow;i++){
    /*    fprintf(stderr,"i=%d,T[%d]=%d\n",i,i,T[i]);*/
    T[i]=(*func_stat)(pdata->d[i],L,pdata->ncol,pdata->na,extra);
  }
}
 
/********************************************************************************/
/*         get1pvalue                                                           */
/********************************************************************************/
/*  we'll do complete resampling with the function next_sample, which is determined by
    the function next_sample when it returns false.
     in the next_sample, you need to decide to choose complet resampling or not or 
     whatever you like

  L: is the labelling of each experiment
  T: is that test statistics
  P: unadjtesed P-values

  To use the function first_sample, and next_sample, they're needed to write into a separate 
  file, where it provides the create_sampling to do initialization such as allocate the space
  (before use the sampling) and delete_sampling after we've done the sampling in the main().

  int first_sample(int *L)
       get the first sample of the labelling.
  int next_sample(int* L)
       get the next sample, if it's done all the sampling, then it returns 0, 
       otherwise it returns 1.


  input: pdata, L, B,next_sample,func_stat
  output: T,P
 */

void  get1pvalue(GENE_DATA* pdata,int* L,float* T,float* P,
		 FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
		 FUNC_SAMPLE func_next_sample,FUNC_CMP func_cmp,const void* extra)
{
  int b=0,*bL,i,is_next,*total;
  float *bT, *count;
  int ncol=pdata->ncol;
  int nrow=pdata->nrow;
  int B=(*func_first_sample)(NULL);
  /*allocate the space and initialziation*/
  assert(bT=(float*)Calloc(nrow,float));
  assert(bL=(int*)Calloc(ncol,int));
  assert(count=(float*)Calloc(nrow,float));
  memset(count,0,sizeof(float)*nrow); 
  assert(total=(int*)Calloc(nrow,int));
  memset(total,0,sizeof(int)*nrow);

  /*comuter the original one first*/
  compute_test_stat(pdata,L,T,func_stat,extra);

  /*iteration for permutaion*/
  (*func_first_sample)(bL);
  is_next=1;
  b=0;
  while(is_next){
    compute_test_stat(pdata,bL,bT,func_stat,extra);
    for(i=0;i<nrow;i++){
      if(bT[i]==NA_FLOAT) continue;
      if(T[i]==NA_FLOAT) continue;
      /* right now I only implements the 3 cases, which are pretty common*/
      if((func_cmp==cmp_high) &&(bT[i]>=T[i]-EPSILON)){
	count[i]+=1;
      }else if((func_cmp==cmp_low) &&(bT[i]<=T[i]+EPSILON)){
	count[i]+=1;
      }else if ((func_cmp==cmp_abs) &&(fabs(bT[i])>=fabs(T[i])-EPSILON)){
	count[i]+=1;
      }	      
      total[i]++;
     }
    b++;
    print_b(b,B,"b=");
    is_next=(*func_next_sample)(bL);
  }

  /*summarize the results*/
  for(i=0;i<nrow;i++){ 
    if(total[i]==0) 
      P[i]=NA_FLOAT;
    else P[i]=count[i]*1.0/total[i];
  }

  /*free the spaces*/
  Free(bT);
  Free(count);
  Free(total);
  Free(bL);
}

/********************************************************************************/
/*                    sort_gene_data                                            */
/********************************************************************************/
/* Description:
   sort the rows of gene_data such that row R[i] of is the first row, i=0,...,m-1,
   wher R[0],...,R[m-1] is a permutation of (0,...,m-1)
*/
										 
void sort_gene_data(GENE_DATA* pdata,int*R)
{
  int i,nrow=pdata->nrow;
  char** old_id; /*the old addresses of the gene id*/
  float** old_d;  /*th old addresses of the gene data*/
  assert(old_d=(float**)Calloc(nrow,float*));
  assert(old_id=(char**)Calloc(nrow,char*));
  /*store the original pointers from pdata*/ 
  for(i=0;i<nrow;i++)
    {
      old_d[i]=pdata->d[i];
      old_id[i]=pdata->id[i];
    }
  /*rearrange the data so that it's ordered according to R*/
  for(i=0;i<nrow;i++)
    {
      pdata->d[i]=old_d[R[i]];
      pdata->id[i]=old_id[R[i]];
    }
  Free(old_id);
  Free(old_d);
}
/********************************************************************************/
/*                    sort_vector                                               */
/********************************************************************************/
/* Desciption
      sort the vector V according to the order R with n elemnets
      where R[0],...,R[n-1] is a permutation of 0,...n-1
*/
void sort_vector(float* V,int*R,int n)
{
  float* old_V;
  int i;
  assert(old_V=(float*)Calloc(n,float));
  for(i=0;i<n;i++)
    old_V[i]=V[i];
  for(i=0;i<n;i++)
    V[i]=old_V[R[i]];
  Free(old_V);
}
/********************************************************************************/
/*                    get_all_samples_P                                         */
/********************************************************************************/
/*  Descriptions: Try to get all the unadjusted p-values for a gene with
    experssion values at V with n experiemtns. 

    int first_sample(int *L)
       get the first sample of the labelling.
       if L==NULL, then it returns all the possible simulations, which depends on
       the initial function create_sampling.
    int next_sample(int* L)
       get the next sample, if it's done all the sampling, then it returns 0, 
       otherwise it returns 1.
   output is P*/
void get_all_samples_P(float* V, int n,float* P,float na, 
		       FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
		       FUNC_SAMPLE func_next_sample,FUNC_CMP func_cmp,const void* extra)
{
  int  *L,*R,i,oldb,is_next,b=0,B_new,B;
  float* T=P,oldf;/*it first stores T, then switch to P*/
 
  B=(*func_first_sample)(NULL);
  /*allocate the spaces*/
  assert(L=(int*)Calloc(n,int));
  assert(R=(int*)Calloc(B,int));

  /*compute all the test_stat*/
  (*func_first_sample)(L);
  is_next=1;
  B_new=0;
  while(is_next){
    T[b]=func_stat(V,L,n,na,extra);
    if(T[b]!=NA_FLOAT)
      B_new++;
    b++;
    is_next=(*func_next_sample)(L);
  }
  if(B!=b){
    fprintf(stderr,"Error we have b(%d)!=B(%d)\n",b,B);
    return;/*exit(1)*/;
  }
  if(myDEBUG)
      print_farray(stderr,T,B);
  /*order the test_stat*/
  order_data(T,R,B,func_cmp);

  /*note the last elements of B-B_new has NA T-value*/
  /*assign the probabilites*/
  oldb=0;
  oldf=T[R[0]];
  for(b=1;b<B_new;b++){
    if((func_cmp==cmp_high)&(T[R[b]]>=oldf-EPSILON)) continue;
    else if ((func_cmp==cmp_low ) &&(T[R[b]]<=oldf+EPSILON)) continue;
    else if((func_cmp==cmp_abs )&& fabs(T[R[b]])>=fabs(oldf)-EPSILON) continue;

    for(i=oldb;i<b;i++)
      P[R[i]]=(b+0.0)/B_new;
    oldb=b;
    if(b<B_new-1) oldf=T[R[b]];
  }
  for(i=oldb;i<b;i++)
    P[R[i]]=1.0;
  
  /*for NA test_stat, assign NA probabilites*/
  for(b=B_new;b<B;b++)
    P[R[b]]=NA_FLOAT;

  /*free the space*/
  Free(L);
  Free(R);
} 

/*get all the samples of T and they're also ordered.It's used only for diagonsis*/
void get_all_samples_T(float* V, int n,float* T,float na, 
		       FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
		       FUNC_SAMPLE func_next_sample,const void* extra)
{
  int  *L,*R,is_next,b=0,B;
 
  B=(*func_first_sample)(NULL);
  /*allocate the spaces*/
  assert(L=(int*)Calloc(n,int));
  assert(R=(int*)Calloc(B,int));

  /*compute all the test_stat*/
  (*func_first_sample)(L);
  is_next=1;
  while(is_next){
    T[b]=func_stat(V,L,n,na,extra);
    b++;
    is_next=(*func_next_sample)(L);
  }
  if(B!=b){
    fprintf(stderr,"Error we have b(%d)!=B(%d)\n",b,B);
    return;/*exit(1)*/;
  }
  if(myDEBUG)
    print_farray(stderr,T,B);
  Free(L);
  Free(R);
} 

void adj_pvalue_quick(GENE_DATA* pdata,float*T, float* P, 
		      float* Adj_P,float* Adj_Lower,
		      FUNC_STAT func_stat,FUNC_STAT func_stat_T,
		      FUNC_SAMPLE func_first_sample, 
		      FUNC_SAMPLE func_next_sample,FUNC_CMP func_cmp,const void* extra)
{

  int *L,b,B,B_new,i,*R,neq; /*b for simulation*, neq is for the number of equal signs*/
  float* all_P,*all_Q,count;
  int ncol=pdata->ncol,nrow=pdata->nrow;
   
  /*allocate the space*/
  B=(*func_first_sample)(NULL);
  assert(L=(int*)Calloc(ncol,int)); 
  assert(R=(int*)Calloc(nrow,int));
  assert(all_P=(float*)Calloc(B,float));
  assert(all_Q=(float*)Calloc(B,float));

  /*get the original unadjusted p-values first
   we'll use the normalized t-statistics*/
  get1pvalue(pdata,pdata->L,T,P,func_stat_T,func_first_sample,func_next_sample,func_cmp,extra);
  if(myDEBUG)
    {
      print_farray(stderr,T,pdata->nrow);
      print_farray(stderr,P,pdata->nrow);
    }
  /*sort the test_stat*/
  order_mult_data(R,nrow,2,P,cmp_low,T,func_cmp);
  /*order_data(P,R,nrow,func_cmp);*/

  /*rearrange the data according the unadjusted p-values*/
  sort_gene_data(pdata,R);
  sort_vector(T,R,nrow);
  sort_vector(P,R,nrow);
  
  /*initialze all_Q[]=NA_FLOAT*/
  for(b=0;b<B;b++)
    all_Q[b]=NA_FLOAT;

  /*loop for each gene*/
  for(i=nrow-1;i>=0;i--){
    get_all_samples_P(pdata->d[i],ncol,all_P,pdata->na,
		      func_stat,func_first_sample,func_next_sample,func_cmp,extra);
    if(myDEBUG)
      print_farray(stderr,all_P,B);

    /*update all_Q*/
    count=0;
    B_new=0;
    neq=0;

    for(b=0;b<B;b++){
      if (all_P[b]==NA_FLOAT) break;/*we don't need care about NA pvlaues*/
      if(all_Q[b]>all_P[b])
	all_Q[b]=all_P[b];/*update q* by the value p*/
      if(all_Q[b]==NA_FLOAT) continue;/*skip NA q*/
      if(all_Q[b]<P[i]){
	count+=1;
      }else if (all_Q[b]<=P[i]+EPSILON)/*it'd already > */
	 neq++;
      B_new++;
    }

    if(myDEBUG)
      {
	print_farray(stderr,all_Q,B);
	fprintf(stderr,"P[%d]=%5.3f,count=%5.2f,neq=%d\n",i,P[i],count,neq);
      }

    /*assign the Adj_P and Adj_Lower for gene i */
    if(B_new!=0) {
      Adj_P[i]=(count+neq)/B_new;

      if(neq==0) 
	Adj_Lower[i]=count/B_new;
      else Adj_Lower[i]=(count+1)/B_new;
    }
    else {
      Adj_P[i]=NA_FLOAT;
      Adj_Lower[i]=NA_FLOAT; 
    }
    /*************************** */
    print_b((nrow-i),nrow,"r="); 
  }

  /* to make monotone of Adj_P and Adj_Lower*/
  for(i=1;i<nrow;i++)
    if(Adj_P[i]<Adj_P[i-1])
      Adj_P[i]=Adj_P[i-1];

  for(i=1;i<nrow;i++)
    if(Adj_Lower[i]<Adj_Lower[i-1])
      Adj_Lower[i]=Adj_Lower[i-1];

  /*free the spaces*/
  Free(L);
  Free(R);
  Free(all_P);
  Free(all_Q);
}
	
/********************************************************************************/
/*         adj-by_t                                                          */
/********************************************************************************/
/*  we'll do complete resampling with the function next_sample, which is determined by
    the function next_sample when it returns false.
     in the next_sample, you need to decide to choose complet resampling or not or 
     whatever you like

  L: is the labelling of each experiment
  T: is that test statistics
  P: unadjtesed P-values
  Adj_P:ajusted p-values by using the max|T|
  To use the function first_sample, and next_sample, they're needed to write into a separate 
  file, where it provides the create_sampling to do initialization such as allocate the space
  (before use the sampling) and delete_sampling after we've done the sampling in the main().

  int first_sample(int *L)
       get the first sample of the labelling.
  int next_sample(int* L)
       get the next sample, if it's done all the sampling, then it returns 0, 
       otherwise it returns 1.


  input: pdata, L, B,next_sample,func_stat
  output: T,P
 */

void  adj_by_T(GENE_DATA* pdata,float* T,float* P,float*Adj_P,
		 FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
		 FUNC_SAMPLE func_next_sample,FUNC_CMP func_cmp,const void* extra)
{
  int b=0,*bL,i,is_next,*total1,*R,*total2;
  float *bT, *count1,*count2,qT;/*qT is the successiv maxima*/
  int ncol=pdata->ncol;
  int nrow=pdata->nrow;
  int B=(*func_first_sample)(NULL);
  /*allocate the space and initialziation*/
  assert(bT=(float*)Calloc(nrow,float));
  assert(bL=(int*)Calloc(ncol,int));
  assert(count1=(float*)Calloc(nrow,float));
  memset(count1,0,sizeof(float)*nrow); 
  assert(total1=(int*)Calloc(nrow,int));
  memset(total1,0,sizeof(int)*nrow);
  assert(count2=(float*)Calloc(nrow,float));
  memset(count2,0,sizeof(float)*nrow); 
  assert(total2=(int*)Calloc(nrow,int));
  memset(total2,0,sizeof(int)*nrow);

  assert(R=(int*)Calloc(nrow,int));
   /*comuter the original t-statfirst*/

  compute_test_stat(pdata,pdata->L,T,func_stat,extra);

  /*sort the T*/  
  order_data(T,R,nrow,func_cmp);
  sort_gene_data(pdata,R);
  sort_vector(T,R,nrow);
  
  /*iteration for permutaion*/
  (*func_first_sample)(bL);

  /*changed to the orignal stat, which is monotone of t and centered*/
  is_next=1;
  b=0;
  while(is_next){
    compute_test_stat(pdata,bL,bT,func_stat,extra);
    /*deal with unajdused value first*/
    for(i=0;i<nrow;i++){
      if(T[i]==NA_FLOAT) continue;
      if(bT[i]!=NA_FLOAT){
	if((func_cmp==cmp_high)&&(bT[i]+EPSILON>=T[i])) count2[i]++;
	if((func_cmp==cmp_low)&&(bT[i]<=T[i]+EPSILON)) count2[i]++;
	if((func_cmp==cmp_abs)&&(fabs(bT[i])>=fabs(T[i])-EPSILON)) count2[i]++;
	total2[i]++;
      }
    }

    /*deal with adjusted values*/
    qT=NA_FLOAT;/*intitalize the qT*/
    for(i=nrow-1;i>=0;i--){ /*looping the row reversely*/
      if(T[i]==NA_FLOAT) continue;
        /* right now I only implements the 3 cases, which are pretty common*/
      if(func_cmp==cmp_high){
	if((bT[i]!=NA_FLOAT)&&(qT!=NA_FLOAT)&&(bT[i]>qT))
	  qT=bT[i];
	if((bT[i]!=NA_FLOAT)&&(qT==NA_FLOAT))
	  qT=bT[i];
	if((qT!=NA_FLOAT)&&(qT>=T[i]-EPSILON)) count1[i]+=1;
      }else if(func_cmp==cmp_low){
	if((bT[i]!=NA_FLOAT)&&(qT!=NA_FLOAT)&&(bT[i]<qT))
	  qT=bT[i];
	if((bT[i]!=NA_FLOAT)&&(qT==NA_FLOAT))
	  qT=bT[i];
	if((qT!=NA_FLOAT)&&(qT<=T[i]+EPSILON)) count1[i]+=1;
      }else if (func_cmp==cmp_abs) {
	if((bT[i]!=NA_FLOAT)&&(qT!=NA_FLOAT)&&(fabs(bT[i])>qT))
	  qT=fabs(bT[i]);
	if((bT[i]!=NA_FLOAT)&&(qT==NA_FLOAT))
	  qT=fabs(bT[i]);
	if((qT!=NA_FLOAT)&&(qT>=fabs(T[i])-EPSILON)) count1[i]+=1;
      }	      
      if(qT!=NA_FLOAT) total1[i]++;
    }
    b++;
    print_b(b,B,"b=");
    is_next=(*func_next_sample)(bL);
  }

  /*summarize the results*/
  /*unadjusted one*/
  for(i=0;i<nrow;i++){ 
    if(total2[i]==0) 
      P[i]=NA_FLOAT;
    else P[i]=count2[i]*1.0/total2[i];
  }
  /*adjused one*/
  for(i=0;i<nrow;i++){ 
    if(total1[i]==0) 
      Adj_P[i]=NA_FLOAT;
    else Adj_P[i]=count1[i]*1.0/total1[i];
  }
  /*enforce the montonicity*/
  for(i=1;i<nrow;i++)
    if(Adj_P[i]<Adj_P[i-1])
      Adj_P[i]=Adj_P[i-1];
  /*free the spaces*/
  Free(bT);
  Free(count1);
  Free(total1);
  Free(count2);
  Free(total2);
  Free(bL);
  Free(R);
}      
void set_seed_sampling(long int seed){
  g_random_seed=seed;
}
