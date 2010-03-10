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
#include "mt.h"
#define mtT 1
#define mtF 2
#define mtPairT 3
#define mtBlockF 4
#define mtWilcoxon 5
#define mtTequalVar 6
#define mtFixedSampling 7
typedef float (*FUNC_NUM_DENUM)(const float *, const int* ,const int,
			const float , float *, float*,const void *);
typedef void (*FUNC_CREATE)(int, int*,int);
typedef void (*FUNC_DELETE)();
typedef struct tagSAMPLING_DATA{
  FUNC_STAT fn_maxT; 
  /*the computing for maxT*, mostly needs to be standardlized*/
  FUNC_STAT fn_minP;  
  /*used to speed up the computation;mostly will be set as fn_stat;*/

  FUNC_NUM_DENUM fn_num_denum;
  /*the numerator and denumerator of maxT*/
  FUNC_STAT fn_stat;/*the centered of the original definition, 
		      no further modification
		     mostly will be fn_minP or fn_maxT, e.g.
		     in Wiloxon test, ranksum- mean, which is also fn_minP
		     in two sample t-test, it will be the t, also fn_maxT
		    */
  FUNC_CMP  fn_cmp;
  FUNC_SAMPLE fn_first;
  FUNC_SAMPLE fn_next;
  FUNC_CREATE fn_create;
  FUNC_DELETE fn_delete;
  int test;
  int is_fixed_seed;
} SAMPLING_DATA;
int type2sample(char** options,SAMPLING_DATA* sd);
int type2test(char* ptest,SAMPLING_DATA* sd);
void create_gene_data(double*d,  int*pnrow, int*pncol, int*L, double*pna,GENE_DATA* pdata,int PrintIDX)
{
  int i,j;
  pdata->nrow=*pnrow;
  pdata->ncol=*pncol;
  pdata->na=*pna;
  malloc_gene_data(pdata);
  for (j=0; j<pdata->ncol; j++) 
    pdata->L[j]=L[j];
  
  for (i=0; i<pdata->nrow; i++) {
    if(PrintIDX) sprintf(pdata->id[i],"%d",i+1); /*used for the indexes*/
    else sprintf(pdata->id[i],"0");
    for (j=0; j<pdata->ncol; j++) {
      pdata->d[i][j]=d[j*pdata->nrow+i];
      /*using the R tradition, which store the data column by column*/
    }
  }
}

void data2vec(double** data,double*d,  int nrow, int ncol)
{
  int i,j;
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      d[j*nrow+i]=data[i][j];
      /*using the R tradition, which store the data column by column*/
    }
  }
}
void get_gene_indexes(GENE_DATA* pdata, int * indexes)
{
  int i;
  for(i=0;i<pdata->nrow;i++){
    indexes[i]=atoi(pdata->id[i]);
  }
}
/*is computing fn_stat*/
void get_stat(double*d, int*pnrow, int* pncol, int*L,double *pna, float *T,char** options,int* extra)
{
  GENE_DATA data;
  SAMPLING_DATA sd;
  if(type2test(options[0],&sd)==0)
    return;
  create_gene_data(d,pnrow,pncol,L,pna,&data,0);
  compute_test_stat(&data,data.L,T,sd.fn_stat,extra);
  free_gene_data(&data);
}
void get_stat_num_denum(double*d, int*pnrow, int* pncol, int*L,double *pna, float *Tnum,float*Tdenum,char**options, int*extra)
{
  GENE_DATA data;
  SAMPLING_DATA sd;
  int i;
  if(type2test(options[0],&sd)==0)
    return;
  create_gene_data(d,pnrow,pncol,L,pna,&data,0);
  for(i=0;i<data.nrow;i++)
    (*sd.fn_num_denum)(data.d[i],data.L,data.ncol,data.na,Tnum+i,Tdenum+i,extra);
  free_gene_data(&data);
}
  
void get_maxT(double*d, int*pnrow, int* pncol, int*L,double *pna, float* T, float* P,float *adjP, int*pB, int *index,char**options, int*extra )
{
  GENE_DATA data;
  SAMPLING_DATA sd;
  if(type2sample(options,&sd)==0)
    return;
  create_gene_data(d,pnrow,pncol,L,pna,&data,1);
  (sd.fn_create)(data.ncol,data.L,*pB);
    adj_by_T(&data,T,P,adjP,sd.fn_maxT,sd.fn_first,sd.fn_next,sd.fn_cmp,extra);
  get_gene_indexes(&data,index);
  free_gene_data(&data);
  sd.fn_delete();
}
void get_minP(double*d, int*pnrow, int* pncol, int*L,double *pna, float* T, float* P,float *adjP, float* adj_lower,int*pB, int *index,char**options, int*extra) 
{
  GENE_DATA data;
  SAMPLING_DATA sd;
  if(type2sample(options,&sd)==0)
    return;
  create_gene_data(d,pnrow,pncol,L,pna,&data,1);
  Rprintf("B=%d\n",*pB);
  sd.fn_create(data.ncol,data.L,*pB);
  adj_pvalue_quick(&data,T,P,adjP,adj_lower,sd.fn_minP,sd.fn_maxT,sd.fn_first,sd.fn_next,sd.fn_cmp,extra);
  get_gene_indexes(&data,index);
  free_gene_data(&data);
  sd.fn_delete();
}
void get_samples_T(float*V, int* pn,int* L,float* T,float *pna,int* pB,char**options, int*extra) 
{
  int n=*pn;
  int B=*pB;
  SAMPLING_DATA sd;
  if(type2sample(options,&sd)==0)
    return;
  sd.fn_create(n,L,B);
  get_all_samples_T(V,n,T,*pna,sd.fn_maxT,sd.fn_first,sd.fn_next,(void*)extra);
  sd.fn_delete();
}
void get_samples_P(float*V, int* pn,int* L,float* P,float *pna,int* pB,char**options, int*extra)
{
  int n=*pn;
  int B=*pB;
  SAMPLING_DATA sd;
  if(type2sample(options,&sd)==0)
    return;
 sd.fn_create(n,L,B);
  get_all_samples_P(V,n,P,*pna,sd.fn_minP,sd.fn_first,sd.fn_next,sd.fn_cmp,(void*)extra);
  sd.fn_delete();
}

void get_sample_labels(int*pn,int*L,int*pB,int* S,char**options, int*extra)
{
  int n=*pn;
  int B=*pB;
  int is_next=1;
  int nb=0;
  int i;
  SAMPLING_DATA sd;
  if(type2sample(options,&sd)==0)
    return;
  sd.fn_create(n,L,B);
  sd.fn_first(L);
  while(is_next){
    for(i=0;i<n;i++)
      S[nb+i]=L[i];
    nb+=n;
    is_next=sd.fn_next(L);
    
  }
  sd.fn_delete();
}
int type2test(char* ptest,SAMPLING_DATA* sd)
{
  int test=0;
  if(strcmp(ptest,"t")==0){
    test=mtT;
    sd->fn_stat=two_sample_tstat;
    sd->fn_maxT=sd->fn_stat;
    sd->fn_minP=sd->fn_stat;
    sd->fn_num_denum=two_sample_tstat_num_denum;
  }else  if(strcmp(ptest,"f")==0){
    test=mtF;
    sd->fn_stat=Fstat;
    sd->fn_num_denum=Fstat_num_denum;
    sd->fn_maxT=sd->fn_stat;
    sd->fn_minP=sd->fn_stat;
  }else  if(strcmp(ptest,"pairt")==0){
    test=mtPairT;
    sd->fn_stat=sign_tstat;
    sd->fn_num_denum=sign_tstat_num_denum;
    sd->fn_maxT=sd->fn_stat;
    sd->fn_minP=sign_sum;/*changed to montone*/
  }else  if(strcmp(ptest,"blockf")==0){
    test=mtBlockF;
    sd->fn_stat=Block_Fstat;
    sd->fn_num_denum=Block_Fstat_num_denum;
    sd->fn_maxT=sd->fn_stat;
    sd->fn_minP=sd->fn_stat;
  }
  else if(strcmp(ptest,"wilcoxon")==0){
    test=mtWilcoxon;
    sd->fn_stat=Wilcoxon_T;
    sd->fn_num_denum=Wilcoxon_num_denum;
    sd->fn_maxT=sd->fn_stat;/*changed to normalize*/
    sd->fn_minP=Wilcoxon_stat;
  }else if(strcmp(ptest,"t.equalvar")==0){
    test=mtTequalVar;
    sd->fn_stat=two_sample_t1stat;
    sd->fn_num_denum=two_sample_t1stat_num_denum;
    sd->fn_maxT=sd->fn_stat;
    sd->fn_minP=ave_diff;/*changed to montone*/
  }else
    return 0;
  sd->test=test;
  return 1;
}
int type2sample(char** options,SAMPLING_DATA* sd)
{
  char *ptest,*pfixed_seed,*pside;
  int test=0;
  int is_fixed_sampling=0;
  int side=-2;
  /************************/  
  ptest=options[0];
  pside=options[1];
  pfixed_seed=options[2];

  /***********************/
  type2test(ptest,sd);
  test=sd->test;

  /***********************/
  if(strcmp(pside,"upper")==0)
    side=1;
  if(strcmp(pside,"lower")==0)
    side=-1;
  if(strcmp(pside,"abs")==0)
    side=0;
  sd->fn_cmp=side2cmp(side);
  /**************/
  if(strcmp(pfixed_seed,"y")==0)
    is_fixed_sampling=mtFixedSampling;
  else
    is_fixed_sampling=0;
  sd->is_fixed_seed=is_fixed_sampling;
  /***************/
  switch(test){
  case mtT: case mtF: case mtWilcoxon: case mtTequalVar:
    if(is_fixed_sampling){
      sd->fn_first=first_sample_fixed;
      sd->fn_next=next_sample_fixed;
      sd->fn_create=create_sampling_fixed;
      sd->fn_delete=delete_sampling_fixed;
    }else{
      sd->fn_first=first_sample;
      sd->fn_next=next_sample;
      sd->fn_create=create_sampling;
      sd->fn_delete=delete_sampling;
    }
    break;
  case mtPairT: 
    if(is_fixed_sampling){
      sd->fn_create=create_sampling_pairt_fixed;
      sd->fn_delete=delete_sampling_pairt_fixed;
      sd->fn_first=first_sample_pairt_fixed;
      sd->fn_next=next_sample_pairt_fixed;
    }else{
	 sd->fn_create=create_sampling_pairt;
	 sd->fn_delete=delete_sampling_pairt;
	 sd->fn_first=first_sample_pairt;
	 sd->fn_next=next_sample_pairt;
    }
    break;
  case mtBlockF:/*have not implemented the solutuion for storing the permutation yet, as it is very memory instensive*/
    sd->fn_create=create_sampling_block;
    sd->fn_delete=delete_sampling_block;
    sd->fn_first=first_sample_block;
    sd->fn_next=next_sample_block;
    break;
  default:
    fprintf(stderr,"Can not recogize the parameter\n");
    return 0;
  }
  return 1;
}
      
    

    
/*test*/
/*main()
{
  #define N 6
  #define NUMB 6
  int n=N;
  int L[N]={0,0,1,1,2,2};
  int B=NUMB;
  int S[NUMB*N];
  int i;
  get_sample_labels(&n,L,&B,S);
  for(i=0;i<B;i++)
    print_narray(stderr,S+i*n,n);
}*/


  



