/*****************************************************************/
/*          some options                                         */
/*****************************************************************/
/*#define USEDOUBLE*/
/*to transport to R, fix it to be a double*/
/*for windows'R*/
/*#define WINDOWS is used to build R package for windows system*/
#include "R.h"
#ifdef WINDOWS
/*#include "R_ext\print.h"*/
#define fprintf win_print
void win_print(FILE*, char*,...);
#endif
/* these options can be incorporated in the command argments if you like*/  
extern int myDEBUG;  /*this variable needs to be declared in the main file*/
extern long int g_random_seed;
#define PRINT_VAR_NUM 10  /*the ouput of the length of each line in the debug session
			  when myDEBUG is 1, you need to try small data and small simulations
			  otherwise, it just has too much printing*/
#define PROMPT_LEN 0        
/*  if PROMPT_LEN>1, then we'll use PROMPT_LEN+1 to report the progess
    in permutations, otherwise, we'll report it when finish 1% of 
    permutations, this only applies to the permutations in 
    get the unadjusted p-values. For calculating adjusted p-values, it always 
    prompt after finish every genes, as the total number of
    genes is typically small, around 6000*/ 

#include <float.h>

#define MAX_ID 40                  /*the max number of characters allowed for ID*/
#define MAX_WARN 256               /*the max of chars allowed in the warning message*/
#define NA_FLOAT FLT_MAX          /*the default NA representation for float number*/

#define NA_DATA   1e30             /*the default NA representation for a gene value*/


/*using the double data*/
#ifdef USEDOUBLE
#define EPSILON (120*DBL_EPSILON)
#define float double
#else
#define EPSILON (12*FLT_EPSILON)
#endif


typedef struct tagGENE_DATA{
  char** id; /*the gene index id*/
  float** d; /*the gene values matrix, mxn*/
  float na;  
  int nrow; /*nrow is the number of the genes*/
  int ncol; /*ncol is the number of the experiments*/
  int* L;   /*the status labelling of each experiment*/
  char name[MAX_ID];/*the name of the status*/
}GENE_DATA;

typedef int (*FUNC_SAMPLE)(int *);
typedef int (*FUNC_CMP)(const void*,const void*);
typedef float (*FUNC_STAT)(const float*,const int*,const int,const float,const void*);
/********************************************************************************/
/*               multiple testing                                              */
/********************************************************************************/
void get_all_samples_P(float* V, int n,float* P,float na, 
		       FUNC_STAT func_stat,FUNC_SAMPLE first_sample, 
		       FUNC_SAMPLE next_sample,FUNC_CMP func_cmp,const void* extra);
void get_all_samples_T(float* V, int n,float* T,float na, 
		       FUNC_STAT func_stat,FUNC_SAMPLE first_sample, 
		       FUNC_SAMPLE next_sample,const void* extra);

void adj_pvalue_quick(GENE_DATA* pdata,float*T, float* P, 
		      float* Adj_P,float* Adj_Lower,
		      FUNC_STAT func_stat,FUNC_STAT func_stat_maxT,FUNC_SAMPLE first_sample, 
		      FUNC_SAMPLE next_sample,FUNC_CMP func_cmp,const void* extra);
void  get1pvalue(GENE_DATA* pdata,int* L,float* T,float* P,
		 FUNC_STAT func_stat,FUNC_SAMPLE first_sample, 
		 FUNC_SAMPLE next_sample,FUNC_CMP func_cmp,const void* extra);
void  adj_by_T(GENE_DATA* pdata,float* T,float* P,float*Adj_P,
		 FUNC_STAT func_stat,FUNC_SAMPLE func_first_sample, 
	       FUNC_SAMPLE func_next_sample,FUNC_CMP func_cmp,const void* extra);	       
		 
/********************************************************************************/
/*              processing with the gene_data                                   */
/********************************************************************************/
void read_infile(char *filename,GENE_DATA *pdata);
void write_outfile(FILE* fp,GENE_DATA* pdata,float*T, float*P,float*Adj_P,float* Adj_Lower);
void malloc_gene_data(GENE_DATA* pdata);
void free_gene_data(GENE_DATA* pdata);
void print_gene_data(GENE_DATA* pdata);
void sort_gene_data(GENE_DATA* pdata,int*R);
void sort_vector(float* V,int*R,int n);

/********************************************************************************/
/*               sampling good for two sample t and F-stat                      */
/********************************************************************************/
void create_sampling(int n,int*L,int B);
void delete_sampling();
int first_sample(int *L);
int next_sample(int* L);

void create_sampling_fixed(int n,int*L,int B);
void delete_sampling_fixed();
int first_sample_fixed(int *L);
int next_sample_fixed(int* L);
void set_seed_sampling(long int seed);

void create_sampling_block(int n,int*L,int B);
void delete_sampling_block();
int first_sample_block(int *L);
int next_sample_block(int* L);

void create_sampling_pairt(int n,int*L,int B);
void delete_sampling_pairt();
int first_sample_pairt(int *L);
int next_sample_pairt(int* L);

void create_sampling_pairt_fixed(int n,int*L,int B);
void delete_sampling_pairt_fixed();
int first_sample_pairt_fixed(int *L);
int next_sample_pairt_fixed(int* L);
/********************************************************************************/
/*            data_sorting                                                      */
/********************************************************************************/
void order_mult_data(int* R,int n,int k,...);
void order_data(float* V,int*R,int n,FUNC_CMP func_cmp);

int cmp_high(const void *v1, const void *v2);
int cmp_low(const void *v1, const void *v2);
int cmp_abs(const void *v1, const void *v2);


/*micesslay functions*/
void print_farray(FILE* fh,float* p_arr,int n);
void print_narray(FILE*fh,int* p_arr,int n);

/********************************************************************************/
/*            common used statistics                                            */
/********************************************************************************/
void compute_test_stat(GENE_DATA* pdata, int* L,float* T,
		       FUNC_STAT func_stat,const void* extra);
float two_sample_tstat(const float *Y, const int* L,const int n,
		       const float na,const void *extra); 
float two_sample_tstat_num_denum(const float *Y, const int* L,const int n, 
			       const float na,float* num, float* denum,const void* extra);

/*t1stat is dealing with two sample t-statistics with equal variance*/
float ave_diff(const float *Y, const int* L,const int n, const float na,const void* extra); /*used to speed up the minP as ave_diff is monotone of the t1stat*/
float two_sample_t1stat(const float *Y, const int* L,const int n,
		       const float na,const void *extra); 
float two_sample_t1stat_num_denum(const float *Y, const int* L,const int n, 
			       const float na,float* num, float* denum,const void* extra);

/* Wilkoxon test*/
float Wilcoxon_stat(const float *Y, const int* L,const int n, const float na,const void* extra); /*T-ET, where ET=1/2 n_0(n_0|+n_1+1), T is the sum of rank*/
float Wilcoxon_T(const float *Y, const int* L,const int n,
		       const float na,const void *extra); 
/* is computing (T-ET)/var(T), wher var(T)=1/12*n_0*n_1*(n_0+n_1+1) */
float Wilcoxon_num_denum(const float *Y, const int* L,const int n, 
			       const float na,float* num, float* denum,const void* extra);

float sign_sum(const float *Y, const int* L,const int n, const float na,const void* extra);
float sign_tstat_num_denum(const float *Y, const int* L,const int n,
			 const float na, float *num, float*denum,const void *extra);
float sign_tstat(const float *Y, const int* L,const int n, 
		 const float na, const void* extra);
float Fstat_num_denum(const float *Y, const int* L,const int n, const float na,
		      float *num, float*denum,const void* extra);
float Fstat(const float *Y, const int* L,const int n,
	    const float na,const void* extra);
float Block_Fstat(const float *Y, const int* L,const int n, 
		  const float na,const void* extra);
float Block_Fstat_num_denum(const float *Y, const int* L,const int n,
			    const float na, float *num, float*denum,const void* extra);
/********************************************************************************/
/*           some useful tools                                                  */
/********************************************************************************/

int bin2int(int*V,int n);/*integrate the bits as an integer*/  
void int2bin(int x,int*V,int n);/*decompose an integar as bits stored in V*/
int bincoeff(int n, int k);
double logbincoeff(int n, int k);
double logfactorial(int n, int k);
void init_label(int n, int k, int*nk, int*L);
void init_label_block(int *L, int n,int m);
int next_label_block(int* L, int n, int m);
void sample_block(int *L, int n,int m);
void sample2label(int n, int k, int* nk,int *permun, int*L);
void label2sample(int n, int k, int* nk,int*L,int *permun);
int next_label(int n, int k, int* nk, int*L);
int next_lex(int* A, int n, int k);
void A2L(int* A,int* L,int n,int k);
FUNC_CMP side2cmp(int side);
void sample(int *V, int n, int m);
float get_rand();
void set_seed(long int seed);
int next_mult_permu(int* V, int n, int k, int* nk);
int next_two_permu(int* V, int n, int k);
int next_permu(int*V,int n);
void data2vec(double** data,double*d,  int nrow, int ncol);
void set_seed_sampling(long int seed);

void get_maxT(double*, int*, int*, int*, double*, float*, float*,
	      float*, int*, int*, char**, int*);
