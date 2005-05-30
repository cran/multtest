#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "math.h"
#include "assert.h"
#include "string.h"
#include "mt.h"

/*This file is used to collect some useful statistics functions*/
/********************************************************************************/
/*                        two_sample_tstat                                      */
/********************************************************************************/
    /* Computes the value of the two sample t-statistic, allowing for missing values. 
       Missing values are represented by na.  (At least two values per group should 
       be present.) 
       if return == NA_FLOAT, then it has some problems
       to calculate the t-stat,such as variance is 0, 
       or the count of one class is less than 2
      Y: the vector of one gene across experiments
      n: the number of experiments
      L: the class labelling of each experiments
      na: the NA representation of gene values.
      extra: the additional information, not used here
    */
float two_sample_tstat(const float *Y, const int* L,const int n, const float na,const void *extra) 
{
  float num,denum,res;
  res=two_sample_tstat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num/denum;
}

float two_sample_tstat_num_denum(const float *Y, const int* L,const int n, const float na,float* num, float* denum,const void* extra) 
{
  float mean_na[2]={0,0},ss_na[2]={0,0},devi;
  float c0,c1;
  int i,count[2]={0,0},class;
  /*compute the mean and count first*/
  /* count is the number of objects in each class*/
  for (i=0; i<n; i++) {
    if (Y[i]==na) 
      continue;
    class=L[i];
    mean_na[class] += Y[i];
    count[class]++; 
  }
 
  mean_na[0]/=(count[0]*1.0);
  mean_na[1]/=(count[1]*1.0);

  /*compute the variance in each group*/
  for (i=0; i<n; i++) {
    if (Y[i]==na) 
       continue;
    class=L[i];
    devi=(Y[i]-mean_na[class]);
   ss_na[class] += devi*devi;
  }
  /* if(ss_na[0]==0) return NA_FLOAT;
  if(ss_na[1]==0) return NA_FLOAT;
  */
  if(ss_na[0]+ss_na[1]<EPSILON) return NA_FLOAT;
  c0=(count[0]*(count[0]-1));
  c1=(count[1]*(count[1]-1));
  *num=mean_na[1]-mean_na[0];
  *denum=sqrt(ss_na[0]/c0+ss_na[1]/c1);
  return 1;
}
float ave_diff(const float *Y, const int* L,const int n, const float na,const void* extra)
{
  float mean_na[2]={0,0};
  int i,count[2]={0,0},class;
  /*compute the mean and count first*/
  /* count is the number of objects in each class*/
  for (i=0; i<n; i++) {
    if (Y[i]==na) 
      continue;
    class=L[i];
    mean_na[class] += Y[i];
    count[class]++; 
  }

  mean_na[0]/=(count[0]*1.0);
  mean_na[1]/=(count[1]*1.0);
  if((count[0]==0)||(count[1]==0)) return NA_FLOAT;
  return mean_na[1]-mean_na[0];
}
  
float two_sample_t1stat(const float *Y, const int* L,const int n, const float na,const void *extra) 
{
  float num,denum,res;
  res=two_sample_t1stat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num/denum;
}
float two_sample_t1stat_num_denum(const float *Y, const int* L,const int n, const float na,float* num, float* denum,const void* extra) 
{
  float mean_na[2]={0,0},ss_na[2]={0,0},devi;
  float c0,c1;
  int i,count[2]={0,0},class;
  /*compute the mean and count first*/
  /* count is the number of objects in each class*/
  for (i=0; i<n; i++) {
    if (Y[i]==na) 
      continue;
    class=L[i];
    mean_na[class] += Y[i];
    count[class]++; 
  }

  mean_na[0]/=(count[0]*1.0);
  mean_na[1]/=(count[1]*1.0);

  /*compute the variance in each group*/
  for (i=0; i<n; i++) {
    if (Y[i]==na) 
       continue;
    class=L[i];
    devi=(Y[i]-mean_na[class]);
   ss_na[class] += devi*devi;
  }
  /*  if(ss_na[0]==0) return NA_FLOAT;
  if(ss_na[1]==0) return NA_FLOAT;
  */
  if(ss_na[0]+ss_na[1]<EPSILON) return NA_FLOAT;
  c0=1/(count[0]*1.0)+1/(count[1]*1.0);
  c1=count[0]+count[1]-2.0;
  *num=mean_na[1]-mean_na[0];
  *denum=sqrt((ss_na[0]+ss_na[1])*c0/c1);
  return 1;
}

float Wilcoxon_stat(const float *Y, const int* L,const int n, const float na,const void* extra)
{
  int i,count=0,n1=0;
  float s=0;
  for(i=0;i<n;i++){
    if(Y[i]==na) continue;
    if(L[i]){
      s+=Y[i];
      n1++;
    }
    count++;
  }
  s-=(1+count)*n1/2.0;
  return s;
}
      
float Wilcoxon_T(const float *Y, const int* L,const int n,
		       const float na,const void *extra)
{
  float num,denum,res;
  res=Wilcoxon_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num/denum;
} 
float Wilcoxon_num_denum(const float *Y, const int* L,const int n, 
			       const float na,float* num, float* denum,const void* extra)
{
  int i,count=0,n1=0;
  float s=0;
  for(i=0;i<n;i++){
    if(Y[i]==na) continue;
    if(L[i]){
      s+=Y[i];
      n1++;
    }
    count++;
  }
  *num=s-(1+count)*n1/2.0;
  *denum=sqrt(n1*(count-n1)*(count+1)/12.0);
  if((*denum)<EPSILON) return NA_FLOAT;
  return 1; 
}
/* The signed sum of Y[i], i.e. compute the sum of Y[i]*L[i], where L[i] is 1 or 0, 
1 is for treatment, 0 is for control,
It is used in the paired t-stat, this test is monotone of the paired t-stat in the 
2xB block design, where we have B blocks, and each block has two experiments, while Y[i]
is already the difference within one block, we don't consider the na to less complicate
the function*/
float sign_sum(const float *Y, const int* L,const int n, const float na,const void* extra)
     /* extra not used*/
{
  float ret=0;
  int i;
  int count=0;
  for(i=0;i<n;i++){
    if(Y[i]==0) continue;
    if(L[i]){
      ret+=Y[i];
    }else{
      ret-=Y[i];
    };
    count++;
  }
  return ret;
}
/*the function is compute the one sample t-statistics,
used to compute the paired t-stat for 2xB block design, where while Y[i]
is already the difference within one block,*/ 
float sign_tstat_num_denum(const float *Y, const int* L,const int n, const float na,float *num, float*denum,const void *extra)
{
  float ss=0;
  float mean=0;
  float devi;
  int i,count;
  count=0;
  for(i=0;i<n;i++){
    if(Y[i]==na) continue;
    if(L[i]){
      mean+=Y[i];
    }else{
      mean-=Y[i];
    };
    count++;
  }
  mean/=(count*1.0);
  for(i=0;i<n;i++){
    if(L[i]){
      devi=Y[i]-mean;
    }else{
      devi=-Y[i]-mean;
    };
    ss+=(devi * devi);
  }
  *num=mean;
  *denum=sqrt(ss/(count*(count-1.0)));
  if((*denum)<EPSILON) return NA_FLOAT;
  return 1;
}
float sign_tstat(const float *Y, const int* L,const int n, const float na, const void* extra)
     /*extra not used*/
{
  float num,denum,res;
  res=sign_tstat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  return num/denum;
}  

/*compute the F-stat for k samples, where L[i] is the labelling of object i,
note for F-test, (int*)*extra has the information of the number of groups*/
float Fstat(const float *Y, const int* L,const int n, const float na,const void* extra)
{
  float num,denum,res;
  res=Fstat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  if(denum<EPSILON) return NA_FLOAT;
  return num/denum;
}
float Fstat_num_denum(const float *Y, const int* L,const int n, const float na,float *num, float*denum,const void* extra)
{
  float wss=0,bss=0;/*within sum square, and between sum square*/
  float mean=0,*meani,*ssi;/*the mean for for the whole and the mean for group i*/
  float dev;
  int i,class,k,*ni,N=0;/* ni the number of objects group i, 
		   k is the number of groups,
		   N is the total number of validate objects*/
  k=*(int*)extra;
  meani=(float *)Calloc(k,float);
  ssi=(float *)Calloc(k,float);
  ni=(int *)Calloc(k,int);
  for(i=0;i<k;i++){
    meani[i]=0;/*intialize to be zero*/
    ssi[i]=0;
    ni[i]=0;
  }
  for(i=0;i<n;i++){
    if(Y[i]==na) 
      continue;
    class=L[i];
    meani[class]+=Y[i];
    ni[class]++;
    N++;
    mean+=Y[i];
  }
  /*summarize the data*/
  mean/=(N*1.0);
  for(i=0;i<k;i++){
    meani[i]/=(ni[i]*1.0);/*get the mean for each group*/
  }
  /*compute bss and wss*/
  for(i=0;i<n;i++){
    if(Y[i]==na) 
      continue;
    class=L[i];
    dev=Y[i]-meani[class];
    ssi[class]+=dev*dev;
  }
  for(i=0;i<k;i++){
    /*    if(ssi[i]==0) {
	     ret=NA_FLOAT;
	     break;
	     }; */ /*each group needs to have non zero variance*/
    wss+=ssi[i];
    dev=meani[i]-mean;
    bss+=dev*dev*ni[i];
  }
  /*summarize the num and denum*/
  *num=bss/(k-1.0);
  *denum=wss/(N-k-0.0);
  Free(meani);
  Free(ni);
  Free(ssi);
  return 1;
} 
/*compute the Block F-stat for k samples, where L[i] is the labelling of object i,
note for block F-test, (int*)*extra has the information of the number of Blocks*/
/* currently, we don't implement the na, as we need the rectangular design*/
float Block_Fstat(const float *Y, const int* L,const int n, const float na,const void* extra)
{
  float num,denum,res;
  res=Block_Fstat_num_denum(Y,L,n,na,&num,&denum,extra);
  if(res==NA_FLOAT) return NA_FLOAT;
  if(denum<EPSILON) return NA_FLOAT;
  return num/denum;
}
float Block_Fstat_num_denum(const float *Y, const int* L,const int n, const float na,float *num, float*denum,const void* extra)
{
  float wss=0,bss=0;/*within sum square, and between sum square*/
  float mean=0,*meani,*meanj;
  /*the mean for for the whole and the mean for group i,j: i is for block, j is for treatment*/
  float dev;
  int treat,block,i,j,B,m,h;/* ni the number of objects group i, 
		   k is the number of groups,*/

  m=*(int*)extra;
  B=n/m;
  if(B*m!=n){
    fprintf(stderr,"The design is not balanced as B(%d)xm(%d)!=n(%d)\n",B,m,n);
    return NA_FLOAT;
  }
  meani=(float *)Calloc(B,float);
  meanj=(float *)Calloc(m,float);
  for(i=0;i<B;i++){
    meani[i]=0;/*intialize to be zero*/
    for(j=0;j<m;j++){
      meani[i]+=Y[j+m*i];
    }
  }
  for(j=0;j<m;j++){
    meanj[j]=0;/*intialize to be zero*/
  }
  for(h=0;h<n;h++){
    treat=L[h];
    meanj[treat]+=Y[h];
    mean+=Y[h];
  }
   
  /*summarize the data*/
  mean/=(n*1.0);
  for(i=0;i<B;i++){
    meani[i]/=(m*1.0);/*get the mean for each block*/
  }
  for(j=0;j<m;j++)
    {
    meanj[j]/=(B*1.0);/*get the mean for each treatments*/
  }
  /*compute bss and wss*/
  for(i=0;i<n;i++){
    block=i/m;    
    dev=Y[i]-meani[block];
    treat=L[i];
    dev-=meanj[treat]-mean;
    wss+=dev*dev;
  }
  for(j=0;j<m;j++){
    dev=meanj[j]-mean;
    bss+=dev*dev*B;
  }
  /*summarize the num and denum*/
  *num=bss/(m-1.0);
  *denum=wss/((m-1.0)*(B-1.0));
  Free(meani);
  Free(meanj);
  return 1;
} 
void int2bin(int r,int*V,int n)
{ int i;
  for(i=n-1;i>=0;i--){
    V[i]=r&1;
    r>>=1; /*divide by 2 so that we can look-at the next digit*/
  }
}
int bin2int(int*V,int n)
{ int i,ret=0;
  for(i=0;i<n-1;i++){
    ret+=V[i];
    ret<<=1; /*multiply by 2 so that we can look at the next digit*/
  }
  ret+=V[n-1];
  return ret;
}
/********************************************************************************/
/*                      bincoeff                                                */
/********************************************************************************/
/*Descriptions: return the binomial coefficient of n choosing k*/
int bincoeff(int n, int k)
{
  float f=n;
  int i;
  for(i=1;i<k;i++)
    f*=(n-i)/(i+1.0);
  return (int)(f+0.5);
}
/*compute the log of the coefficient*/
double logbincoeff(int n, int k)
{
  double f=log(n);
  int i;
  for(i=1;i<k;i++)
    f+=log((n-i)/(i+1.0));
  return f;
}
double logfactorial(int n, int k)
{
  double f=log(n);
  int i;
  for(i=1;i<k;i++)
    f+=log((n-i)*1.0);
  return f;
}
/********************************************************************************/
/*                 A2L                                                      */
/********************************************************************************/
/*Descriptions:
    Assume we have n objects, of which k of them are labeled with 0, the rest of them
    are labeling with 1. A is the subset of the objects which have label 0, L is the 
    labelling of each object. This function transforms A to Label
*/ 
void A2L(int* A,int* L,int n,int k)
{
  int i;
  for(i=0;i<k;i++)
    L[i]=0;
  for(i=k+1;i<n;i++)
    L[i]=1;
}     

/********************************************************************************/
/*                     next_lex                                                 */
/********************************************************************************/
   /* Given a list a of k numbers a_0<a_1<...<a_(k-1) with a_i in {0, 1,..., 
      n-1}, returns the list b: b_0<b_1<...<b_(k-1) which immediately follows a in 
      lexicographic order. It can be used to determine all subsets of size k in
      {0, 1,..., n-1}. 

      note array is in A and after this function, is modified
      and store the new array B back in array A
     
   */
int next_lex(int* A, int n, int k)
{
  int l=k-1, s=n-1,i,old;/*l is for the location of A to increase*/

  /*look for a location to increase*/  
  while (A[l]==s &&l>=0) {
    l--;
    s--;
  }
  if(l<0) {
    if (myDEBUG)
      {
	fprintf(stderr,"%s%s","We've achieved the maximum permutation already\n",
		"We can not find the next one in next_lex\n");
      }

    return 0;/*note we can not generate the next permutations*/
  }
  /*we increase every number by 1*/
  old=A[l];
  for(i=l;i<k;i++)
    A[i]=old+i-l+1;
  return 1;
}

/*intialize the label such that the first nk[0] is labelled as 0, etc.*/
void init_label(int n, int k, int*nk, int*L)
{
  int l,s,j;
  s=0;/*s is for starting*/
  for(l=0;l<k;l++){
    for(j=0;j<nk[l];j++){
      L[s]=l;
      s++;
    }
   }
}
void sample2label(int n, int k, int* nk,int *permun, int*L)
{
  int l,s,j;
  s=0;/*s is for starting*/
  for(l=0;l<k;l++){
    for(j=0;j<nk[l];j++){
      L[permun[s]]=l;
      s++;
    }
  }
}
void label2sample(int n, int k, int* nk,int*L,int *permun)
{ 
  int l,j;
  int *s;/*s is for starting*/

  /*initialize the beginning*/
  assert(s=(int*)Calloc(k,int));
  s[0]=0;
  for(l=1;l<k;l++){
    s[l]=s[l-1]+nk[l-1];
  }
  for(j=0;j<n;j++){
    l=L[j];
    permun[s[l]]=j;
    s[l]++;
  }
  Free(s);
}
  
int next_label(int n, int k, int* nk, int*L)
{
  int *permun,ret;
  assert(permun=(int*)Calloc(n,int));
  label2sample(n,k,nk,L,permun);
  ret=next_mult_permu(permun,n,k,nk);
  sample2label(n,k,nk,permun,L);
  Free(permun);
  return ret;
}

/*V has n elements, the first k elements(referred as array A) are ordered, 
  and the rest n-k elements (referred as array B)are ordered increasing order, 
  this is the problem to list all binomial permutation from n choosing k.
  after the next_permu, if it's the last one, return false, and the whole array
  V will be ordered, just swap the array A and B. otherwise, V will be next permutation and keep the same order, the algorithm has computation O(N)*/

int next_two_permu(int* V, int n, int k)
{
  int i,j;
  int old,maxb=V[n-1];
  int* A=V;
  int* B=V+k;
  int* tempV,*cpyV;
  assert(tempV=(int*)Calloc(n,int));

  i=k-1;
  while(i>=0&& A[i]>maxb){
    i--;
  }
  /*there's no next_permu as all the elements of array A 
    is greater than array 
    B*/
  if(i<0){
    /*rearrange the output so that V be ordered for the whole array.*/
    memcpy(tempV,B,sizeof(int)*(n-k));
    memcpy(tempV+(n-k),A,sizeof(int)*k);
    /*using the tempV to swap the array A and array B*/
    memcpy(V,tempV,sizeof(int)*n);
    /*coppying back to V*/
    Free(tempV);
    return 0;
  }
  /*else to find the next permutation*/
  /*first to find how many elements in B are between A[i] and A[i+1]*/
  j=n-k-2;
  old=A[i];
  while(j>=0&&(B[j]>old)){
    j--;
  }
  /*keep the original A[0..(i-1)] elements to tempV*/
  memcpy(tempV,A,sizeof(int)*i);
  /*keep the original B[0..j] elements to tempV+k*/
  if(j+1>0)
    memcpy(tempV+k,B,sizeof(int)*(j+1));
  /*copy the (k-i) elements from array 
    (A[i]<)B[j+1],...B[n-k-1],A[i+1],..A[k-1]*/
  /*copy the ((n-k)-(j+1)) elements from array 
    B[j+1],...B[n-k-1],..,A[i+1],..A[k-1]*/
  /*construct the array B[j+1],...B[n-k-1],A[i+1],..A[k-1]*/
  assert(cpyV=(int*)Calloc(n,int));
  memcpy(cpyV,B+j+1,sizeof(int)*((n-k)-(j+1)));
  if(k>(i+1))
    memcpy(cpyV+(n-k)-(j+1),A+i+1,sizeof(int)*(k-(i+1)));
  memcpy(tempV+i,cpyV,sizeof(int)*(k-i));
  tempV[k+j+1]=A[i];
  if((n-k)>(j+2))
    memcpy(tempV+k+j+2,cpyV+(k-i),sizeof(int)*((n-k)-(j+2)));
  /*copy back to V*/
  memcpy(V,tempV,sizeof(int)*n);
  Free(cpyV);
  Free(tempV);
  return 1;
}

/*the next_permutation for multiple classes*/
int next_mult_permu(int* V, int n, int k, int* nk)
{
  int olds,s,l;/*s is for starting location*/
  int next=0;
  /*initialize the begining*/
  s=nk[0];
  for(l=1;l<k;l++){
    olds=s;
    s+=nk[l];
    next=next_two_permu(V,s,olds);
    if(next) return 1;
  }
  return 0;/*we couldn't find the next permutation*/
}
FUNC_CMP side2cmp(int side)
{
  FUNC_CMP func_cmp;
  if(side==0){
    func_cmp=cmp_abs;
  }else if(side==-1){
    func_cmp=cmp_low;
  }else{
    func_cmp=cmp_high;
  }
  return func_cmp;
}
/*V[0],V[1],..,V[n] is a permutation of h_1<h_2<...<h_n, where h_i may be i in
most case, but here to be general in order to be mos use*/
/*it returns the next permutation in V*/ 
int next_permu(int*V,int n) /* n has to be at least 2*/
{
  int i,j,old,*cpyV,l;
  i=n-2;
  while(i>=0){
    if(V[i]<V[i+1]) break;
    i--;
  }
  if(i<0){
    if (myDEBUG){
   	fprintf(stderr,"%s%s","We've achieved the maximum permutation already\n",
		"We can not find the next one-in next_permu\n");
    }
    return 0;/*note we can not generate the next permutations*/
  }
  /*find the location of j, V[i]<V[i+1]>...>V[n-1]*/
  /*i.e. V[n-1]<V[n-2]<...V[j+1]<V[i]=old<V[j]<...V[i+1]*/
  old=V[i];
  j=n-1;
  while(j>i){
    if(V[j]>old) break;
    j--;
  }
  assert(cpyV=(int*)Calloc(n,int));
  memcpy(cpyV,V,sizeof(int)*n);
  V[i]=cpyV[j];
  cpyV[j]=old;
  for(l=i+1;l<n;l++){
    V[l]=cpyV[n+i-l];
  }
  Free(cpyV);
  return 1;
}

/********************************************************************************/
/*                      print_farray                                             */
/********************************************************************************/
/*Description:
      used to print an array with n elements.
      we can specify*/

void print_narray(FILE* fh,int* p_arr,int n)
{
  int i;
  for(i=0;i<n;i++){
    fprintf(fh," %7d ",p_arr[i]);
    if((PRINT_VAR_NUM)&&((i+1)%PRINT_VAR_NUM==0))
      fprintf(fh,"\n");
  }
  fprintf(fh,"\n");
}
void print_farray(FILE* fh,float* p_arr,int n)
{
  int i;
  for(i=0;i<n;i++){
    fprintf(fh," %9g ",p_arr[i]);{
      if((PRINT_VAR_NUM)&&((i+1)%PRINT_VAR_NUM==0))
	fprintf(fh,"\n");
    }
  }
  fprintf(fh,"\n");
}


  
/********************************************************************************/
/*                       read_infile                                            */
/********************************************************************************/
/*read the file into the struct GENE_DATA.
  1). The first row of the file should contain the status*
  2)  The space of pdata should be allocated already by malloc_gene_data(),
  otherwise it might have loss of memory problem
*/
void read_infile(char *filename,GENE_DATA *pdata) {
  FILE *fh;
  int i, j;
  double ftemp;
  assert(fh=fopen(filename,"r"));
  /*read the labelling first*/
  assert(fscanf(fh, "%s", pdata->name));
  for (j=0; j<pdata->ncol; j++) 
    assert(fscanf(fh, "%d", pdata->L+j));
  
/*read the mxn matrix of the gene values data*/
  for (i=0; i<pdata->nrow; i++) {
    /*read gene id and the first gene values*/
    assert(fscanf(fh, "%s", pdata->id[i]));
    /*read the rest of it*/
    for (j=0; j<pdata->ncol; j++) {
      /*deal with the double data*/
      assert(fscanf(fh, "%lg",&ftemp));
      pdata->d[i][j]=ftemp;
    }
  }
  fclose(fh);
}

/********************************************************************************/
/*             print_gene_data                                                  */
/********************************************************************************/
/*print the gene_data to stderr, useful in debug*/
void print_gene_data(GENE_DATA* pdata)
{
  int i, j;
  for (i=0; i<pdata->nrow; i++){
    fprintf(stderr,"%20s", pdata->id[i]);
    for (j=0; j<pdata->ncol; j++) 
      fprintf(stderr," %5.3f", pdata->d[i][j]);
    fprintf(stderr,"\n");
  }
}
/********************************************************************************/
/*                     write_outfile                                            */
/********************************************************************************/
/*Descriptions:
       write the test-statistics, unadjusted p-values, adjusted pvalues and Adjusted
       p-values lower to the file.
  input parameters:
      filename: the file to write
      pdata:    the pointer of the whole data
      T,P,Adj_P,Adj_Lower: the array stores the test-statistics,
      unadjusted p-values, adjusted pvalues and adjusted p-values lower, respectively
      if Adj_Lower==NULL, it will not print this item
*/
void write_outfile(FILE* fh,GENE_DATA* pdata,float*T, float*P,float*Adj_P,float* Adj_Lower)
{
  int i,nrow;
  /*float num,denum;*/
  nrow=pdata->nrow;  /*the length of the array T,P,etc.*/

  if(myDEBUG)
  {
    fprintf(stderr,"\nThe results of T,P Adj_P and Adj_Lower");
    print_farray(stderr,T,nrow);      
    print_farray(stderr,P,nrow);
    print_farray(stderr,Adj_P,nrow);
    if(Adj_Lower) print_farray(stderr,Adj_Lower,nrow);
  };
  fprintf(stderr,"\nWe're writing the output\n");
  fprintf(fh,"%20s %10s %10s %10s","gene_id","test-stat",
	  "unadj-p","adjusted-p");
  if(Adj_Lower)
    fprintf(fh,"%10s","p-lower");
  fprintf(fh,"\n");

  for (i=0; i<nrow; i++) {
    /*  t_stat_num_den(pdata->d[i],pdata->L, pdata->ncol,pdata->na,&num,&denum);*/
    fprintf(fh, "%20s %10.6f    %7g    %7g", 
	    pdata->id[i],T[i],P[i],Adj_P[i]);
    if(Adj_Lower){
      fprintf(fh,"    %7g",Adj_Lower[i]);
    }
    fprintf(fh,"\n");
  }
}

/*testing */
/*int main()
{
  #define N 5
  int V[N];
  int i,is_next=1;
  for(i=0;i<N;i++)
    V[i]=i;
  while(is_next){
    print_narray(stderr,V,N);
    is_next=next_permu(V,N);
  }
}
*/
int next_label_block(int* L, int n, int m)
{
  int s,l,i,j,block;/*s is for starting location*/
  int next=0;
  /*initialize the begining*/
  block=n/m;
  s=0;
  for(l=0;l<block;l++){
    next=next_permu(L+s,m);/*find the permutation within one block*/
    if(next) {
      /*resetting of the previous to the original*/
      for(i=0;i<l;i++){
	for(j=0;j<m;j++){
	  L[i*m+j]=j;
	}
      }
      return 1;
    }
    s+=m;
  }
  return 0;/*we couldn't find the next permutation*/
}
void init_label_block(int *L, int n,int m)
{
  int i,b,block;
  block=n/m;
  for(b=0;b<block;b++){
    for(i=0;i<m;i++){
      L[b*m+i]=i;
    }
  }
}
void sample_block(int *L, int n,int m)
{
  int block,b,s;
  s=0;
  block=n/m;
  for(b=0;b<block;b++){
    sample(L+s,m,m);
    s+=m;
  }
}
  
/*int main()
{ 
  GENE_DATA data;
  #define ROW 100
  float T[ROW];
  int L[38];
  int i,k=2;

  data.na=NA_FLOAT;
  data.nrow=ROW;  
  data.ncol=38;
  malloc_gene_data(&data);
  read_infile("data",&data); */
  /*print_gene_data(&data);
    compute_test_stat(&data,data.L,T,two_sample_tstat,NULL); checked
    compute_test_stat(&data,data.L,T,Fstat,(const void *)&k);
    compute_test_stat(&data,L,T,Block_Fstat,(const void *)&k);
  print_farray(stderr,T,ROW);
  free_gene_data(&data);
  }*/

