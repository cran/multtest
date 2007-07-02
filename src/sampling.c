/*the l is for local global variable in this file*/
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
typedef struct tagPERMU_ARRAY{
  int n; /*the number of original observations (samples) needs to permute*/
  int k; /* the number of classes, labelled from 0..(k-1)
	  which functions as the base of the integar representation*/
  int* nk; /* the number of groups in class 0..(k-1)*/
  int B;/* the number of permutations samples*/
  int len; /*len= floor(log(imax,k)), where imax the maximum of integers */
  int sz;/*the number of integars for each permutation needed sz=ceil(n/len))*/
  unsigned int * v;/* the array, which has size of B*sz)
 		      unsigned integars*/
}PERMU_ARRAY;
static int init_permu_array(PERMU_ARRAY* pa,int *L,int n, int B);
static int get_permu(PERMU_ARRAY* pa, int h,int *L);
/* get the h-th permutation of the permu_array. L needs to be array
   of length pa->n*/
static int set_permu(PERMU_ARRAY* pa, int h,int *L);
static void delete_permu_array(PERMU_ARRAY* pa);


static int l_b=0; /* the number of permutations are done*/
static int l_B=0; /*the number of all permutations */

static PERMU_ARRAY l_pa;

/*store all the samples in random case, the first one needs to be from the original data*/
void create_sampling(int n,int*L,int B)
{  
  int i,rest,maxB=0;
  int imax;
  double f;
  /*initiate the prelim computation*/
  init_permu_array(&l_pa,L,n,0);
  
  /*setting the value of f=log(maxB)*/
  f=0;
  rest=n;
  for(i=0;i<l_pa.k;i++){
    f+=logbincoeff(rest,l_pa.nk[i]);
    rest-=l_pa.nk[i];
  }

  /*setting the maximum B*/
  imax=(unsigned int)(~0)>>1;/*divide by 2 to avoid the negative number*/
  if(fabs(f)<log(imax)){
    maxB=1;
    rest=n;
    for(i=0;i<l_pa.k;i++){
      maxB*=bincoeff(rest,l_pa.nk[i]);
      rest-=l_pa.nk[i];
    }
  }else{/*we can set the only maximum of B*/
    maxB=imax;
  }

  /*to check random or complete*/
  if((B<=0) || (B>=maxB)){
    /* checking if complete permutation doable*/
    if (fabs(f)>log(imax)){
      fprintf(stderr,"as B(log(B)=%5.2lf) is too big,we can not do the complete permutations\n",f);
      return;/*exit(0);*/
    }
    /*when exceeding the maximum numbers, we'll use the complete permutaions*/
    l_B=maxB;
/*    fprintf(stderr,"\nWe're doing %d complete permutations\n",l_B);*/
    Rprintf("\nWe're doing %d complete permutations\n",l_B);
  }else{
    /*doing random permutation*/
    int * ordern,* permun,*myL;
    l_B=B;
    /*fprintf(stderr,"\nWe're doing %d random permutations\n",l_B);*/
    Rprintf("\nWe're doing %d random permutations\n",l_B);
    /*reintiailize the permu_array*/
    delete_permu_array(&l_pa);
    init_permu_array(&l_pa,L,n,B);
    assert(permun=(int*)Calloc(l_pa.n,int));
    assert(ordern=(int*)Calloc(l_pa.n,int));
    assert(myL=(int*)Calloc(l_pa.n,int));
    for(i=0;i<n;i++){
      ordern[i]=i;
    }
    /*allocate and assign the values for l_first_sample*/
    set_permu(&l_pa,0,L);
    set_seed(g_random_seed);
    for(i=1;i<B;i++){
      memcpy(permun,ordern,sizeof(int)*n);
      sample(permun,n,n);
      /*change to labbeling*/
      sample2label(n,l_pa.k,l_pa.nk,permun,myL);
      set_permu(&l_pa,i,myL);
    }
    Free(myL);
    Free(ordern);
    Free(permun);
  }
}
void delete_sampling()
{
  delete_permu_array(&l_pa);
}

int first_sample(int *L)
{
  if(L==NULL)
    return l_B;

  /*if is random, we'll choose the original L as the first 
    sample*/
  if(l_pa.B > 0){
    get_permu(&l_pa,0,L);
  }else{ 
    init_label(l_pa.n,l_pa.k,l_pa.nk,L);
  }
  l_b=1;/*resetting the the number of permuatins done*/
  /*print_narray(L,16);*/

  return 1;
}

int next_sample(int* L)
{
  if(l_b>=l_B) return 0;

  if(l_pa.B > 0){
    get_permu(&l_pa,l_b,L);
  }    
  else{
    next_label(l_pa.n,l_pa.k,l_pa.nk,L);
  }
  l_b++;
  return 1;
}

static int init_permu_array(PERMU_ARRAY* pa, int *L,int n, int B)
{
  int i;
  unsigned imax;
  pa->n=n;
  pa->B=B;
  pa->nk=NULL;
  pa->v=NULL;

  /* compute the k*/
  pa->k=0;
  for(i=0;i<n;i++)
    if(L[i]>pa->k)
      pa->k=L[i];
  (pa->k)++;
  
  /*compue nk*/
  assert(pa->nk=(int*)Calloc(pa->k,int));
  memset(pa->nk,0,sizeof(int)*pa->k);
  for(i=0;i<n;i++)
    pa->nk[L[i]]++;
  
  /*computer imax, len*/
  imax=~0; /*get all bits are 1 for the integars*/
  pa->len=floor(log(imax+1.0)/log(pa->k)); 
  pa->sz=ceil(n/(pa->len*1.0));
  /*allocate the space for v*/
  assert(pa->v=(unsigned int*)Calloc(B*pa->sz,int));
  return 1;
}



static int get_permu(PERMU_ARRAY* pa, int h, int *L)
{
  int i,j;
  unsigned val;
  memset(L,0,sizeof(unsigned int)*pa->n);
  if((h+1)> pa->B) return 0;
  for(j=0;j<pa->sz;j++){
    i=j*pa->len; /*starting from the last bit*/
    val=pa->v[h*pa->sz+j];
    while(val>0){
      /*this code maybe faster if necessary*/
      L[i]=val%(unsigned int)(pa->k);
      i++;
      val/=(unsigned int)(pa->k);/*to move another bit*/
    }
  }
  return 1;
}
      
      
static int set_permu(PERMU_ARRAY* pa, int h,int *L)
{
  int i,j,nextbound;
  unsigned val,pow;
  if((h+1)> pa->B) return 0;
  i=0; /*starting from the last bit*/
  for(j=0;j<pa->sz;j++){
    nextbound=(j+1)*pa->len;
    if(nextbound> (pa->n))
      nextbound=pa->n;
    pow=1;
    val=0;
    while(i<nextbound){
      val+=(unsigned int)(L[i])*pow;
      pow*=(unsigned int)pa->k;
      i++;
    }
    pa->v[h*pa->sz+j]=val;
  }
  return 1;
}
static void delete_permu_array(PERMU_ARRAY* pa)
{
  Free(pa->nk);
  pa->nk=NULL;
  if(pa->B!=0){
    Free(pa->v);
    pa->v=NULL;
  }
}


  
  
