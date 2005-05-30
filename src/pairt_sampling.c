/*the l is for local global variable in this file*/
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"

static int l_n=0;/*the number of samples for permutations*/
static int l_B=0;/*the number of total simultaions*/
static int l_b=0;/* the number of permutations are done*/
static int l_is_random=1;/* the permuation is random or not*/
static unsigned int* l_all_samples=NULL;
/*store all the samples in random case*/
static int l_sz=0; /*the number of bytes for per permutation*/
static int l_len=0;
static int get_binpermu(int h,int n,int sz,int len,int *L,int hMax,unsigned int *V)/*sz=ceiling(n/sizeof(int)*8)*/;
static int set_binpermu(int *L,int h,int n,int sz,int len,int hMax,unsigned int *V);
void create_sampling_pairt(int n,int*L,int B)
{
  int i,maxB;
  unsigned int imax;
  l_n=n;
  l_b=0;
  imax=(unsigned int)(~0);
  l_len=floor(log(imax+1.0)/log(2));
  l_sz=ceil(n/(l_len*1.0));
  
  /*setting the maximum B*/
  if(fabs(n*log(2))<log(imax>>1)){ /*to be safe, moved two bits*/
    maxB=1<<n;    
  }else{/*we can set the only maximum of B*/
    maxB=imax>>1;
  }
  
  if((B==0) ||(B>=maxB)){
    if(n>=(l_len-1)){
      fprintf(stderr,"as n=%d is very large, we can not do complete permutation\n, Please try random permutation\n",n);
      return;
    }
    l_is_random=0;
    l_B=maxB;
    /*when exceeding the maximum numbers, we'll use the complete permutaions*/
    /*fprintf(stderr,"\nWe're doing %d complete permutations\n",l_B);*/
    Rprintf("\nWe're doing %d complete permutations\n",l_B);
  }
  else{
    int* myL;
    assert(myL=(int*)Calloc(n,int));
    l_B=B;
    l_is_random=1;
    /*fprintf(stderr,"\nWe're doing %d random permutations\n",l_B);*/
    Rprintf("\nWe're doing %d random permutations\n",l_B);
    set_seed(g_random_seed);

    assert(l_all_samples=(unsigned int*)Calloc(l_B*l_sz,int)); 
    /*setting the first sample as the original data*/
    set_binpermu(L,0,n,l_sz,l_len,l_B,l_all_samples);
    /*the extra as a buffer*/
    for(i=1;i<l_B;i++){
      int j;
      float tmp;
      for(j=0;j<n;j++){
	tmp=get_rand();
	if(tmp>0.5)
	  myL[j]=1;
	else
	  myL[j]=0; 
      }
      set_binpermu(myL,i,n,l_sz,l_len,l_B,l_all_samples);
    }
    Free(myL);
    if(myDEBUG)
    {
      fprintf(stderr,"the samples are\n");
      for(i=0;i<l_B;i++)
	fprintf(stderr,"%d ",l_all_samples[i]);
    }
  }
}

void delete_sampling_pairt()
{
  if(l_is_random){
    if(l_B!=0){
      Free(l_all_samples);
      l_all_samples=NULL;
    }
  }
}
int first_sample_pairt(int *L)
{
  /*return the number of real samples*/ 
  if(L==NULL)
    return l_B;
  
  /*call different sampling function*/  
  if(l_is_random){
    get_binpermu(0,l_n,l_sz,l_len,L,l_B,l_all_samples);
  }
  else
    int2bin(0,L,l_n);

  l_b=1;/*resetting the the number of permuatins done*/

  return 1;
}
int next_sample_pairt(int* L)
{ 
  if(l_b >=l_B)
    return 0; /*no next sample*/

  /*call different sampling function*/  
  if(l_is_random)
    get_binpermu(l_b,l_n,l_sz,l_len,L,l_B,l_all_samples);
  else
    int2bin(l_b,L,l_n); /* note for the complete resampling, we can not
			   do more than 2^32 times*/

  l_b++; 

  return 1; 
}
static int get_binpermu(int h,int n,int sz,int len,int *L,int hMax,unsigned int *V)/*sz=ceiling(n/sizeof(int)*8)*/
{
  int i,j;
  unsigned val;
  memset(L,0,sizeof(unsigned int)*n);
  if((h+1)> hMax) return 0;
  for(j=0;j<sz;j++){
    i=j*len; /*starting from the last bit*/
    val=V[h*sz+j];
    while(val>0){
      /*this code maybe faster if necessary*/
      L[i]=val&1;
      i++;
      val>>=1;/*to move another bit*/
    }
  }
  return 1;
}
      
      
static int set_binpermu(int *L,int h,int n,int sz,int len,int hMax,unsigned int *V)
{

  int i,j,nextbound;
  unsigned val,pow;
  if((h+1)> hMax) return 0;
  i=0; /*starting from the last bit*/
  for(j=0;j<sz;j++){
    nextbound=(j+1)*len;
    if(nextbound> n)
      nextbound=n;
    pow=1;
    val=0;
    while(i<nextbound){
      val+=(unsigned int)(L[i])*pow;
      pow<<=1;
      i++;
    }
    V[h*sz+j]=val;
  }
  return 1;
}
  
  





