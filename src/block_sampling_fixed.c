/*This file is used to do the sampling for block sampling*/
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"

static int l_n=0;
static int l_B=0;/*the number of total simultaions*/
static int l_b=0;/* the number of permutations are done*/
static int l_is_random=1;/* the permuation is random or not*/
static int* l_L=NULL;
static int l_m=0;/*the number of treaments*/
static int* l_order_block=NULL;
void create_sampling_block(int n,int*L,int B)
{
  int i,maxB,Nblock,m,imax,fac;/*m is the number of treatments*/
  double logfac;
  m=0;
  for(i=0;i<n;i++)
    if(L[i]>m){
      m++;
    }
  m++;
  Nblock=n/m;
  logfac=logfactorial(m,m)*Nblock;
  imax=(unsigned int)(~0)>>1;/*divide by 2 to avoid the negative number*/
  if(fabs(logfac)<log(imax)){
    fac=1;
    for(i=1;i<m+1;i++)
      fac*=i;
    maxB=fac;
    for(i=1;i<Nblock;i++)
      maxB*=fac;
  }else{
    maxB=imax;
  }
  if((B<=0) || (B>=maxB)){
    /* checking if complete permutation doable*/
    if (fabs(logfac)>log(imax)){
      fprintf(stderr,"as B(log(B)=%5.2f) is too big,we can not do the complete permutations\n",logfac);
      return; /*exit(0)*/
    }
    l_B=maxB;
    fprintf(stderr,"\nWe're doing %d complete permutations\n",l_B);
    l_is_random=0;
  }else{
    /*doing random permutation*/
    l_B=B;
    l_is_random=1;
    set_seed(g_random_seed);
  }
    l_n=n;
    l_b=0;
    l_m=m;
    assert(l_L=(int*)Calloc(n,int));
    memcpy(l_L,L,sizeof(int)*n);
    assert(l_order_block=(int*)Calloc(n,int));
    init_label_block(l_order_block,n,m);
}
      
  
void delete_sampling_block()
{
  free(l_L);
  l_L=NULL;
  free(l_order_block);
}
int next_sample_block(int* L)
{
  if(l_b>=l_B) return 0;

  if(l_is_random){
    memcpy(L,l_order_block,sizeof(int)*l_n);
    sample_block(L,l_n,l_m);
  } else{
    next_label_block(L,l_n,l_m);
  }
  l_b++;
  return 1;
}
int first_sample_block(int *L)
{
  if(L==NULL)
    return l_B;
  if(l_is_random){
    memcpy(L,l_L,sizeof(int)*l_n);
  }else{ 
    init_label_block(L,l_n,l_m);
  }
  l_b=1;
  return 1;
}





