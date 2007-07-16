/*the l is for local global variable in this file*/
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mt.h"

static int l_n=0;/*the number of samples for permutations*/
static int l_B=0;/*the number of total simultaions*/
static int l_b=0;/* the number of permutations are done*/
static int* l_L=NULL;
void create_sampling_pairt_fixed(int n,int*L,int B)
{
  l_n=n;
  l_B=B;
  l_b=0;
  if(B<=0){
    fprintf(stderr,"B needs to be positive\n");
    return;/*exit(0)*/;
  }
  assert(l_L=(int*)malloc(sizeof(int)*n));
  memcpy(l_L,L,sizeof(int)*n);
}


void delete_sampling_pairt_fixed()
{
  free(l_L);
  l_L=NULL;
}
int first_sample_pairt_fixed(int *L)
{
  if(L==NULL)
    return l_B;
  else{
    memcpy(L,l_L,sizeof(int)*l_n);
  }
  l_b=1;
  set_seed(g_random_seed);
  return 1;  
}
int next_sample_pairt_fixed(int* L)
{ 
  int n=l_n,i;
  float tmp;
  if(l_b>=l_B) return 0;
  for(i=0;i<n;i++){
    tmp=get_rand();
    if(tmp>0.5)
      L[i]=1;
    else
      L[i]=0;       
  l_b++;
  }
  return 1;
}

  
  





