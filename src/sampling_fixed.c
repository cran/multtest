/*the l is for local global variable in this file*/
#include "assert.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "mt.h"
static int l_n=0; /*the length of L*/
static int l_k=0; /*the number of groups*/
static int* l_nk=NULL;/* the number of objects in each groups*/
static int* l_L=NULL;/*the storrage of first label*/
static int l_b=0; /* the number of permutations are done*/
static int l_B=0; /*the number of all permutations */
static int* l_permun=NULL;
static int* l_ordern=NULL;

void create_sampling_fixed(int n,int*L,int B)
{  
  int i,k;
  l_n=n;
  l_B=B;
  l_b=0;
  if(B<=0){
    fprintf(stderr,"B needs to be positive\n");
    return;/*exit(0)*/;
  }
  assert(l_L=(int*)Calloc(n,int));
  memcpy(l_L,L,sizeof(int)*n);
  
  k=0;
  for(i=0;i<n;i++)
    if(L[i]>k)
      k=L[i];
  k++;
  l_k=k;
  assert(l_nk=(int*)Calloc(k,int));
  memset(l_nk,0,sizeof(int)*k);
  for(i=0;i<n;i++)
    l_nk[L[i]]++;

  assert(l_permun=(int*)Calloc(n,int));
  assert(l_ordern=(int*)Calloc(n,int));
  for(i=0;i<n;i++){
      l_ordern[i]=i;
    }
}
void delete_sampling_fixed()
{
  Free(l_L);
  l_L=NULL;
  Free(l_nk);
  l_nk=NULL;
  Free(l_permun);
  l_permun=NULL;
  Free(l_ordern);
  l_ordern=NULL;
}

int first_sample_fixed(int *L)
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

int next_sample_fixed(int* L)
{
  int n=l_n;
  if(l_b>=l_B) return 0;
  memcpy(l_permun,l_ordern,sizeof(int)*n);
  sample(l_permun,n,n);
  /*change to labbeling*/
  sample2label(n,l_k,l_nk,l_permun,L);
  l_b++;
  return 1;
}

  
  
