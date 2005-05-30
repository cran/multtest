/* mt/random.c
   Copyright (C)  Yongchao Ge, Berkeley, USA
   This program is part of the the free software of the multtest, 
   which is a C stand alone package, the multtest software has also 
   been wrapped to R package for people easy to use.
 */
  

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mt.h"
/*-------------------------------ran2()---------------------------*/
/* This function is implementation of random number generating algorithm 
   based on ran2 function in the book of Numerical Recipes in C pp282(1992), 2nd edition

   It generates a unifrom random number from 0 to 1 (not including 0 and 1) using 
   random number generator of L'Ecuyer with Bayes_Durham Shuffle and added safeguards.
*/
#define A1 40014
#define Q1 53668
#define R1 12211
#define M1 2147483563
#define A2 40692
#define Q2 52774
#define R2 3791
#define M2 2147483399 

#define ONE (1.0- 1.2e-7)
#define N_SHUFF 32
#define N_WARMUP 8

/*prestoring the results to increase the computing speed*/
static long int  N_DIV=(1+(M1-1)/N_SHUFF); 
static float  M1inv=(1.0/M1);

/*---------------------------begin of RNG------------------------------------*/
typedef struct tagRNG{
  long int z1; /*the generatore 1*/
  long int z2; /*the generator 2*/
  long int y;
  long int V[N_SHUFF];
}RNG;
RNG l_rng;
void set_seed(long int seed){
  long int z1,z2,*V;
  int i;
  z1=abs(seed); /*setting the seed to z1 and initialize it*/
  if(z1==0) z1=1;/*be sure to prevent seed=0*/
  z2=z1;/*initializing z2*/
  /*warm up*/
  for(i=0;i<N_WARMUP;i++){
     long int t=z1/Q1;
     z1=A1*(z1-t*Q1)-R1*t;
     if(z1<0) z1+=M1;
  }
  /*initializing the array V*/
  V=l_rng.V;
  for(i=N_SHUFF-1;i>=0;i--){
    long int t=z1/Q1;
     z1=A1*(z1-t*Q1)-R1*t;
     if(z1<0) z1+=M1;
     V[i]=z1;
  }
  l_rng.z1=z1;
  l_rng.z2=z2;
  l_rng.y=z1;
}
float get_rand(){
  int i;
  long t;
  long int z1=l_rng.z1,z2=l_rng.z2,y=l_rng.y,*V=l_rng.V; 
  float res;
  /*generator 1*/
  t=z1/Q1;
  z1=A1*(z1-t*Q1)-R1*t;
  if(z1<0) z1+=M1;
  /*generator 2*/
  t=z2/Q2;
  z2=A2*(z2-t*Q2)-R2*t;
  if(z2<0) z2+=M2;
  /*shuffling*/
  i=y/N_DIV; /*N_DIV=(1+P1/N_SHUFF);to make sure i to be 0..N_SHUFF-1*/
  y=V[i]-z2; /*V[i] is a random number similar to z1, y will be in  {-(m2-1)} .. (m1-1)*/
  if(y<1) y+=(M1-1); /*to make sure y will be 1...(m1-1) note m1 is almost equal to m2,m1>m2*/
  V[i]=z1; /*filling the new generated random number z1 to array V*/ 
  /*normalizing from 0 to 1*/
  res=y*M1inv; /*=y/M1*/
  l_rng.z1=z1;
  l_rng.z2=z2;
  l_rng.y=y;
  if(res>ONE) return ONE;
  else return res;
}
/*-----------------------end of RNG----------------------------*/
/*get the n samples from the n-dim vector V. the results are stored 
in the first m member of vector V*/
void sample(int *V, int n, int m)
{
  int i,j,temp;
  float f;
  for(i=0;i<m;i++){
    /* no need to worry yet    
       if(i==(n-1)) continue; */ /*no need to swap with the last elements*/
    j=n;
    while (j==n){/*skip the border, we only want random
		    numbers from i,i+1,i+2,...,n-1*/
      f=get_rand()*(n-i);
      j=i+floor(f);
    }
    /*swap the nubmer V[i] and V[j] whther if i==j*/
    temp=V[j];
    V[j]=V[i];
    V[i]=temp;
  }   
}
  

/*void main(int argc, char* argv[])
{
  int n;
  long seed=100;
  int i;
  float temp;
 #define DIM 20
  int V[DIM];
  n=atoi(argv[1]);
  seed=atol(argv[2]);
  set_seed(seed);
  for(i=0;i<n;i++){
    temp=get_rand();
    fprintf(stderr,"%f ",temp);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"The sampling of %d\n",DIM);
  for(i=0;i<DIM;i++)
    V[i]=i+1;
  sample(V,DIM,DIM);
  for(i=0;i<DIM;i++)
  fprintf(stderr,"%d ",V[i]);
}*/

  
  
