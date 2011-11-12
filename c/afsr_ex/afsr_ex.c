#include <stdio.h>
#include <limits.h>
#include <math.h>

int main(){
    int i,j,k; //loop variables
    int p=5,r=3,n=140; //define prime, length of register, m, and ouput length
    int *m=calloc(n,sizeof(int)); m[0]=4;
//    int S[]={0,1,2,3,4}; //set S for a's and T for tap q's
//    int T[]={0,1,2,3,4};
    int *a=calloc(n,sizeof(int));   //active register of a's
    int *q=calloc(r+1,sizeof(int)); //permanent q taps
    int tau=0;

    a[0]=1; a[1]=2; a[2]=2; //initial active register
    q[0]=2; q[1]=2; q[2]=1; q[3]=2; //permanent q taps
    int qinv=3; //specific to mod 5

    for (i=r; i<n; i++){
       tau=0;
       for (j=1; j<r+1; j++){
	  tau+=q[j]*a[i-j];
       }
       tau+=m[i-r];
       a[i]=(qinv*tau)%p;
       k=1; m[i-r+1]=0;
       while (p*k<=tau){
	  m[i-r+1]+=1;
	  k++;
       }
    }

    for (i=0; i<n; i++){
       printf("%d",a[i]);
    }
    printf("\n");
    for (i=0; i<n; i++){
       printf("%d",m[i]);
    }
    printf("\n");

    return 0;
}
