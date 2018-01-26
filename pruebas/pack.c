#include<stdio.h>
#include<math.h>

#include<complex.h>
#include <stdlib.h>
void zgeTranspose( double complex *Transposed, double complex *M ,int n);
void MatrixComplexInverse(double complex *invA, double complex *A, int n);


//........................................................................................
void zgeTranspose( double complex *Transposed, double complex *M ,int n)
{

int i,j;
for(i=0;i<n;i++)

for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

//.........................................................................................
void MatrixComplexInverse(double complex *invA, double complex *A, int n)
{

int LWORK=10*n;

int *permutations;

double complex *WORK, *tempA;

tempA = (double complex*) malloc( n*n*sizeof(double complex) );

permutations = (int*) malloc( 2*n*sizeof(int) );

WORK = (double complex *)malloc(LWORK*sizeof(double complex));


int INFO;

zgeTranspose(tempA,A,n);


zgetrf_( &n, &n, tempA , &n, permutations , &INFO );

if (INFO != 0) {
 printf("ComplexMatrixInverse: Error at zgetrf  \n"); exit(0);
 }



zgetri_( &n, tempA , &n, permutations , WORK, &LWORK, &INFO );

if (INFO != 0) {
 printf("ComplexMatrixInverse: Error at zgetri  \n"); exit(0);
 }

zgeTranspose(invA,tempA,n);

free(WORK);

free(tempA);
free(permutations);

}


/////////////////////////////////////////////////////////////////////////////////////////


int main()
{
int i,j;
const int N = 3;


double complex A[] = { 1.+I , 2. ,  3 , 4. , 5.+I , 6. , 7., 8., 9. + I};

double complex invA[N*N];

MatrixComplexInverse(invA,A,N);


for(i=0;i<N;i++){
 for(j=0;j<N;j++) printf(" (%f,%f) \t", invA[i*N + j]);

 printf("\n");
}

printf("------------------------------------------------------------\n");
return 0;

}

