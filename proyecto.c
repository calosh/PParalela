#include<stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>

void makeAB(); //makes the [A] and [B] matrixes

void imprimir(float [][25]); //print the content of output matrix [C];


float determinant(float [][25], float);
void cofactor(float [][25], float);
void transpose(float [][25], float [][25], float);

float A[25][25];
float B[25][25];  
float AF[25][25];
float BF[25][25];
float AFBF[25][25], k, d; // Matriz Inversa
int i, j;


int main()
{
  srand(time(NULL)); //El mayordomo pone a girar la diana

  printf("Enter the order of the Matrix : ");
  scanf("%f", &k);
  printf("Enter the elements of %.0fX%.0f Matrix : \n", k, k);
  
  // Se llena aleatoriamente la matriz A y B
  makeAB();
  printf("\n\n\t\t\t Matriz de A");
  imprimir(A);
  printf("\n\n\t\t\t Matriz de B");
  imprimir(B);

  // La matriz AF y BF se llenan con los factoriales de las Matrices A y B
  for (i = 0;i < k; i++)
  {
    for (j = 0;j < k; j++)
       {
        //scanf("%f", & I[i][j]);
         //I[i][j] = 1+rand()%(10-1);
        // FACTORIAL
        int fact = 1;
        for (int b = A[i][j]; b > 1; b--){
          fact = fact * b;
        }

        AF[i][j] = fact;

        fact = 1;
        for (int b = B[i][j]; b > 1; b--){
          fact = fact * b;
        }

        BF[i][j] = fact;
        }
    }

    printf("\n\n\t\t\t Factoriales de A");
    imprimir(AF);
    printf("\n\n\t\t\t Factoriales de B");
    imprimir(BF);


    // Multiplicacion de la Matriz AF y BF
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            AFBF[i][j]=0;
            for(int l=0;l<k;l++){
                AFBF[i][j]=(AFBF[i][j]+(AF[i][l]*BF[l][j]));
            }
        }
    }

    printf("\n\n\t\t\t Multiplicacion de AF * BF");
    imprimir(AFBF);

  
  // Se saca el determinante de la matriz
  d = determinant( AFBF, k);
  if (d == 0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor( AFBF, k);
}
 
/*For calculating Determinant of the Matrix */
float determinant(float  AFBF[25][25], float k)
{
  float s = 1, det = 0, b[25][25];
  int i, j, m, n, c;
  if (k == 1)
    {
     return ( AFBF[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] =  AFBF[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * ( AFBF[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}
 
void cofactor(float num[25][25], float f)
{
 float b[25][25], fac[25][25];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f);
}
/*Finding transpose of matrix*/ 
void transpose(float num[25][25], float fac[25][25], float r)
{
  int i, j;
  float b[25][25], inverse[25][25], d;
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is : \n");
 
   for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         printf("\t%f", inverse[i][j]);
        }
    printf("\n");
     }
}


// https://stackoverflow.com/questions/8671366/undefined-reference-to-pow-and-floor


void makeAB()
{
  for (i = 0;i < k; i++)
    {
     for (j = 0;j < k; j++)
       {

        A[i][j] = 1+rand()%(5-1);
        B[i][j] = 1+rand()%(5-1);

        AFBF[i][j] = 1+rand()%(5-1);

        }
    }
}


void imprimir(float matrix[25][25]){
  printf("\n");
  for(int i=0;i<k;i++){
      printf("\n\t\t");
      for(int j=0;j<k;j++){
          printf("  %6f  ", matrix[i][j]);
      }
  }
  printf("\n");
}
