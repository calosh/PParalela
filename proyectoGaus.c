#include<stdio.h>
#include<math.h>
#include <time.h>
#include <stdlib.h>



void makeAB(); //makes the [A] and [B] matrixes

void imprimir(long double [][1000]); //print the content of output matrix [C];


void fila_pivote(void);
void col_pivote(void);
void otros(void);
void imprimir2(void);


double pivote;


long double A[1000][1000];
long double B[1000][1000];  
long double AF[1000][1000];
long double BF[1000][1000];
long double AFBF[1000][1000]; 
double entrada; // 

long double R[1000][1000]; //


int i, j;


int ii, k, l;
char opcion;

int main()
{
  srand(time(NULL)); //El mayordomo pone a girar la diana

  printf("Enter the order of the Matrix : ");
  scanf("%lf", &entrada);
  //printf("Enter the elements of %.0fX%.0f Matrix : \n", entrada, entrada);
  
  // Se llena aleatoriamente la matriz A y B
  makeAB();
  printf("\n\n\t\t\t Matriz de A");
  //imprimir(A);
  printf("\n\n\t\t\t Matriz de B");
  //imprimir(B);

  // La matriz AF y BF se llenan con los factoriales de las Matrices A y B
  for (i = 0;i < entrada; i++)
  {
    for (j = 0;j < entrada; j++)
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
    //imprimir(AF);
    printf("\n\n\t\t\t Factoriales de B");
    //imprimir(BF);


    // Multiplicacion de la Matriz AF y BF
    for(int i=0;i<entrada;i++){
        for(int j=0;j<entrada;j++){
            AFBF[i][j]=0;
            for(int l=0;l<entrada;l++){
                AFBF[i][j]=(AFBF[i][j]+(AF[i][l]*BF[l][j]));
            }
        }
    }

    printf("\n\n\t\t\t Multiplicacion de AF * BF");
    //imprimir(AFBF);

    // Inversa
    //printf("\n\n\t\t\t Inversa");
    for(ii=0; ii<entrada; ii++)
    {
        j=ii;
        pivote=AFBF[ii][j];
        R[ii][j]=1/pivote;
        fila_pivote();
        col_pivote();
        otros();
        for(k=0; k<entrada; k++)
            for(l=0; l<entrada; l++)
                AFBF[k][l]=R[k][l];
    }
    //imprimir2();
    printf("\n\n\t\t\t Inversa de AF * BF");
    //imprimir(AFBF);

}
 
// https://stackoverflow.com/questions/8671366/undefined-reference-to-pow-and-floor


void makeAB()
{
  for (i = 0;i < entrada; i++)
    {
     for (j = 0;j < entrada; j++)
       {

        A[i][j] = 1+rand()%(15-1);
        B[i][j] = 1+rand()%(15-1);
        //AFBF[i][j] = 1+rand()%(10-1);

        }
    }
}


void imprimir(long double matrix[1000][1000]){
  printf("\n");
  for(int i=0;i<entrada;i++){
      printf("\n\t\t");
      for(int j=0;j<entrada;j++){
          printf("  %.20Lf  ", matrix[i][j]);
      }
  }
  printf("\n");
}





// Gaus


void fila_pivote(void)
{
    int m;
    for(m=0; m<entrada; m++)
        if(m != ii)
            R[ii][m]=AFBF[ii][m]/pivote;
}

void col_pivote()
{
    int m;
    for(m=0; m<entrada; m++)
        if(m != j)
            R[m][j]=-AFBF[m][j]/pivote;
}

void otros(void)
{
    int x,y;
    for(x=0 ;x<entrada; x++)
        for(y=0; y<entrada; y++)
            if(x!=ii && y!=j)
                R[x][y]=AFBF[x][y]-(AFBF[ii][y]*AFBF[x][j])/pivote;
}

void imprimir2(void)
{
    puts("\nMatriz inversa:\n");
    for(ii=0; ii<entrada; ii++)
    {
        for(j=0; j<entrada; j++)
            printf("%Lf ", AFBF[ii][j]);
        printf("\n");
    }
}
