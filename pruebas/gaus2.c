#include <stdio.h>
#include <stdlib.h>

float pivote, AFBF[50][50], R[50][50];
int ii, j, n, k, l;
char opcion;

void lee_matriz(void);
void fila_pivote(void);
void col_pivote(void);
void otros(void);
void imprimir(void);

int main(void)
{
    lee_matriz();
    for(ii=0; ii<n; ii++)
    {
        j=ii;
        pivote=AFBF[ii][j];
        R[ii][j]=1/pivote;
        fila_pivote();
        col_pivote();
        otros();
        for(k=0; k<n; k++)
            for(l=0; l<n; l++)
                AFBF[k][l]=R[k][l];
    }
    imprimir();
    puts("\n");
    system("PAUSE");
    return EXIT_SUCCESS;
}

void lee_matriz(void)
{
    printf("INVERSA DE UNA MATRIZ nxn\n\nn: ");
    scanf("%d",&n);
    putchar('\n');
    for(ii=0; ii<n; ii++)
        for(j=0; j<n; j++)
        {
            printf("Elemento A[%d][%d]: ", ii+1, j+1);
            scanf("%f", &AFBF[ii][j]);
        }
}

void fila_pivote(void)
{
    int m;
    for(m=0; m<n; m++)
        if(m != ii)
            R[ii][m]=AFBF[ii][m]/pivote;
}

void col_pivote()
{
    int m;
    for(m=0; m<n; m++)
        if(m != j)
            R[m][j]=-AFBF[m][j]/pivote;
}

void otros(void)
{
    int x,y;
    for(x=0 ;x<n; x++)
        for(y=0; y<n; y++)
            if(x!=ii && y!=j)
                R[x][y]=AFBF[x][y]-(AFBF[ii][y]*AFBF[x][j])/pivote;
}

void imprimir(void)
{
    puts("\nMatriz inversa:\n");
    for(ii=0; ii<n; ii++)
    {
        for(j=0; j<n; j++)
            printf("%4.2f ", AFBF[ii][j]);
        printf("\n");
    }
}
