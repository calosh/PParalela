#include <stdio.h>
#include <stdlib.h>

float pivote, a[50][50], b[50][50];
int i, j, n, k, l;
char opcion;

void lee_matriz(void);
void fila_pivote(void);
void col_pivote(void);
void otros(void);
void imprimir(void);

int main(void)
{
    lee_matriz();
    for(i=0; i<n; i++)
    {
        j=i;
        pivote=a[i][j];
        b[i][j]=1/pivote;
        fila_pivote();
        col_pivote();
        otros();
        for(k=0; k<n; k++)
            for(l=0; l<n; l++)
                a[k][l]=b[k][l];
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
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
        {
            printf("Elemento A[%d][%d]: ", i+1, j+1);
            scanf("%f", &a[i][j]);
        }
}

void fila_pivote(void)
{
    int m;
    for(m=0; m<n; m++)
        if(m != i)
            b[i][m]=a[i][m]/pivote;
}

void col_pivote()
{
    int m;
    for(m=0; m<n; m++)
        if(m != j)
            b[m][j]=-a[m][j]/pivote;
}

void otros(void)
{
    int x,y;
    for(x=0 ;x<n; x++)
        for(y=0; y<n; y++)
            if(x!=i && y!=j)
                b[x][y]=a[x][y]-(a[i][y]*a[x][j])/pivote;
}

void imprimir(void)
{
    puts("\nMatriz inversa:\n");
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
            printf("%4.2f ", a[i][j]);
        printf("\n");
    }
}
