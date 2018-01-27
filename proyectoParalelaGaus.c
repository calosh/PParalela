/**********************************************************************************************
* Matrix Multiplication Program using MPI.
*
* Viraj Brian Wijesuriya - University of Colombo School of Computing, Sri Lanka.
* 
* Works with any type of two matrixes [A], [B] which could be multiplied to produce a matrix [c].
*
* Master process initializes the multiplication operands, distributes the muliplication 
* operation to worker processes and reduces the worker results to construct the final output.
*  
************************************************************************************************/
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#include<mpi.h>

#define constante 1000 //rows of input [A]
#define NUM_ROWS_A 1000 //rows of input [A]
#define NUM_COLUMNS_A 1000 //columns of input [A]
#define NUM_ROWS_B 1000 //rows of input [B]
#define NUM_COLUMNS_B 1000 //columns of input [B]


#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 4 //tag for messages sent from slaves to master
void makeAB(); //makes the [A] and [B] matrixes
void makeFacyorial();
void printArray(); //print the content of output matrix [C];
void imprimir(double [][constante]); //print the content of output matrix [C];


int rank; //process rank
int size; //number of processes
int i, j, k; //helper variables

int ii, kk, l;
char opcion;


int fact;

// Matriz Inversa
void fila_pivote(void);
void col_pivote(void);
void otros(void);
void imprimir2(void);

void makeInversa(void);

double entrada=constante; // 

double pivote;

// Matriz Inicial
double A[NUM_ROWS_A][NUM_COLUMNS_A]; //declare input [A]
double B[NUM_ROWS_B][NUM_COLUMNS_B]; //declare input [B]

// Matriz Factorial
double mat_a[NUM_ROWS_A][NUM_COLUMNS_A]; //declare input [A]
double mat_b[NUM_ROWS_B][NUM_COLUMNS_B]; //declare input [B]
// Resualtado
double mat_result[NUM_ROWS_A][NUM_COLUMNS_B]; //declare output [C]

// Inversa
long double R[NUM_ROWS_A][NUM_COLUMNS_B]; //


double start_time; //hold start time
double end_time; // hold end time
int low_bound; //low bound of the number of rows of [A] allocated to a slave
int upper_bound; //upper bound of the number of rows of [A] allocated to a slave
int portion; //portion of the number of rows of [A] allocated to a slave
MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend
int main(int argc, char *argv[])
{   
    srand(time(NULL)); //El mayordomo pone a girar la diana

    MPI_Init(&argc, &argv); //initialize MPI operations
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
    /* master initializes work*/
    if (rank == 0) {
        printf("\n\n\t\t\t Matriz  A y B");
        makeAB();
        printf("\n\n\t\t\t Factorial de A y B");
        makeFacyorial();
        start_time = MPI_Wtime();

        printf("\n\n\t\t\t Multiplicacion");
        for (i = 1; i < size; i++) {//for each slave other than the master
            portion = (NUM_ROWS_A / (size - 1)); // calculate portion without master
            low_bound = (i - 1) * portion;
            if (((i + 1) == size) && ((NUM_ROWS_A % (size - 1)) != 0)) {//if rows of [A] cannot be equally divided among slaves
                upper_bound = NUM_ROWS_A; //last slave gets all the remaining rows
            } else {
                upper_bound = low_bound + portion; //rows of [A] are equally divisable among slaves
            }
            //send the low bound first without blocking, to the intended slave
            MPI_Isend(&low_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            //next send the upper bound without blocking, to the intended slave
            MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &request);
            //finally send the allocated row portion of [A] without blocking, to the intended slave
            MPI_Isend(&mat_a[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_A, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
        }
    }
    //broadcast [B] to all the slaves
    MPI_Bcast(&mat_b, NUM_ROWS_B*NUM_COLUMNS_B, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* work done by slaves*/
    if (rank > 0) {
        //receive low bound from the master
        MPI_Recv(&low_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        //next receive upper bound from the master
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        //finally receive row portion of [A] to be processed from the master
        MPI_Recv(&mat_a[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_A, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        for (i = low_bound; i < upper_bound; i++) {//iterate through a given set of rows of [A]
            for (j = 0; j < NUM_COLUMNS_B; j++) {//iterate through columns of [B]
                for (k = 0; k < NUM_ROWS_B; k++) {//iterate through rows of [B]
                    mat_result[i][j] += (mat_a[i][k] * mat_b[k][j]);
                }
            }
        }
        //send back the low bound first without blocking, to the master
        MPI_Isend(&low_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        //send the upper bound next without blocking, to the master
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        //finally send the processed portion of data without blocking, to the master
        MPI_Isend(&mat_result[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_B, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
    }
    /* master gathers processed work*/
    if (rank == 0) {
        for (i = 1; i < size; i++) {// untill all slaves have handed back the processed data
            //receive low bound from a slave
            MPI_Recv(&low_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            //receive upper bound from a slave
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            //receive processed data from a slave
            MPI_Recv(&mat_result[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_B, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        }
        
        //printArray();
        printf("\n\n\t\t\t Inversa");

        // Inversa
        makeInversa();
        //imprimir(mat_result);

        end_time = MPI_Wtime();
        printf("\nRunning Time = %f\n\n", end_time - start_time);

    }
    MPI_Finalize(); //finalize MPI operations
    return 0;
}

void makeAB()
{
    for (i = 0; i < NUM_ROWS_A; i++) {
        for (j = 0; j < NUM_COLUMNS_A; j++) {
            A[i][j] = 1+rand()%(5-1);
        }
    }
    for (i = 0; i < NUM_ROWS_B; i++) {
        for (j = 0; j < NUM_COLUMNS_B; j++) {
            B[i][j] = 1+rand()%(5-1);
        }
    }
}



void makeFacyorial(){
    // La matriz AF y BF se llenan con los factoriales de las Matrices A y B

    k = NUM_ROWS_A;
    for (i = 0;i < k; i++)
    {
      for (j = 0;j < k; j++)
         {
          fact = 1;
          for (int b = A[i][j]; b > 1; b--){
            fact = fact * b;
          }

          mat_a[i][j] = fact;

          fact = 1;
          for (int b = B[i][j]; b > 1; b--){
            fact = fact * b;
          }

          mat_b[i][j] = fact;
          }
      }
}

void makeInversa(){
  // Inversa
  //printf("\n\n\t\t\t Inversa");
  for(ii=0; ii<entrada; ii++)
  {
      j=ii;
      pivote=mat_result[ii][j];
      R[ii][j]=1/pivote;
      fila_pivote();
      col_pivote();
      otros();
      for(kk=0; kk<entrada; kk++)
          for(l=0; l<entrada; l++)
              mat_result[kk][l]=R[kk][l];
  }
}

void printArray()
{
    printf("\n\n\t\t\t Factorial de A");
    for (i = 0; i < NUM_ROWS_A; i++) {
        printf("\n");
        for (j = 0; j < NUM_COLUMNS_A; j++)
            printf("%8.2f  ", mat_a[i][j]);
    }
    printf("\n\n\n");
    printf("\n\n\t\t\t Factorial de B");
    for (i = 0; i < NUM_ROWS_B; i++) {
        printf("\n");
        for (j = 0; j < NUM_COLUMNS_B; j++)
            printf("%8.2f  ", mat_b[i][j]);
    }
    printf("\n\n\n");
    printf("\n\n\t\t\t Multiplicacion de Fac A*B");
    for (i = 0; i < NUM_ROWS_A; i++) {
        printf("\n");
        for (j = 0; j < NUM_COLUMNS_B; j++)
            printf("%8.2f  ", mat_result[i][j]);
    }
    printf("\n\n");
}



void fila_pivote(void)
{
    int m;
    for(m=0; m<entrada; m++)
        if(m != ii)
            R[ii][m]=mat_result[ii][m]/pivote;
}

void col_pivote()
{
    int m;
    for(m=0; m<entrada; m++)
        if(m != j)
            R[m][j]=-mat_result[m][j]/pivote;
}

void otros(void)
{
    int x,y;
    for(x=0 ;x<entrada; x++)
        for(y=0; y<entrada; y++)
            if(x!=ii && y!=j)
                R[x][y]=mat_result[x][y]-(mat_result[ii][y]*mat_result[x][j])/pivote;
}


void imprimir(double matrix[constante][constante]){
  printf("\n");
  for(int i=0;i<entrada;i++){
      printf("\n\t\t");
      for(int j=0;j<entrada;j++){
          printf("  %f  ", matrix[i][j]);
      }
  }
  printf("\n");
}