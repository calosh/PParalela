#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 16 /*assumption: SIZE a multiple of number of nodes*/
#define FROM_MASTER 1/*setting a message type*/
#define FROM_WORKER 2/*setting a message type*/
#define DEBUG 1/*1 = debug on, 0 = debug off*/

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];
static double b_to_trans[SIZE][SIZE];
static void init_matrix(void)
{
    int i, j;
    for (i = 0; i < SIZE; i++)
    {
        for (j = 0; j < SIZE; j++) {
            a[i][j] = 1.0;
            if(i >= SIZE/2) a[i][j] = 2.0;
            b_to_trans[i][j] = 1.0;
            if(j >= SIZE/2) b[i][j] = 2.0;
//          c[i][j] = 1.0;
        }
    }
}

static void print_matrix(void)
{
    int i, j;
    for(i = 0; i < SIZE; i++) {
        for(j = 0; j < SIZE; j++) {
            printf("%7.2f", c[i][j]);
        }
    printf("\n");
    }
}

static void transpose_matrix()
{
    int i, j;
    for(i = 0; i<SIZE; i++)
        for(j = 0; j<SIZE;j++)
            b[i][j] = b_to_trans[j][i];
}

int main(int argc, char **argv)
{
    int myrank, nproc;
    int rows; /*amount of work per node (rows per worker)*/
    int mtype; /*message type: send/recv between master and workers*/
    int dest, src, offseta, offsetb;
    int runthrough, runmod;
    double start_time, end_time;
    int i, j, k, l;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    rows = SIZE/nproc;
    mtype = FROM_MASTER;

    if (myrank == 0) {
        /*Initialization*/
        printf("SIZE = %d, number of nodes = %d\n", SIZE, nproc);
        init_matrix();
        transpose_matrix();
        start_time = MPI_Wtime();

        if(nproc == 1) { /*In case we only run on one processor, the master will simply do a regular matrix-matrix multiplacation.*/
            for(i = 0; i < SIZE; i++) {
                for(j = 0; j < SIZE; j++) {
                    for(k = 0; k < SIZE; k++)
                        c[i][j] = c[i][j] + a[i][k]*b[j][k];
                }
            }
            end_time = MPI_Wtime();
            if(DEBUG) /*Prints the resulting matrix c*/
                print_matrix();
            printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);
        }
        else {

            for(l = 0; l < nproc; l++){
                offsetb = rows*l;
                offseta = rows;
                mtype = FROM_MASTER;

                for(dest = 1; dest < nproc; dest++){
                    MPI_Send(&offseta, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                    MPI_Send(&offsetb, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                    MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                    MPI_Send(&a[offseta][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                    MPI_Send(&b[offsetb][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
                    offseta += rows;
                    offsetb = (offsetb+rows)%SIZE;
                }

                offseta = rows;
                offsetb = rows*l;
                //printf("Rank: %d, offseta: %d, offsetb: %d\n", myrank, offseta, offsetb);
                //printf("Offseta: %d\n", offseta);
                //printf("Offsetb: %d\n", offsetb);
                for(i = 0; i < offseta; i++) {
                    for(j = offsetb; j < offsetb+rows; j++) {
                            for(k = 0; k < SIZE; k++){
                                c[i][j] = c[i][j] + a[i][k]*b[j][k];
                        }
                    }
                }
                mtype = FROM_WORKER;
                for(src = 1; src < nproc; src++){
                    MPI_Recv(&offseta, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                    MPI_Recv(&offsetb, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                    MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
                    for(i = 0; i < rows; i++) {
                        MPI_Recv(&c[offseta+i][offsetb], offseta, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status); /*returns answer c(1,1)*/
                    }
                }
            }


            end_time = MPI_Wtime();
            if(DEBUG) /*Prints the resulting matrix c*/
                print_matrix();
            printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);
        }
    }
    else{
        if(nproc > 1) {
            for(l = 0; l < nproc; l++){
                mtype = FROM_MASTER;
                MPI_Recv(&offseta, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&offsetb, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&a[offseta][0], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&b[offsetb][0], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);

                for(i = offseta; i < offseta+rows; i++) {
                    for(j = offsetb; j < offsetb+rows; j++) {
                        for(k = 0; k < SIZE; k++){
                            c[i][j] = c[i][j] + a[i][k]*b[j][k];
                        }
                    }
                }

                mtype = FROM_WORKER;
                MPI_Send(&offseta, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                MPI_Send(&offsetb, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
                for(i = 0; i < rows; i++){
                    MPI_Send(&c[offseta+i][offsetb], offseta, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
                }
            }
        }
    }
    MPI_Finalize();
    return 0;
}