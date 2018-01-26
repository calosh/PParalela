#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>

#define MATRIX_MSG 0
#define VECTOR_MSG 1
#define COLLECT_MSG 2
#define TIME_MSG 3

int me,
    nnodes, // number of MPI nodes in the computation
    chunk, // number of vertices handled by each node
    startv,endv; // start, end vertices for this node  // my node number ;
float T2,T1;
int *matrix,*vector,*result;
float *TNODE;

int randint(int n) {
        if ((n - 1) == RAND_MAX) {
                return rand();
        } else {
                // Chop off all of the values that would cause skew...
                long end = RAND_MAX / n; // truncate skew
                //assert (end > 0L);
                end *= n;

                // ... and ignore results from rand() that fall above that limit.
                // (Worst case the loop condition should succeed 50% of the time,
                // so we can expect to bail out of this loop pretty quickly.)
                int r;
                while ((r = rand()) >= end) ;

                return r % n;
        }
}

void init(int ac, char **av, int nv){
        MPI_Init(&ac,&av);
        MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
        MPI_Comm_rank(MPI_COMM_WORLD,&me);
        chunk = nv/nnodes;
        startv = me * chunk;
        endv = startv + chunk;
        matrix = malloc(nv*nv*sizeof(int));
        vector = malloc(nv*sizeof(int));
        result = malloc(nv*sizeof(int));
        TNODE = malloc(nnodes*sizeof(float));
        if(me==0) {
                for (int i = 0; i < nv; i++) {
                        for (int j = 0; j < nv; j++)  {
                                matrix[nv*i+j] = randint(nv);
                        }
                }
                for (int i = 0; i < nv; i++) {
                        vector[i]=randint(nv);
                }
        }
        for (int i = 0; i < nv; i++) {
                result[i]=0;
        }
        for (int i = 0; i < nnodes; i++) {
                TNODE[i]=0;
        }
}


void dowork(int nv, int dbg){
        MPI_Status status;
        int *result_temp;
        int *tnode_temp;
        result_temp = malloc(nv*sizeof(int));
        tnode_temp = malloc(nnodes*sizeof(float));
        if(me == 0) {
                for(int i=1; i <nnodes; i++) {
                        MPI_Send(matrix+nv*(i*chunk), chunk*nv, MPI_INT, i, MATRIX_MSG, MPI_COMM_WORLD);
                        MPI_Send(vector,nv, MPI_INT, i, VECTOR_MSG, MPI_COMM_WORLD);
                }
        }else{
                MPI_Recv(matrix+nv*(me*chunk), chunk*nv,MPI_INT,0,MATRIX_MSG,MPI_COMM_WORLD,&status);
                MPI_Recv(vector,nv,MPI_INT,0,VECTOR_MSG,MPI_COMM_WORLD,&status);
        }
        if(me==0&&dbg) {
                printf("%d-%d-%d\n", startv,endv, nnodes);
        }
        T1 = MPI_Wtime();
        for(int i=startv; i < endv; i++) {
                result[i] = 0;
                int temp_result = 0;
                for(int j=0; j<nv; j++) {
                        if(me==0&&dbg) {
                                printf("%d %d %d\n", matrix[nv*i+j], vector[j], (matrix[nv*i+j]*vector[j]));
                        }
                        temp_result+=(matrix[nv*i+j]*vector[j]);
                }
                result[i] =  temp_result;
                if(me==0&&dbg) {
                        printf("%d\n", result[i]);
                }
        }
        T2=MPI_Wtime();
        if(dbg) {
                printf("%.15f : %.15f : %.15f\n", T2,T1,T2-T1);
        }
        TNODE[me] = T2-T1;
        if(me!=0) {
                MPI_Send(result+startv,chunk,MPI_INT,0,COLLECT_MSG,MPI_COMM_WORLD);
                MPI_Send(TNODE+me,1,MPI_FLOAT,0,TIME_MSG,MPI_COMM_WORLD);
        }else{
                for(int i = 1; i< nnodes; i++) {
                        MPI_Recv(result+i*chunk,chunk,MPI_INT,i,COLLECT_MSG,MPI_COMM_WORLD,&status);
                        MPI_Recv(TNODE+i,1,MPI_FLOAT,i,TIME_MSG,MPI_COMM_WORLD,&status);
                }
        }
}

int main(int ac, char **av)
{
        int nv = atoi(av[1]);
        int print = atoi(av[2]);
        int dbg = atoi(av[3]);
        srand(1000);
        init(ac,av,nv);
        dowork(nv,dbg);
        if (print && me == 0)  {
                printf("graph weights:\n");
                for (int i = 0; i < nv; i++)  {
                        for (int j = 0; j < nv; j++)
                                printf("%u  ",matrix[nv*i+j]);
                        printf("\n");
                }
                printf("vector weights:\n");
                for (int i = 0; i < nv; i++)  {
                        printf("%u\n",vector[i]);
                }
                printf("result weights:\n");
                for (int i = 0; i < nv; i++)  {
                        printf("%u\n",result[i]);
                }
                printf("time spends:\n");
                for (int i = 0; i < nnodes; i++)  {
                        printf("%d : %.15f\n",i,TNODE[i]);
                }
        }
        //if (me == 0) printf("time at node 0: %f\n",(float)(T2-T1));
        MPI_Finalize();
}