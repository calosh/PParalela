#include<mpi.h>
#include <stdio.h>
#include <time.h>
#define TIME_THIS(X)

#define NRA 500                     //Numero de filas  en la matriz de A 
#define NCA 500                         //Numero de columnas  en la matriz de A 
#define NCB 500                         //Numero de columnas  en la matriz de B 
#define MASTER 0                    //Id. de tarea de la primera tarea 
#define FROM_MASTER 1           //el establecimiento de un tipo de mensaje 
#define FROM_WORKER 2           //el establecimiento de un tipo de mensaje 

int main(argc,argv)
int argc;
char *argv[];
{
int numtasks,                    /* número de tareas en la partición*/
    taskid,                         /* un identificador de tarea*/
    numworkers,                 /* número de tareas de trabajo */
    source,                     /* tarea Identificación de la fuente de mensaje */
    dest,                           /* id tarea de destino del mensaje */
    mtype,                       /* Tipo de mensaje */
    rows,                           /* filas de la matriz A enviada a cada trabjador*/
    averow, extra, offset,              /* para determinar filas enviadas a cada trabajador */
    i, j, k, y, rc;                 /* Varios */
long    a[NRA][NCA],                /* la matriz A, que se multiplica */
    b[NCA][NCB],                /* la matriz B que se   multiplica */
    c[NRA][NCB];                 /* resultado de la matriz C */

struct timespec ti,tf;

MPI_Status status;
   rc = MPI_Init(&argc,&argv);                  /* Inizializando MPI */
   rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);    /* Número de procesos */
   rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);      /* Identificador del proceso*/
   if (rc != 0)
      printf ("Error al inicializar el MPI y la obtención de información de la tarea ID \n");
   else
      printf ("task ID = %d\n", taskid);
   numworkers = numtasks-1;
/**************************** Tarea del Master ************************************/
// Inicializa las matrices
   if (taskid == MASTER)
   {
      printf("Number of worker tasks = %d\n",numworkers);

      for (i=0; i<NRA; i++)
         for (j=0; j<NCA; j++)
            a[i][j]= i*j;
      for (i=0; i<NCA; i++)
         for (j=0; j<NCB; j++)
            b[i][j]= i*j;

      // Envía los datos de la matriz con las tareas de los trabajadores
clock_gettime( CLOCK_REALTIME, &ti );           // Tiempo Inicial de la ejecución
      averow = NRA/numworkers;
      extra = NRA%numworkers;
      offset = 0;
      mtype = FROM_MASTER;
      for (dest=1; dest<=numworkers; dest++)
      {
         rows = (dest <= extra) ? averow+1 : averow;    
         printf("   sending %d rows to task %d\n",rows,dest);
         MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&a[offset][0], rows*NCA, MPI_LONG, dest, mtype,
                   MPI_COMM_WORLD);
         MPI_Send(&b, NCA*NCB, MPI_LONG, dest, mtype, MPI_COMM_WORLD);
         offset = offset + rows;
      }
     
      /* Esperar los resultados de todas las tareas de los trabajadores */
      mtype = FROM_WORKER;
      for (i=1; i<=numworkers; i++)
      {
         source = i;
         MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&c[offset][0], rows*NCB, MPI_LONG, source, mtype, MPI_COMM_WORLD, &status);
      }
    clock_gettime( CLOCK_REALTIME, &tf );       //Tiempo final de la ejecución


      /* print results 
      printf("Aquí está la matriz de resultados \n");
      for (i=0; i<NRA; i++)
      {
         printf("\n"); 
         for (j=0; j<NCB; j++) 
            printf("%6.2f   ", c[i][j]);
      } */

/* Calculo del tiempo que demora la multiplicación de Matrices */
printf( " Tiempo Final: %f\n", \
      (float) ( 1.0*(1.0*tf.tv_nsec - ti.tv_nsec*1.0)*1e-9 \
      + 1.0*tf.tv_sec - 1.0*ti.tv_sec ) ); \
   printf ("\n");
   }

/**************************** Tarea del Esclavo ************************************/
   if (taskid > MASTER)
   {
      mtype = FROM_MASTER;
      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&a, rows*NCA, MPI_LONG, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&b, NCA*NCB, MPI_LONG, MASTER, mtype, MPI_COMM_WORLD, &status);

// Multiplica las Matrices, a través de sumas, por ejemplo 2*4 = 2+2+2+2
      for (k=0; k<NCB; k++)
         for (i=0; i<rows; i++)
         {
            c[i][k] = 0;
            for (j=0; j<NCA; j++)
    for(y=0;y<b[k][j];y++){
        c[i][k] += a[i][j];}
         }
      mtype = FROM_WORKER;
      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&c, rows*NCB, MPI_LONG, MASTER, mtype, MPI_COMM_WORLD);
   }
   MPI_Finalize();
}
