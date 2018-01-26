#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include<mpi.h>

#include "blas.h"
#include "blacs.h"
#include "scalapack.h"

extern void pdsyev_( char *jobz, char *uplo, int *n,
                double *a, int *ia, int *ja, int *desca, double *w, double *z, i                                                                                   nt *iz, int *jz, int *descz,
                double *work, int *lwork, int *info );


static int max( int a, int b ){
        if (a>b) return(a); else return(b);
}
static int min( int a, int b ){
        if (a<b) return(a); else return(b);
}

int main(int argc, char **argv) {
        int iam, nprocs;
        int myrank_mpi, nprocs_mpi;
        int ictxt, nprow, npcol, myrow, mycol;
        int np, nq, nb, n;
        int mpA, nqA;
        int i, j, k, info, itemp, seed, lwork, min_mn;
        int descA[9], descZ[9];
        double *A, *Z, *work, *W;
        int izero=0,ione=1;
        double mone=(-1.0e0),pone=(1.0e0),dzero=(0.0e0);
/**/
        double MPIt1, MPIt2, MPIelapsed, GFLOPS, GFLOPS_per_proc ;
        char jobz, uplo;
/**/
        MPI_Init( &argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
/**/
        n = 100; nprow = 1; npcol = 1; nb = 64; jobz= 'V'; uplo='U';
        for( i = 1; i < argc; i++ ) {
                if( strcmp( argv[i], "-jobz" ) == 0 ) {
                        if (i+1<argc) {
                                if( strcmp( argv[i+1], "V" ) == 0 ){ jobz = 'V';                                                                                    i++; }
                                else if( strcmp( argv[i+1], "N" ) == 0 ){ jobz =                                                                                    'N'; i++; }
                                else if( strcmp( argv[i+1], "A" ) == 0 ){ jobz =                                                                                    'A'; i++; }
                                else printf(" ** warning: jobu should be set to                                                                                    V, N or A in the command line ** \n");
                        }
                        else
                                printf(" ** warning: jobu should be set to V, N                                                                                    or A in the command line ** \n");
                }
                if( strcmp( argv[i], "-n" ) == 0 ) {
                        n      = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-p" ) == 0 ) {
                        nprow  = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-q" ) == 0 ) {
                        npcol  = atoi(argv[i+1]);
                        i++;
                }
                if( strcmp( argv[i], "-nb" ) == 0 ) {
                        nb     = atoi(argv[i+1]);
                        i++;
                }
        }
/**/
        if (nb>n)
                nb = n;
        if (nprow*npcol>nprocs_mpi){
                if (myrank_mpi==0)
                        printf(" **** ERROR : we do not have enough processes av                                                                                   ailable to make a p-by-q process grid ***\n");
                        printf(" **** Bye-bye                                                                                                                                                            ***\n");
                MPI_Finalize(); exit(1);
        }
/**/
        Cblacs_pinfo( &iam, &nprocs ) ;
        Cblacs_get( -1, 0, &ictxt );
        Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
        Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
/**/
        if (iam==0)
                printf("\n");
/**/
/*
 *      if (iam==0)
 *              printf("\tn = %d\tnrhs = %d\t(%d,%d)\t%dx%d\n",n,nrhs,nprow,npco                                                                                   l,nb,nb);
 *      printf("Hello World, I am proc %d over %d for MPI, proc %d over %d for B                                                                                   LACS in position (%d,%d) in the process grid\n",
 *                      myrank_mpi,nprocs_mpi,iam,nprocs,myrow,mycol);
 */
/*
*
*     Work only the process in the process grid
*
*/
        if ((myrow < nprow)&(mycol < npcol)){
/*
*
*     Compute the size of the local matrices (thanks to numroc)
*
*/
                mpA    = numroc_( &n     , &nb, &myrow, &izero, &nprow );
                nqA    = numroc_( &n     , &nb, &mycol, &izero, &npcol );
/*
*
*     Allocate and fill the matrices A and B
*
*/
                A = (double *)calloc(mpA*nqA,sizeof(double)) ;
                if (A==NULL){ printf("error of memory allocation A on proc %dx%d                                                                                   \n",myrow,mycol); exit(0); }
                Z = (double *)calloc(mpA*nqA,sizeof(double)) ;
                if (Z==NULL){ printf("error of memory allocation VT on proc %dx%                                                                                   d\n",myrow,mycol); exit(0); }
                W = (double *)calloc(min_mn,sizeof(double)) ;
                if (W==NULL){ printf("error of memory allocation S on proc %dx%d                                                                                   \n",myrow,mycol); exit(0); }
/**/
                seed = iam*(mpA*nqA*2); srand(seed);
/**/
                k = 0;
                for (i = 0; i < mpA; i++) {
                        for (j = 0; j < nqA; j++) {
                                A[k] = ((double) rand()) / ((double) RAND_MAX) -                                                                                    0.5 ;
                                k++;
                        }
                }
/*
*
*     Initialize the array descriptor for the distributed matrices xA, U and VT
*
*/
                itemp = max( 1, mpA );
                descinit_( descA,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &it                                                                                   emp, &info );
                descinit_( descZ,  &n, &n, &nb, &nb, &izero, &izero, &ictxt, &it                                                                                   emp, &info );

                work = (double *)calloc(2,sizeof(double)) ;
                if (work==NULL){ printf("error of memory allocation for work on                                                                                    proc %dx%d (1st time)\n",myrow,mycol); exit(0); }
                lwork=-1;
                pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione,                                                                                    &ione, descZ, work, &lwork, &info );
                lwork= (int) work[0];
                free(work);
/**/
                work = (double *)calloc(lwork,sizeof(double)) ;
                if (work==NULL){ printf("error of memory allocation work on proc                                                                                    %dx%d\n",myrow,mycol); exit(0); }
/**/
                MPIt1 = MPI_Wtime();

                pdsyev_( &jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione,                                                                                    &ione, descZ, work, &lwork, &info );

                MPIt2 = MPI_Wtime();
                MPIelapsed=MPIt2-MPIt1;
/**/
                free(work);
                if ( iam==0 ){
                        printf("n=%d\t(%d,%d)\t%d\tjobz=%c\t%8.2fs \n",n,nprow,n                                                                                   pcol,nb,jobz,MPIelapsed);
                }
/**/
                free(W);
                free(Z);
                free(A);
        }
/*
*     Print ending messages
*/
        if ( iam==0 ){
                printf("\n");
        }
/**/
        Cblacs_gridexit( 0 );
        MPI_Finalize();
        exit(0);
}