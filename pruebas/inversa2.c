 #include<stdio.h>  
 #include<stdlib.h>  
 #define Max 10  
 int main()  
 {  
   float a[Max][Max+1],t,det=1;  
   int i,j,k,N;  
   printf("Enter the number of unknowwns : ");  
   scanf("%d",&N);

    printf("\nEnter the elements of the augmented matrix :\n");  
   for(i=0;i<N;i++)  
     for(j=0;j<N+1;j++)  
       scanf("%f",&a[i][j]);  
   for(i=0;i<N;i++)  
     for(j=0;j<N;j++)  
       if(i!=j)  
         {  
           t=a[j][i]/a[i][i];  
           for(k=0;k<N+1;k++)  
             a[j][k]-=a[i][k]*t;  
         }  
   for(i=0;i<N;i++)  
     det*=a[i][i];  
   printf("\nDeterminant = %.4f\n",det);  
   if(det==0){  
     printf("\nThe matrix is singular .\n");  
     exit(1);  
     }  
   printf("\nThe Gauss-Jordan Matrix is :\n\n");  
   for(i=0;i<N;i++)  
   {  
     for(j=0;j<N+1;j++)  
       printf("%.4f ",a[i][j]);  
     printf("\n");  
   }  
   printf("\nThe solution is :\n\n");  
   for(i=0;i<N;i++)  
     printf("x[ %d]=%.4f\n",i+1,a[i][N]/a[i][i]);  
   return 0;  
 }  



// http://code.nepalijivan.com.np/2013/05/gauss-jordan-method-implementation-with.html