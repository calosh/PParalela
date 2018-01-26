#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// Ejemplo aprenderaprogramar.com
int main() {
srand(time(NULL)); //El mayordomo pone a girar la diana
int test = rand(); //Primer disparo del robot
int num=1+rand()%(10-1);
printf ("El numero prueba %d\n", num);
printf ("El numero aleatorio de test vale %d\n", test);
printf ("Otros numeros aleatorios son: %d, %d, %d\n",rand(),rand(),rand()) ;
printf("La constante RAND_MAX vale %d\n", RAND_MAX);
return 0;
}