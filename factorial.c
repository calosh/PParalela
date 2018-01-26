//Calcular el factorial de un numero www.masterhacks.net
#include <stdio.h>

int main()
{
  int a, b, fact = 1;

  printf("Escribe un numero para calcular su factorial\n");
  scanf("%d", &a);

  for (b = a; b > 1; b--){
    fact = fact * b;
}
  printf("El factorial de %d = %d\n", a, fact);
  getchar();
  getchar();
  return 0;
}

// https://masterhacks.net/programacion/calcular-el-factorial-de-un-numero-en-lenguaje-c/