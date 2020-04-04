#include <stdio.h>  
#include <stdlib.h>  
extern "C"{
int foo(int a, int b)  
{  
  printf("you input %d and %d\n", a, b);  
  return a+b;  
}  
}
