#include<odesolver.h>
#include<odesolver2.h>



int main()
{
   OdeSolver testcase;
   testcase.verlet();
   testcase.rk4();
   return 0;
}
