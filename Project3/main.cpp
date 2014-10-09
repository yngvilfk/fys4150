#include<odesolver.h>
#include<odesolver2.h>



int main()
{
   OdeSolver2 testcase;
   //testcase.verlet();
   testcase.rk4();
   return 0;
}
