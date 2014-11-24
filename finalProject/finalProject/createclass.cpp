#include "createclass.h"



CreateClass::CreateClass()
{
}


System
CreateClass::randomSystem(int N)
{
   double PI = 4*std::atan(1.0);
   System newSystem;
   long value = -2;
   long value2 = -3;
   for (std::size_t i = 0; i < N; ++i)
   {
      double theta = Distribution::ran2(value)*2*PI;
      double phi = Distribution::ran2(value2)*2*PI;
      Object newObject;
      arma::Col<double> positionSun(3), velocitySun(3);
      double massSun = 1.0; //SolarMass
   }
}
