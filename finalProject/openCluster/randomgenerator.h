#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include<system.h>
#include<object.h>
#include<math.h>
#include<cmath>

class RandomGenerator
{
   public:
      RandomGenerator();

      System randomSystem(int N, double R0);

      // ran2 for uniform deviates, initialize with negative seed.
      double ran2(long &idum);

      // function for gaussian random numbers
      double gaussian_deviate(long & idum);
};

#endif // RANDOMGENERATOR_H
