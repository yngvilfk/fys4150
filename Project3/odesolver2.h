#ifndef ODESOLVER2_H
#define ODESOLVER2_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>
#include<system.h>
#include<object.h>



class OdeSolver2
{
   public:
      OdeSolver2( System &mysystem);
      ~OdeSolver2();

      void rk4(double time,
               int nSteps);
      void verlet(double time,
                  int nSteps);

      System mysolarsystem;

   protected:


   private:
      double delta_t;

      int n;






};

#endif // ODESOLVER2_H
