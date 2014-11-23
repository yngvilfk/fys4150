#ifndef ODESOLVER_H
#define ODESOLVER_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>
#include<system.h>
#include<object.h>

class OdeSolver
{
public:
      OdeSolver(System &mysystem);


      void rk4(double time,
               int nSteps,
               std::string filename);
      void verlet(double time,
                  int nSteps,
                  std::string filename);
      int size() const {return mysolarsystem.numberOfObject;}

protected:


private:
      System mysolarsystem;
      double delta_t;
      int n;

};

#endif // ODESOLVER_H
