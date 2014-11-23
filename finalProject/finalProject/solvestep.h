#ifndef SOLVESTEP_H
#define SOLVESTEP_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>
#include<system.h>
#include<object.h>
#include<cmath>

class SolveStep
{
public:
      SolveStep(System &mysystem);

      void solve(double, double);

      void writeToFile(int i);

      void rk4Step(double delta_t,
               int objectNumber,
               Object &mainbody);

      void verlet(double delta_t,
               int objectNumber);
      int size() const {return mysolarsystem_.numberOfObject;}

private:
      System mysolarsystem_;
};


#endif // SOLVESTEP_H
