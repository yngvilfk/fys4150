#ifndef ODESOLVER_H
#define ODESOLVER_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>

class OdeSolver
{
   public:
      OdeSolver();
      ~OdeSolver();

      void derivatives(double&, arma::Col<double>&, arma::Col<double>&);
      void rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout, double& delta_t);
      void rk4();





};

#endif // ODESOLVER_H
