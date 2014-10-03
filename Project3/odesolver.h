#ifndef ODESOLVER_H
#define ODESOLVER_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>

class OdeSolver
{
   public:
      OdeSolver();
      ~OdeSolver();

      derivatives(double&, arma::Col<double>&, arma::Col<double>&);
      rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout, double& delta_t);
      rk4();
};

#endif // ODESOLVER_H
