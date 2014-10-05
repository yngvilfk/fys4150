#ifndef ODESOLVER_H
#define ODESOLVER_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>


class OdeSolver
{
   public:
      OdeSolver( double x_0, double y_0, double v0_x, double v0_y, int timesteps, double delta_t_in);
      OdeSolver(double x_0, double y_0, double v0_x, double v0_y, int timesteps);
      OdeSolver(double x_0, double y_0, double v0_x, double v0_y);
      ~OdeSolver();


      void rk4();
      void verlet();


   protected:
      void derivatives(double&, arma::Col<double>&, arma::Col<double>&);
      void rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout);

   private:
      double delta_t;
      double x0;
      double y0;
      double v0x;
      double v0y;
      int n;






};

#endif // ODESOLVER_H
