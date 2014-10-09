#ifndef ODESOLVER2_H
#define ODESOLVER2_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>


class OdeSolver2
{
   public:
      OdeSolver2( double x_0, double y_0, double v0_x, double v0_y, int timesteps, double delta_t_in);
      OdeSolver2(double x_0, double y_0, double v0_x, double v0_y, int timesteps);
      OdeSolver2(double x_0, double y_0, double v0_x, double v0_y);
      OdeSolver2();
      ~OdeSolver2();


      void rk4();
      void verlet();


   protected:
      void derivatives(double&, arma::Col<double>&, arma::Col<double>&);
      void derivativesEarth(double&, arma::Col<double>&, arma::Col<double>&,double& posJupiter, double& RJE);
      void rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout);
      void rk4_stepEarth(double R, arma::Col<double>& yin, arma::Col<double>& yout,double& posJupiter, double& RJE);

   private:
      double delta_t;
      double x0;
      double y0;
      double v0x;
      double v0y;
      int n;
      double x0J;
      double y0J;
      double v0xJ;
      double v0yJ;






};

#endif // ODESOLVER2_H
