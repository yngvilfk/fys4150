#include "rk4.h"
#include <cmath>
#include <iostream>
#include<fstream>
#include <iomanip>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream ofile;

void derivatives(double t,
                 double *y,
                 double *dydt);
void initialize (double& initial_x,
                 double&initial_y,
                 int& n_steps);
void output (double t,
             double *y,
             double E0);
void runge_kutta_4 (double *y,
                    double *dydx,
                    int n,
                    double x,
                    double h,
                    double yout,
                    void(*dervis)(double, double *,double *));


RK4::RK4()
{
   //based on the algorithm at p 249 in compendium

}



RK4::runge_kutta_4(double *y,
                   double *dydx,
                   int n,
                   double x,
                   double h,
                   double yout,
                   void(*dervis)(double, double *,double *))
{
   int i;
   double xh, hh, h6
}
