#include "odesolver.h"

OdeSolver::OdeSolver()
{
}



OdeSolver::~OdeSolver()
{
}



OdeSolver::derivatives(double& R, arma::Col<double>& in,arma::Col<double>& out)
{
   out(0) = in(1); // out(0) = dx/dt = v_x or dy/dt = v_y
   out(1) = in(0)/(R*R*R);
}


OdeSolver::rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout, double& delta_t)
{
   /*
    R    - holds the distance between two objects
    yin  - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t
    yout - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t+delta_t
    */
   arma::col<double> k1(2), k2(2), k3(2), k4(2), y_k(2);
   //Callculation of k1
   OdeSolver::derivatives(R, yin, yout);
   k1(1) = yout(1)* delta_t;
   k1(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k1(0)*0.5;
   y_k(1) = yin(1)+k1(1)*0.5;

   //Callculation of k2
   OdeSolver::derivatives(R, y_k, yout);
   k2(1) = yout(1)* delta_t;
   k2(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k2(0)*0.5;
   y_k(1) = yin(1)+k2(1)*0.5;

   //Callculation of k3
   OdeSolver::derivatives(R, y_k, yout);
   k3(1) = yout(1)* delta_t;
   k3(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k3(0);
   y_k(1) = yin(1)+k3(1);

   //Callculation of k4
   OdeSolver::derivatives(R, y_k, yout);
   k4(1) = yout(1)* delta_t;
   k4(0) = yout(0)* delta_t;

   //calculate new values for x and v_x
   yout(0) = yin(0) + 1.0/6.0 * (k1(0) + 2.*k2(0) + 2.*k3(0) + k4(0));
   yout(1) = yin(1) + 1.0/6.0 * (k1(1) + 2.*k2(1) + 2.*k3(1) + k4(1));
}



OdeSolver::rk4()
{
   arma::Col<double> yout(2), y_h(2), xout(2), x_h(2);
   double t_h = 0.0;
   y_h(0) = y(0); //y
   y_h(1) = y(1); //v_y
   x_h(0) = x(0); //x
   x_h(1) = x(1); //x_y

   std::ofstream fout("rk4.m");
   fout.setf(ios::scientific);
   fout.precision(20);
   fout << "% [t x v_x y v_y]" << "\n";
   fout << "A = [ " ;
   for (int i = 1 ; i<=n ; ++i)
   {
      rk4_step(R, x_h, xout, delta_t_roof);
      rk4_step(R, y_h, yout, delta_t_roof);
      fout << i*delta_t << "\t\t" << xout(0) << "\t\t" << xout(1) << "\t\t" << yout(0) << "\t\t" << yout(1) << "\n";
      t_h += delta_t_roof;
      y_h(0) = yout(0);
      y_h(1) = yout(1);
   }
   fout << "]" << "\n\n";
   fout.close;
}
