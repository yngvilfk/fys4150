#include "odesolver.h"

OdeSolver::OdeSolver(double x_0,
                     double y_0,
                     double v0_x,
                     double v0_y,
                     int    timesteps,
                     double delta_t_in)
   : delta_t(delta_t_in)
   , x0(x_0)
   , y0(y_0)
   , v0x(v0_x)
   , v0y(v0_y)
   , n(timesteps)
{
    PI = 4*std::atan(1.0);
}

OdeSolver::OdeSolver(double x_0, double y_0, double v0_x, double v0_y, int timesteps)
   : x0(x_0)
   , y0(y_0)
   , v0x(v0_x)
   , v0y(v0_y)
   , n(timesteps)
{
   delta_t = 0.0001;
}

OdeSolver::OdeSolver(double x_0, double y_0, double v0_x, double v0_y)
   : x0(x_0)
   , y0(y_0)
   , v0x(v0_x)
   , v0y(v0_y)
{
   PI = 4*std::atan(1.0);
   delta_t = 0.0001;
   n = 100;
}



OdeSolver::OdeSolver()
{
   PI = 4*std::atan(1.0);
   delta_t = 0.001;
   n = 1000;
   x0 = 1.0;
   y0 = 0.0;
   v0x = 0.0;
   v0y = -2.*PI;
}


OdeSolver::~OdeSolver()
{
}


void
OdeSolver::derivatives(double& R, arma::Col<double>& in,arma::Col<double>& out)
{
   out(0) = -in(1); // out(0) = dx/dt = v_x or dy/dt = v_y
   out(1) = 4*PI*PI*in(0)/(R*R*R);
}

void
OdeSolver::rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout)
{
   /*
    R    - holds the distance between two objects
    yin  - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t
    yout - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t+delta_t
    */
   arma::Col<double> k1(2), k2(2), k3(2), k4(2), y_k(2);
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


void
OdeSolver::rk4()
{
    arma::Col<double> yout(2), xout(2), y(2), x(2);
    x(0) = x0;
    y(0) = y0;
    x(1) = v0x;
    y(1) = v0y;
   double t_h = 0.0;

   std::ofstream fout("rk4.m");
   //fout.setf(std::ios::scientific);
   fout.precision(8);
   //fout.width(8);
   fout << "describe = '[t x v_x y v_y]'" << "\n";
   fout << "A = [ " ;
   arma::Col<double> earth(2), sun(2);
   sun(0) = 0.0;
   sun(1) = 0.0;
   for (int i = 1 ; i<=n ; ++i)
   {
       earth(0) = x(0);
       earth(1) = y(0);
       Distance d;
      double R = d.twoObjects(earth,sun);
      OdeSolver::rk4_step(R, x, xout);
      OdeSolver::rk4_step(R, y, yout);
      fout << i*delta_t << "\t\t" << xout(0) << "\t\t" << xout(1) << "\t\t" << yout(0) << "\t\t" << yout(1) << "\n";
      t_h += delta_t;
      y(0) = yout(0);
      y(1) = yout(1);
      x(0) = xout(0);
      x(1) = xout(1);
      std::cout << "x: " << xout(1) << "y " << yout(1) << std::endl;

   }
   fout << "];" << "\n\n";
   fout << "plot (A(:,2),A(:,4))" << "\n\n";
   fout.close();
}
void
OdeSolver::verlet()
{
   arma::Col<double> y(2), x(2), vy(2), vx(2);
   //startpos
   x(0) = x0;
   y(0) = y0;
   vx(0) = v0x;
   vy(0) = v0y;


   arma::Col<double> earth(2), sun(2);
   sun(0) = 0.0;
   sun(1) = 0.0;
   earth(0) = x(0);
   earth(1) = y(0);


   std::ofstream fout("verlet1.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(8);
   fout << "describe = '[t x y]'" << "\n";
   fout << "A = [ " ;
   fout << 0.0 << "\t\t" << x(0) << "\t\t" << y(0) << "\n";



   for (int i = 1 ; i<=n ; ++i)
   {
       earth(0) = x(0);
       earth(1) = y(0);
       Distance d;
      const double R = d.twoObjects(earth,sun);
      std::cout << R << std::endl;

      x(1) = x(0)+vx(0)*delta_t - 0.5*delta_t*delta_t*4*PI*PI*x(0)/(R*R*R); //x_i+1
      y(1) = y(0)+vy(0)*delta_t - 0.5*delta_t*delta_t*4*PI*PI*y(0)/(R*R*R); //y_i+1
      vx(1) = vx(0) - 0.5*(4*PI*PI*x(0)/(R*R*R)+4*PI*PI*x(1)/(R*R*R))*delta_t; //vx_i+1
      vy(1) = vy(0) - 0.5*(4*PI*PI*y(0)/(R*R*R)+4*PI*PI*y(1)/(R*R*R))*delta_t; //vy_i+1

      fout << i*delta_t << "\t\t" << x(1) << "\t\t" << y(1) << "\n";
      x(0) = x(1);
      vx(0) = vx(1);
      y(0) = y(1);
      vy(0) = vy(1);


   }
   fout << "];" << "\n\n";
   fout << "plot (A(:,2),A(:,3))" << "\n\n";
   fout.close();
}
