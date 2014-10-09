#include "odesolver2.h"

OdeSolver2::OdeSolver2(double x_0,
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
}

OdeSolver2::OdeSolver2(double x_0, double y_0, double v0_x, double v0_y, int timesteps)
   : x0(x_0)
   , y0(y_0)
   , v0x(v0_x)
   , v0y(v0_y)
   , n(timesteps)
{
   delta_t = 0.0001;
}

OdeSolver2::OdeSolver2(double x_0, double y_0, double v0_x, double v0_y)
   : x0(x_0)
   , y0(y_0)
   , v0x(v0_x)
   , v0y(v0_y)
{
   delta_t = 0.0001;
   n = 100;
}



OdeSolver2::OdeSolver2()
{
   delta_t = 0.01;
   n = 1000;
   x0 = 1.0;
   y0 = 0.0;
   v0x = 0.0;
   v0y = -2.*3.14159;
   x0J = 5.2; //AU
   y0J = 0.0;
   v0xJ = 0.0;
   v0yJ = -2.*3.14159*0.44;
}


OdeSolver2::~OdeSolver2()
{
}


void
OdeSolver2::derivatives(double& R, arma::Col<double>& in,arma::Col<double>& out)
{
   out(0) = -in(1); // out(0) = dx/dt = v_x or dy/dt = v_y
   out(1) = 4*3.14159*3.14159*in(0)/(R*R*R);
}
void
OdeSolver2::derivativesEarth(double& R, arma::Col<double>& in, arma::Col<double>& out,double& posJupiter, double& RJE)
{
   out(0) = -in(1); // out(0) = dx/dt = v_x or dy/dt = v_y
   double indiff= std::sqrt((in(0)-posJupiter)*(in(0)-posJupiter));
   out(1) = 4*3.14159*3.14159*((in(0)/(R*R*R))+(indiff/(RJE*RJE*RJE*10)));
}
void
OdeSolver2::rk4_step(double R, arma::Col<double>& yin, arma::Col<double>& yout)
{
   /*
    R    - holds the distance between two objects
    yin  - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t
    yout - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t+delta_t
    */
   arma::Col<double> k1(2), k2(2), k3(2), k4(2), y_k(2);
   //Callculation of k1
   OdeSolver2::derivatives(R, yin, yout);
   k1(1) = yout(1)* delta_t;
   k1(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k1(0)*0.5;
   y_k(1) = yin(1)+k1(1)*0.5;

   //Callculation of k2
   OdeSolver2::derivatives(R, y_k, yout);
   k2(1) = yout(1)* delta_t;
   k2(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k2(0)*0.5;
   y_k(1) = yin(1)+k2(1)*0.5;

   //Callculation of k3
   OdeSolver2::derivatives(R, y_k, yout);
   k3(1) = yout(1)* delta_t;
   k3(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k3(0);
   y_k(1) = yin(1)+k3(1);

   //Callculation of k4
   OdeSolver2::derivatives(R, y_k, yout);
   k4(1) = yout(1)* delta_t;
   k4(0) = yout(0)* delta_t;

   //calculate new values for x and v_x
   yout(0) = yin(0) + 1.0/6.0 * (k1(0) + 2.*k2(0) + 2.*k3(0) + k4(0));
   yout(1) = yin(1) + 1.0/6.0 * (k1(1) + 2.*k2(1) + 2.*k3(1) + k4(1));
}

void
OdeSolver2::rk4_stepEarth(double R, arma::Col<double>& yin, arma::Col<double>& yout,double& posJupiter, double& RJE)
{
   /*
    R    - holds the distance between two objects
    yin  - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t
    yout - a two dimension array who holds the value of: yin(0)=x, yin(1)=dx/dt = v_x at time t+delta_t
    */
   arma::Col<double> k1(2), k2(2), k3(2), k4(2), y_k(2);
   //Callculation of k1
   OdeSolver2::derivativesEarth(R, yin, yout, posJupiter, RJE);
   k1(1) = yout(1)* delta_t;
   k1(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k1(0)*0.5;
   y_k(1) = yin(1)+k1(1)*0.5;

   //Callculation of k2
   OdeSolver2::derivativesEarth(R, y_k, yout, posJupiter, RJE);
   k2(1) = yout(1)* delta_t;
   k2(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k2(0)*0.5;
   y_k(1) = yin(1)+k2(1)*0.5;

   //Callculation of k3
   OdeSolver2::derivativesEarth(R, y_k, yout, posJupiter, RJE);
   k3(1) = yout(1)* delta_t;
   k3(0) = yout(0)* delta_t;
   y_k(0) = yin(0)+k3(0);
   y_k(1) = yin(1)+k3(1);

   //Callculation of k4
   OdeSolver2::derivativesEarth(R, y_k, yout, posJupiter, RJE);
   k4(1) = yout(1)* delta_t;
   k4(0) = yout(0)* delta_t;

   //calculate new values for x and v_x
   yout(0) = yin(0) + 1.0/6.0 * (k1(0) + 2.*k2(0) + 2.*k3(0) + k4(0));
   yout(1) = yin(1) + 1.0/6.0 * (k1(1) + 2.*k2(1) + 2.*k3(1) + k4(1));
}


void
OdeSolver2::rk4()
{
    arma::Col<double> yout(2), xout(2), youtJ(2), xoutJ(2), y(2), x(2), yJ(2), xJ(2);
    x(0) = x0;
    y(0) = y0;
    x(1) = v0x;
    y(1) = v0y;
    xJ(0) = x0J;
    yJ(0) = y0J;
    xJ(1) = v0xJ;
    yJ(1) = v0yJ;
   double t_h = 0.0;

   std::ofstream fout("rk4Two.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(16);
   fout << "describe = '[t xE yE xJ vJ]'" << "\n";
   fout << "A = [ " ;
   arma::Col<double> earth(2), sun(2), jupiter(2);
   sun(0) = 0.0;
   sun(1) = 0.0;
   for (int i = 1 ; i<=n ; ++i)
   {
       earth(0) = x(0);
       earth(1) = y(0);
       jupiter(0) = xJ(0);
       jupiter(1) = yJ(0);
       Distance dSunEarth;
      double RSunEarth = dSunEarth.twoObjects(earth,sun);
      Distance dSunJupiter;
     double RSunJupiter = dSunJupiter.twoObjects(jupiter,sun);
     Distance dEarthJupiter;
    double RJE = dEarthJupiter.twoObjects(jupiter,earth);

      OdeSolver2::rk4_stepEarth(RSunEarth, x, xout, xJ(0), RJE);
      OdeSolver2::rk4_stepEarth(RSunEarth, y, yout, yJ(0), RJE);

//      OdeSolver2::rk4_step(RSunEarth, x, xout);
//      OdeSolver2::rk4_step(RSunEarth, y, yout);

      OdeSolver2::rk4_step(RSunJupiter, xJ, xoutJ);
      OdeSolver2::rk4_step(RSunJupiter, yJ, youtJ);

      fout << i*delta_t << "\t\t" << xout(0) << "\t\t" << yout(0) << "\t\t" << xoutJ(0) << "\t\t" << youtJ(0) << "\n";
      t_h += delta_t;
      y(0) = yout(0);
      y(1) = yout(1);
      x(0) = xout(0);
      x(1) = xout(1);
      yJ(0) = youtJ(0);
      yJ(1) = youtJ(1);
      xJ(0) = xoutJ(0);
      xJ(1) = xoutJ(1);

   }
   fout << "];" << "\n\n";
   fout << "plot (A(:,2),A(:,3))" << "\n\n";
   fout << "hold on;" << "\n\n";
   fout << "plot (A(:,4),A(:,5))" << "\n\n";
   fout << "legend ('Earth','Jupiter')" << "\n\n";
   fout.close();
}



void
OdeSolver2::verlet()
{
   arma::Col<double> y(3), x(3), yin(2), xin(2), yout(2), xout(2);
   //startpos
   x(0) = 1.004;
   y(0) = -0.001;
   x(1) = x0;
   y(1) = y0;


//   //use rk4 to compute first step:
//   yin(0) = y(0);
//   yin(1) = v0y;
//   xin(0) = x(0);
//   xin(1) = v0x;
   arma::Col<double> earth(2), sun(2);
   sun(0) = 0.0;
   sun(1) = 0.0;
   earth(0) = x(0);
   earth(1) = y(0);
//   Distance d;
//   double R = d.twoObjects(earth,sun);
//   OdeSolver::rk4_step(R, xin, xout);
//   OdeSolver::rk4_step(R, yin, yout);
//   x(1) = xout(0);
//   y(1) = yout(0);
//   std::cout << "x: " << xout(1) << "y " << yout(1) << std::endl;

   std::ofstream fout("verlet.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(8);
   fout << "describe = '[t x y]'" << "\n";
   fout << "A = [ " ;
   fout << 0.0 << "\t\t" << x(0) << "\t\t" << y(0) << "\n";
   fout << 1*delta_t << "\t\t" << x(1) << "\t\t" << y(1) << "\n";


   for (int i = 2 ; i<=n ; ++i)
   {
       earth(0) = x(1);
       earth(1) = y(1);
       Distance d;
      const double R = d.twoObjects(earth,sun);
      std::cout << R << std::endl;
      x(2) = 2*x(1)-x(0)+(delta_t*delta_t*4.0*3.14159*3.14159*x(1)*x(1))/(R*R*R);
      y(2) = 2*y(1)-y(0)+(delta_t*delta_t*4.0*3.14159*3.14159*y(1)*y(1))/(R*R*R);
      fout << i*delta_t << "\t\t" << x(2) << "\t\t" << y(2) << "\n";
      x(0) = x(1);
      x(1) = x(2);
      y(0) = y(1);
      y(1) = y(2);
   }
   fout << "];" << "\n\n";
   fout << "plot (A(:,2),A(:,3))" << "\n\n";
   fout.close();
}
