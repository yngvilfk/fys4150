#include<armadillo>
#include<distance.h>
#include<object.h>
#include<odesolver.h>
#include<system.h>
#include<solvestep.h>
#include<gaussiandeviate.h>
#include<createclass.h>


int main()
{
   double PI = 4*std::atan(1.0);

//   initialize sun
   std::string nameSun = "sun";
   arma::Col<double> positionSun(3), velocitySun(3);
   double massSun = 1.0; //SolarMass
   positionSun.zeros();
   velocitySun.zeros();



   //initialize earth
   std::string nameEarth = "earth";
   arma::Col<double> positionEarth(3), velocityEarth(3);
   double massEarth = 3e-6; //solarmass
   positionEarth.zeros();
   velocityEarth.zeros();

   //initialize moon
   std::string nameMoon = "moon";
   arma::Col<double> positionMoon(3), velocityMoon(3);
   double massMoon = 3e-6*0.012; //solarmass
   positionMoon.zeros();
   velocityMoon.zeros();

//declare solar system
   System solarsystem;

//   positionSun(0) = -massEarth/(massSun+massEarth);
   positionEarth(0) = 1.0;//+positionSun(0);
   velocityEarth(1) = -2*PI;
   positionMoon(0) = 1.0;
   positionMoon(1) = -0.00257;
   velocityMoon(1) = -2*PI;
   velocityMoon(0) = 0.1668*PI;
//   velocitySun(1) = -(velocityEarth(1)*positionEarth(0)*massEarth)/(-positionSun(0)*massSun);


//Declaration of objects
   Object earth(positionEarth,velocityEarth, massEarth, nameEarth);
   Object sun(positionSun, velocitySun, massSun, nameSun);
   Object moon(positionMoon, velocityMoon, massMoon, nameMoon);


//add objects to system
   solarsystem.addObject(sun);
   solarsystem.addObject(earth);
   solarsystem.addObject(moon);

   SolveStep solver(solarsystem);
   solver.solve(5.0,5.0);

//   System solarsystemRK4 = solarsystem;
//   System solarsystemVerlet = solarsystem;


//   double year = 80.0;
//   double timestep = year*12;

//   std::string rk4Name = "rk43D.m";
//   OdeSolver solveSystemRK4(solarsystemRK4);
//   solveSystemRK4.rk4(year, timestep, rk4Name);

//   std::string verletName = "verlet3D.m";
//   OdeSolver solveSystemVerlet(solarsystemVerlet);
//   solveSystemVerlet.verlet(year, timestep, verletName);

   //SolveStep solve(solarsystemRK4, rk4Name);

//std::size_t
   long value = -2;
   for ( int i = 0; i < 100; ++i)
   {
      double A = Distribution::gaussian_deviate(value);
      std::cout << A << ", " << value <<  std::endl;
   }


//   OdeSolver mysolver;
//   mysolver.rk4();
return 0;
}

