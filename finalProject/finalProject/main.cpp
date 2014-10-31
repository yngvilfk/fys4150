#include <armadillo>
#include<distance.h>
#include<object.h>
#include<odesolver.h>
#include<system.h>

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

//declare solar system
   System solarsystem;

//   positionSun(0) = -massEarth/(massSun+massEarth);
   positionEarth(0) = 1.0;//+positionSun(0);
   velocityEarth(1) = -2*PI;
//   velocitySun(1) = -(velocityEarth(1)*positionEarth(0)*massEarth)/(-positionSun(0)*massSun);


//Declaration of objects
   Object earth(positionEarth,velocityEarth, massEarth, nameEarth);
   Object sun(positionSun, velocitySun, massSun, nameSun);


//add objects to system
   solarsystem.addObject(sun);
   solarsystem.addObject(earth);



   //cant make the for-loop work, make it "manually":
//       Object tempsun = sun;
//       arma::Col<double> newpos(3);
//       newpos(0)=tempsun.position(0)-(earth.position(0)*earth.mass/(sun.mass+earth.mass));
//       newpos(1) = tempsun.position(1);
//       newpos(2) = tempsun.position(2);
//       sun.update(newpos,tempsun.velocity);
   std::cout << earth.position << std::endl;
   std::cout << sun.position << std::endl;
   std::cout << earth.velocity << std::endl;
   std::cout << sun.velocity << std::endl;


   std::string rk4Name = "rk4Time4.m";
   OdeSolver solveSystemRK4(solarsystem);
   solveSystemRK4.rk4(1.0, 1.0*365.0*24.0, rk4Name);

   std::string verletName = "verletTime4.m";
   OdeSolver solveSystemVerlet(solarsystem);
   solveSystemVerlet.verlet(1.0, 1.0*365.0*24.0, verletName);



//   OdeSolver mysolver;
//   mysolver.rk4();
return 0;
}

