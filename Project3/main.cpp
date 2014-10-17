#include<odesolver.h>
#include<odesolver2.h>
#include<object.h>
#include<system.h>
#include<armadillo>
#include<distance.h>



int main()
{
   double PI = 4*std::atan(1.0);

   std::string nameEarth = "earth";
   std::string nameSun = "sun";
   std::string nameJupiter = "jupiter";

   System solarsystem;


   Object earth(1.0,0.0,0.0,-2.*PI, 3e-6, nameEarth);
   Object sun(0.0,0.0,0.0,0.0, 1.0, nameSun);
   Object jupiter(5.2,0.0,0.0,0.9*PI, 0.001, nameJupiter);


   solarsystem.addObject(sun);
   solarsystem.addObject(earth);
//   solarsystem.addObject(jupiter);


   OdeSolver2 solves(solarsystem);

   solves.rk4(5.0, 15000);


   OdeSolver mysolver;
   mysolver.rk4();
return 0;
}
