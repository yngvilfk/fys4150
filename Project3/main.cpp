#include<odesolver.h>
#include<odesolver2.h>
#include<odesolver3.h>
#include<object.h>
#include<system.h>
#include<armadillo>
#include<distance.h>



int main()
{
   double PI = 4*std::atan(1.0);


//name of objects in my solar system
   std::string nameEarth = "earth";
   std::string nameSun = "sun";
   std::string nameJupiter = "jupiter";


//declare solar system
   System solarsystem;

//Declaration of planets
   Object earth(1.0,0.0,0.0,2.*PI, 3e-6, nameEarth);
   Object sun(0.0,0.0,0.0,0.0, 1.0, nameSun);
   Object jupiter(5.2,0.0,0.0,0.9*PI, 0.001, nameJupiter);



//add planets to system
   solarsystem.addObject(sun);
   solarsystem.addObject(earth);
   solarsystem.addObject(jupiter);



//simulate system using RK4 method
//   OdeSolver2 solveSystem(solarsystem);
//   solveSystem.rk4(4.0, 1000);


   //initial position of sun:
   arma::Col<double> sum(2);
   sum.zeros();
   for (int i = 1 ; i <= solarsystem.numberOfObject; ++i)
   {
       Object planet = solarsystem.objectlist[i];
       sum(0) += (planet.mass * planet.position(0));
       sum(1) += (planet.mass * planet.position(1));
   }

   Odesolver3 solveSystem(solarsystem);
   solveSystem.rk4(10, 10000);



//   OdeSolver mysolver;
//   mysolver.rk4();
return 0;
}
