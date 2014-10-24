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
   std::string nameMars = "mars";
   std::string nameVenus = "venus";
   std::string nameSaturn = "saturn";


//declare solar system
   System solarsystem;

//Declaration of objects
   Object earth(1.0,0.0,0.0,2.0*PI, 3e-6, nameEarth);
   Object sun(0.0,0.0,0.0,0.0, 1.0, nameSun);
   Object jupiter(5.2,0.0,0.0,0.9*PI, 0.001, nameJupiter);
   Object mars(1.52,0.0,0.0,2.0*PI, 3.3e-7, nameMars);
   Object venus(0.72,0.0,0.0,0.72*PI,2.45e-6,nameVenus);
   Object saturn(9.54,0.0,0.0,0.8*PI, 2.75e-4, nameSaturn);



//add objects to system
   solarsystem.addObject(sun);
   solarsystem.addObject(earth);
   solarsystem.addObject(jupiter);
   solarsystem.addObject(mars);
   solarsystem.addObject(venus);
   solarsystem.addObject(saturn);



//simulate system using RK4 method
//   OdeSolver2 solveSystem(solarsystem);
//   solveSystem.rk4(4.0, 1000);


   //initial position of sun:
   arma::Col<double> sum(2);
   sum.zeros();
   for (int i = 1 ; i < solarsystem.numberOfObject; ++i)
   {
       Object planet = solarsystem.objectlist[i];
       sum(0) += (planet.mass * planet.position(0));
       sum(1) += planet.mass;
   }
   double initialpos = -sum(0)/(sum(1)+sun.mass);

//   for (int i = 0 ; i < solarsystem.numberOfObject; ++i)
//   {
//       Object &planet = solarsystem.objectlist[i];
//       Object &tempPlanet = planet;
//       arma::Col<double> newpos(2);
//       newpos(0)=tempPlanet.position(0)+initialpos;
//       newpos(1) = tempPlanet.position(1);
//       planet.update(newpos,tempPlanet.velocity);

//   }
   //cant make the for-loop work, make it "manually":
       Object tempsun = sun;
       arma::Col<double> newpos(2);
       newpos(0)=tempsun.position(0)+initialpos;
       newpos(1) = tempsun.position(1);
       sun.update(newpos,tempsun.velocity);


       Object &tempearth = earth;
       newpos(0)=tempearth.position(0)+initialpos;
       newpos(1) = tempearth.position(1);
       earth.update(newpos,tempearth.velocity);


       Object tempPlanet = jupiter;
       newpos(0)=tempPlanet.position(0)+initialpos;
       newpos(1) = tempPlanet.position(1);
       jupiter.update(newpos,tempPlanet.velocity);

       tempPlanet = mars;
       newpos(0)=tempPlanet.position(0)+initialpos;
       newpos(1) = tempPlanet.position(1);
       mars.update(newpos,tempPlanet.velocity);

       tempPlanet = venus;
       newpos(0)=tempPlanet.position(0)+initialpos;
       newpos(1) = tempPlanet.position(1);
       venus.update(newpos,tempPlanet.velocity);

       tempPlanet = saturn;
       newpos(0)=tempPlanet.position(0)+initialpos;
       newpos(1) = tempPlanet.position(1);
       saturn.update(newpos,tempPlanet.velocity);



   double sumMom = 0.0;
   for (int i = 1 ; i < solarsystem.numberOfObject; ++i)
   {
       Object &planet = solarsystem.objectlist[i];
       sumMom += planet.mass*planet.position(0)*planet.velocity(1);
   }

   double newVely= sumMom/sun.position(0)*sun.mass;
   arma::Col<double> newVel(2);
   newVel(0) = 0.0;
   newVel(1) = newVely;
   sun.update(sun.position,newVel);

   std::cout << sun.position << sun.velocity << std::endl;

//   std::string rk4Name = "rk4Time16twoplanets1.m";
//   OdeSolver2 solveSystemRK4(solarsystem);
//   solveSystemRK4.rk4(16.0, 16.0*365.0*24.0, rk4Name);

   std::string verletName = "newfinalverlet20.m";
   Odesolver3 solveSystemVerlet(solarsystem);
   solveSystemVerlet.rk4(20.0, 20.0*365.0*6, verletName);



//   OdeSolver mysolver;
//   mysolver.rk4();
return 0;
}
