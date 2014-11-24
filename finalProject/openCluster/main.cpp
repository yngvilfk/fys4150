#include<randomgenerator.h>
#include<system.h>
#include<object.h>
#include<solvestep.h>


int main()
{
   double PI = 4*std::atan(1.0);

//-------------------------------------------------------------
   /* problem a)
    * additional functionality as parallelisation
    * and adaptive timesteps in other classes was commented out
    * while generating plots and time for this part*/
//-------------------------------------------------------------
   //initialize sun
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

   positionEarth(0) = 1.0;
   velocityEarth(1) = -2*PI;


   //Declaration of objects
   Object earth(positionEarth,velocityEarth, massEarth, nameEarth);
   Object sun(positionSun, velocitySun, massSun, nameSun);

//   //add objects to system
//   solarsystem.addObject(sun);
//   solarsystem.addObject(earth);

//   std::string nameAdd = "rk4";
//   earth.newFile(nameAdd);
//   sun.newFile(nameAdd);
//   SolveStep solver(solarsystem);

//   solver.solve(0.01,0.05,"rk4_","rk4");

//   earth.closeFile(nameAdd);
//   sun.closeFile(nameAdd);

//-------------------------------------------------------------
   /* problem b)
    * additional functionality in other classes was commented out
    * while generating plots for this part*/
//-------------------------------------------------------------

   //reset sun
   positionSun.zeros();
   velocitySun.zeros();
   earth.setPosition(positionSun);
   earth.setVelocity(velocitySun);

   //reset earth
   positionEarth.zeros();
   velocityEarth.zeros();
   positionEarth(0) = 1.0;
   velocityEarth(1) = -2*PI;
   earth.setPosition(positionEarth);
   earth.setVelocity(velocityEarth);

   //initialize moon
   std::string nameMoon = "moon";
   arma::Col<double> positionMoon(3), velocityMoon(3);
   double massMoon = 3e-6*0.0123; //solarmass
   positionMoon.zeros();
   velocityMoon.zeros();

//   positionMoon(0) = 1.0;
//   positionMoon(1) = -0.00257;
//   velocityMoon(1) = -2*PI;
//   velocityMoon(0) = 0.1668*PI;
   positionMoon(0) = 1.0;
   positionMoon(1) = -2e-3;
   velocityMoon(1) = -2*PI;
   velocityMoon(2) = 0.21544;
   Object moon(positionMoon, velocityMoon, massMoon, nameMoon);

   //declare solar system
   System solarsystemAdapt;

   //add objects to system
   solarsystemAdapt.addObject(sun);
   solarsystemAdapt.addObject(earth);
   solarsystemAdapt.addObject(moon);

   std::string nameAdd = "verlet_";
   earth.newFile(nameAdd);
   sun.newFile(nameAdd);
   moon.newFile(nameAdd);
   SolveStep solverAdapt(solarsystemAdapt);
   solverAdapt.solve((1.0/(365*24)),2.0,nameAdd,"verlet");

   earth.closeFile(nameAdd);
   sun.closeFile(nameAdd);
   moon.closeFile(nameAdd);
   //open cluster
//------------------------------------------------------------------------------------------------------
//   RandomGenerator generate;
//   double R0 = 20;
//   System mysystem = generate.randomSystem(100, R0);
//   mysystem.setEpsilon(0.2);
//   SolveStep solver(mysystem);
//   solver.solve(100, 1000, "solveR_", "rk4");

//   std::ofstream fout("start.m");
//   fout << "A = [";
//   for (int i = 0 ; i < mysystem.numberOfObject ; ++i)
//   {
//      Object &mainbody = mysystem.objectlist[i];

//      fout << mainbody.getPosition()(0) << "\t\t" << mainbody.getPosition()(1)
//           << "\t\t" << mainbody.getPosition()(2) << "\n";
//   }
//   fout << "] \n";
//   fout << "plot3(A(:,1), A(:,2),A(:,3), 'o')";
//   fout.close();
//   std::cout << mysystem.numberOfObject << std::endl;
}
