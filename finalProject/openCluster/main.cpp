#include<randomgenerator.h>
#include<system.h>
#include<object.h>
#include<solvestep.h>
#include <time.h>

int main()
{
   double PI = 4*std::atan(1.0);
   clock_t start, finish;
   std::string dimension = "AU"; //can be 'AU' or 'ly'

//-------------------------------------------------------------
   /* problem a)
    * additional functionality as parallelisation
    * and adaptive timesteps in other classes was commented out
    * while generating plots and time for this part*/
//-------------------------------------------------------------
   //initialize sun
//   std::string nameSun = "sun";
//   arma::Col<double> positionSun(3), velocitySun(3);
//   double massSun = 1.0; //SolarMass
//   positionSun.zeros();
//   velocitySun.zeros();

//   //initialize earth
//   std::string nameEarth = "earth";
//   arma::Col<double> positionEarth(3), velocityEarth(3);
//   double massEarth = 3e-6; //solarmass
//   positionEarth.zeros();
//   velocityEarth.zeros();

//   //declare solar system
//   System solarsystem;

//   positionEarth(0) = 1.0;
//   velocityEarth(1) = -2*PI;


//   //Declaration of objects
//   Object earth(positionEarth,velocityEarth, massEarth, nameEarth);
//   Object sun(positionSun, velocitySun, massSun, nameSun);

//   //add objects to system
//   solarsystem.addObject(sun);
//   solarsystem.addObject(earth);

//   std::string nameAdd = "rk4";
//   earth.newFile(nameAdd);
//   sun.newFile(nameAdd);
//   SolveStep solver(solarsystem);

//   solver.solve(0.01,0.05,"rk4_","rk4", dimension);

//   earth.closeFile(nameAdd);
//   sun.closeFile(nameAdd);

//-------------------------------------------------------------
   /* problem b)
    * additional functionality in other classes was commented out
    * while generating plots for this part*/
//-------------------------------------------------------------

//   //reset sun
//   positionSun.zeros();
//   velocitySun.zeros();
//   earth.setPosition(positionSun);
//   earth.setVelocity(velocitySun);

//   //reset earth
//   positionEarth.zeros();
//   velocityEarth.zeros();
//   positionEarth(0) = 1.0;
//   velocityEarth(1) = -2*PI;
//   earth.setPosition(positionEarth);
//   earth.setVelocity(velocityEarth);

//   //initialize moon
//   std::string nameMoon = "moon";
//   arma::Col<double> positionMoon(3), velocityMoon(3);
//   double massMoon = 3e-6*0.0123; //solarmass
//   positionMoon.zeros();
//   velocityMoon.zeros();

////   positionMoon(0) = 1.0;
////   positionMoon(1) = -0.00257;
////   velocityMoon(1) = -2*PI;
////   velocityMoon(0) = 0.1668*PI;
//   positionMoon(0) = 1.0;
//   positionMoon(1) = -2e-3;
//   velocityMoon(1) = -2*PI;
//   velocityMoon(2) = 0.21544;
//   Object moon(positionMoon, velocityMoon, massMoon, nameMoon);

//   //declare solar system
//   System solarsystemAdapt;

//   //add objects to system
//   solarsystemAdapt.addObject(sun);
//   solarsystemAdapt.addObject(earth);
//   solarsystemAdapt.addObject(moon);

//   std::string nameAdd = "verletAdapt_";
//   earth.newFile(nameAdd);
//   sun.newFile(nameAdd);
//   moon.newFile(nameAdd);
//   SolveStep solverAdapt(solarsystemAdapt);

//   double timestep = 1.0;

//   start = clock(); //start timer
//   solverAdapt.solve(timestep,4.0,nameAdd,"verlet", dimension);
//   finish = clock(); //stop timer
//   std::cout << "Simulation time: " <<
//                static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << std::endl;

//   earth.closeFile(nameAdd);
//   sun.closeFile(nameAdd);
//   moon.closeFile(nameAdd);

   //-------------------------------------------------------------
      /* problem c)
       * creation of a simple model of an open cluster. The lines
       * in Solvestep that saves every position for every particle
       * is commented out during these calculations. */
   //-------------------------------------------------------------
   dimension = "ly";
   RandomGenerator generate;
   double R0 = 20;
   double epsilon = 0.0;//6.0e-4;
   int N = 100;
   System mysystem = generate.randomSystem(N, R0);
   mysystem.setEpsilon(epsilon);
   SolveStep solver(mysystem);
   double timestep = 10.0;
   double simulationTime = 100.0;
   std::string fileName = "solve_";
   std::string method = "rk4";
   solver.solve(timestep, simulationTime, fileName, method, dimension);

   std::ofstream fout("start.m");
   fout << "A = [";
   for (int i = 0 ; i < mysystem.numberOfObject ; ++i)
   {
      Object &mainbody = mysystem.objectlist[i];

      fout << mainbody.getPosition()(0) << "\t\t" << mainbody.getPosition()(1)
           << "\t\t" << mainbody.getPosition()(2) << "\n";
   }
   fout << "] \n";
   fout << "plot3(A(:,1), A(:,2),A(:,3), 'o')";
   fout.close();

}
