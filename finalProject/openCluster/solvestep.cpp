#include <sstream>
#include <omp.h>
#include <time.h>
#include "solvestep.h"

SolveStep::SolveStep(System &mysystem)
   :mysolarsystem_(mysystem)
{
}

int
SolveStep::solve(double timestep,
                 double time,
                 std::string name,
                 std::string method,
                 std::string dimension)
{

//   clock_t start, finish;
//   omp_set_num_threads(8);
   System newsystem = mysolarsystem_;
   int n = time/timestep;
   double minTime = timestep;
   double limit = 200;
   Distance d;
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = newsystem.objectlist[i];
      newsystem.acceleration(mainbody, i, 0.0,dimension);

   }
   arma::Col<double> center(3), position(3);
   center.zeros();
   std::ofstream bound("BoundEnergy.m");
   bound << "A = [";
   double t = 0.;
   std::ofstream energy("energy.m");
   energy << "A = [";
   for (int k = 0 ; k < n ; ++k)
   {
      t += timestep;
      name.append("t");
      std::string filename = name;
      filename.append(".m");
      std::ofstream fout(filename.c_str());
      fout << "t = " << t << "\n A = [";
      std::cout << "t = " << t << std::endl;

      double addtime = 0.0;
//#pragma omp parallel for
      for (int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = newsystem.objectlist[i];
         double dt = timestep;
         int j = 0;
         double maxTimestep;

         if (mainbody.maxTimestep() < mysolarsystem_.maxTimestep(i) && mainbody.maxTimestep() > 1e-4)
         {
            maxTimestep = mainbody.maxTimestep();
         }
         else if (mainbody.maxTimestep() > mysolarsystem_.maxTimestep(i) && mysolarsystem_.maxTimestep(i) > 1e-4)
         {
            maxTimestep = mysolarsystem_.maxTimestep(i);
         }
         else
         {
            maxTimestep = 1e-2;
         }

         std::cout << maxTimestep << std::endl;
         while(dt > maxTimestep)
         {
               dt = dt/2.0;

               j += 1;
         }


         for (int kk = 0 ; kk < std::pow(2,j); ++kk)
         {
            if (method == "rk4")
            {
//               start = clock(); //start timer
               rk4Step(dt, i, mainbody, addtime, dimension);
//               finish = clock(); //stop timer
//               std::cout << "time one rk4 step: " <<
//                            static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << std::endl;
//               mainbody.addToFile(name);
               addtime += dt;
            }
            else if(method == "verlet")
            {
//               start = clock(); //start timer
               verlet(dt, i, mainbody, addtime, dimension);
//               finish = clock(); //stop timer
//               std::cout << "time one verlet step: " <<
//                            static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << std::endl;
//               mainbody.addToFile(name);
               addtime += dt;
            }
            else
            {
               std::cout << "method must be 'verlet' or 'rk4'" << std::endl;
            }

//            if(dt < minTime)
//            {
//               minTime = dt;
//            }
         }

      }
      mysolarsystem_ = newsystem;

      for (int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = mysolarsystem_.objectlist[i];
         position = mainbody.getPosition();
         double dist = d.twoObjects(position,center) ;
         if(dist < limit)
         {
            fout << mainbody.getPosition()(0) << "\t\t" << mainbody.getPosition()(1)
                 << "\t\t" << mainbody.getPosition()(2) << "\n";
         }
         else
         {
            mysolarsystem_.removeObject(i);
 //           i-=1;
         }
      }
      fout << "]; \n";
      fout << "plot3(A(:,1), A(:,2),A(:,3), 'o')";
      fout.close();

      bound << t << "\t\t" << mysolarsystem_.boundTotalEnergy(limit) << "\n";

      energy << t << "\t\t" << mysolarsystem_.totalEnergy() << "\n";

   }
   energy << "]; \n";
//   energy << "plot(A(:,1), A(:,2))";
   energy.close();

   bound << "]; \n";
//   bound << "plot(A(:,1), A(:,2))";
   bound.close();
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem_.objectlist[i];
      mysolarsystem_.acceleration(mainbody, i, 0.0,dimension);
   }
   std::cout << mysolarsystem_.numberOfObject <<"or: "<< size() << std::endl;
   return 0;
}

void
SolveStep::rk4Step(double delta_t,
                   int objectNumber,
                   Object &mainbody,
                   double addtime,
                   std::string dimension)
{
   int i = objectNumber;
   System tempSystem = mysolarsystem_;
   Object &tempBody = tempSystem.objectlist[i];
   arma::Col<double> k1Pos(3), k1Vel(3), k2Pos(3), k2Vel(3), k3Pos(3), k3Vel(3), k4Pos(3), k4Vel(3);

   // Calculate k1

   arma::Col<double> dPosdt, dVeldt;
   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k1Pos = dPosdt*delta_t;
   k1Vel = dVeldt*delta_t;

   tempBody.setPosition(mainbody.getPosition() + k1Pos*0.5);
   tempBody.setVelocity(mainbody.getVelocity() + k1Vel*0.5);

   // k2
   tempSystem.acceleration(tempBody, i, addtime,dimension);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k2Pos = dPosdt*delta_t;
   k2Vel = dVeldt*delta_t;
   tempBody.setPosition(mainbody.getPosition() + k2Pos*0.5);
   tempBody.setVelocity( mainbody.getVelocity() + k2Vel*0.5);

   // k3
   tempSystem.acceleration(tempBody, i, addtime,dimension);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k3Pos = dPosdt*delta_t;
   k3Vel = dVeldt*delta_t;
   tempBody.setPosition( mainbody.getPosition() + k3Pos);
   tempBody.setVelocity( mainbody.getVelocity() + k3Vel);


   // k4
   tempSystem.acceleration(tempBody, i, addtime,dimension);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k4Pos = dPosdt*delta_t;
   k4Vel = dVeldt*delta_t;

   tempBody.setPosition( mainbody.getPosition() + 1.0/6.0 * (k1Pos + 2.*k2Pos + 2.*k3Pos + k4Pos));
   tempBody.setVelocity( mainbody.getVelocity() + 1.0/6.0 * (k1Vel + 2.*k2Vel + 2.*k3Vel + k4Vel));

   mainbody.setPosition(tempBody.getPosition());
   mainbody.setVelocity(tempBody.getVelocity());
   mysolarsystem_.acceleration(mainbody, i, addtime,dimension);

}

void
SolveStep::verlet(double delta_t,
                  int objectNumber,
                  Object &mainbody,
                  double addtime,
                  std::string dimension)
{
   int i = objectNumber;

   System tempSystem = mysolarsystem_;
   Object &tempBody = tempSystem.objectlist[i];

   tempBody.setPosition(mainbody.getPosition()+mainbody.getVelocity()*delta_t +
                        mainbody.getAcceleration()*delta_t*delta_t*0.5); //pos_i+1

   tempSystem.acceleration(tempBody, i, addtime,dimension);

   tempBody.setVelocity(mainbody.getVelocity() + 0.5*
                        (mainbody.getAcceleration()+tempBody.getAcceleration())*delta_t); //vel_i+1

   mainbody.setPosition(tempBody.getPosition());
   mainbody.setVelocity(tempBody.getVelocity());
   mysolarsystem_.acceleration(mainbody, i, addtime,dimension);
   }
