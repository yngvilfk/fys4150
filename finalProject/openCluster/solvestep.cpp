#include <sstream>
#include <omp.h>
#include "solvestep.h"

SolveStep::SolveStep(System &mysystem)
   :mysolarsystem_(mysystem)
{
}

int
SolveStep::solve(double timestep,
                 double time,
                 std::string name,
                 std::string method)
{
   omp_set_num_threads(8);
   System newsystem = mysolarsystem_;
   int n = time/timestep;
   double minTime = timestep;
   double limit = 40;
   Distance d;
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = newsystem.objectlist[i];
      newsystem.acceleration(mainbody, i);

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
      name.append("t");
      t += timestep;
      std::string filename = name;
      filename.append(".m");
      std::ofstream fout(filename.c_str());
      fout << "A = [";
#pragma omp parallel for
      for (int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = newsystem.objectlist[i];
         double dt = timestep;
         int j = 0;
//         std::cout << mainbody.getAcceleration() << std::endl;
         while(dt > mainbody.maxTimestep())
         {
               dt = dt/2.0;

               j += 1;
         }
         for (int kk = 0 ; kk < std::pow(2,j); ++kk)
         {
            if (method == "rk4")
            {
               rk4Step(dt, i, mainbody);
            }
            else if(method == "verlet")
            {
               verlet(dt, i, mainbody);
            }
            else
            {
               std::cout << "method must be 'verlet' or 'rk4'" << std::endl;
//               return 1;
            }

            if(dt < minTime)
            {
               minTime = dt;
            }

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
//            i-=1;
         }
      }
      fout << "] \n";
      fout << "plot3(A(:,1), A(:,2),A(:,3), 'o')";
      fout.close();

      bound << t << "\t\t" << mysolarsystem_.boundTotalEnergy(limit) << "\n";

      energy << t << "\t\t" << mysolarsystem_.totalEnergy() << "\n";

   }
   energy << "] \n";
   energy << "plot(A(:,1), A(:,2))";
   energy.close();

   bound << "] \n";
   bound << "plot(A(:,1), A(:,2))";
   bound.close();
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem_.objectlist[i];
      mysolarsystem_.acceleration(mainbody, i);
   }
   std::cout << minTime << std::endl;
   return 0;
}

void
SolveStep::rk4Step(double delta_t,
               int objectNumber,
               Object &mainbody)
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
   tempSystem.acceleration(tempBody, i);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k2Pos = dPosdt*delta_t;
   k2Vel = dVeldt*delta_t;
   tempBody.setPosition(mainbody.getPosition() + k2Pos*0.5);
   tempBody.setVelocity( mainbody.getVelocity() + k2Vel*0.5);

   // k3
   tempSystem.acceleration(tempBody, i);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k3Pos = dPosdt*delta_t;
   k3Vel = dVeldt*delta_t;
   tempBody.setPosition( mainbody.getPosition() + k3Pos);
   tempBody.setVelocity( mainbody.getVelocity() + k3Vel);


   // k4
   tempSystem.acceleration(tempBody, i);

   dPosdt = tempBody.getVelocity();
   dVeldt = tempBody.getAcceleration();

   k4Pos = dPosdt*delta_t;
   k4Vel = dVeldt*delta_t;

   tempBody.setPosition( mainbody.getPosition() + 1.0/6.0 * (k1Pos + 2.*k2Pos + 2.*k3Pos + k4Pos));
   tempBody.setVelocity( mainbody.getVelocity() + 1.0/6.0 * (k1Vel + 2.*k2Vel + 2.*k3Vel + k4Vel));

   mainbody.setPosition(tempBody.getPosition());
   mainbody.setVelocity(tempBody.getVelocity());
   mysolarsystem_.acceleration(mainbody, i);

}

void
SolveStep::verlet(double delta_t,
                  int objectNumber,
                  Object &mainbody)
{
   int i = objectNumber;

   System tempSystem = mysolarsystem_;
   Object &tempBody = tempSystem.objectlist[i];

   tempBody.setPosition(mainbody.getPosition()+mainbody.getVelocity()*delta_t +
                        mainbody.getAcceleration()*delta_t*delta_t*0.5); //pos_i+1

   tempSystem.acceleration(tempBody, i);

   tempBody.setVelocity(mainbody.getVelocity() + 0.5*
                        (mainbody.getAcceleration()+tempBody.getAcceleration())*delta_t); //vel_i+1

   mainbody.setPosition(tempBody.getPosition());
   mainbody.setVelocity(tempBody.getVelocity());
   mysolarsystem_.acceleration(mainbody, i);



   }
