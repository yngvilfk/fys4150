#include <sstream>
#include "solvestep.h"

SolveStep::SolveStep(System &mysystem)
   :mysolarsystem_(mysystem)
{
}

void
SolveStep::solve(double timestep,
                 double time)
{
   int n = time/timestep;
   for (int i = 0 ; i < size() ; ++i)
   {
//      std::ostringstream os;
//      os << filename_ << "_" << i;
      Object &mainbody = mysolarsystem_.objectlist[i];
      mysolarsystem_.acceleration(mainbody, i);
      mainbody.newFile();

//      filename_.append(os.str());
//      std::ofstream filename(os.str());
   }
   for (int k = 0 ; k < n ; ++k)
   {

      for (int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = mysolarsystem_.objectlist[i];
         double dt = timestep;
         int j = 0;
         while(dt > mainbody.maxTimestep())
         {
               dt = dt/2.0;
               j += 1;
         }
         for (int kk = 0 ; kk < std::pow(2,j); ++kk)
         {
            rk4Step(dt, i, mainbody);
            mainbody.addToFile();

         }

      }
   }
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem_.objectlist[i];
      mysolarsystem_.acceleration(mainbody, i);
      mainbody.closeFile();
//      std::string file = filename_.append(static_cast<char>(i));
//      file.close();
   }

}

void
SolveStep::writeToFile(int i)
{
//   Object &mainbody = mysolarsystem_.objectlist[i];
//   std::string file = filename_.append(static_cast<std::string>(i));
//   file << mainbody.getPosition()(0) << "\t\t" << mainbody.getPosition()(1)
//           << "\t\t" << mainbody.getPosition()(2) << "\n";
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
                  int objectNumber)
{
   int i = objectNumber;
   int n = 10;
   //set initial acceleration
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem_.objectlist[i];

      mysolarsystem_.acceleration(mainbody, i);
   }

   for (int k = 1 ; k<=n ; ++k)
   {

      for ( int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = mysolarsystem_.objectlist[i];
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

      } //end for


   }
}
