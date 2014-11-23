#include "odesolver.h"

OdeSolver::OdeSolver(System &mysystem)
    : mysolarsystem(mysystem)
{
}


void
OdeSolver::rk4(double time,
                int nSteps,
                std::string filename)
{
   arma::Col<double> k1Pos(3), k1Vel(3), k2Pos(3), k2Vel(3), k3Pos(3), k3Vel(3), k4Pos(3), k4Vel(3);
   n = nSteps;
   delta_t = time/(n+1);


   std::ofstream fout(filename.c_str());
//   std::ofstream fout("test.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(16);
//   fout << "describe = '[t x y z Ek Ep AngMom]' runge-kutta4" << "\n";
//   fout << "numberOfObjects = " << mysolarsystem.numberOfObject << "\n";
//   fout << "timestep = " << delta_t << "\n";
//   fout << "timelength = " << time << "\n";
   fout << "A = [ " ;

   //set initial acceleration
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem.objectlist[i];

      mysolarsystem.acceleration(mainbody, i);
   }

   for ( int k = 0 ; k < n ; ++k )
   {

      for (int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = mysolarsystem.objectlist[i];
         System tempSystem = mysolarsystem;
         Object &tempBody = tempSystem.objectlist[i];

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
         mysolarsystem.acceleration(mainbody, i);

         double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
         double potentialEnergy = 0.0;
         arma::Col<double> angularMomentum = mysolarsystem.angularMomentum(mainbody);


         for ( int j = 0 ; j < size() ; ++j)
         {
            if (j != i)
            {
               Object body = mysolarsystem.objectlist[j];
               double energy = mysolarsystem.potentialEnergy(mainbody, body);
               potentialEnergy += energy;
            }

         } //end for

//         if (i == 0)
//         {
//            arma::Col<double> temp(3);
//            temp.zeros();
//            tempBody.setVelocity( temp);
//            tempBody.setPosition( temp);
//         }
         double ttotalEnergy = kineticEnergi+potentialEnergy;

         fout << k*delta_t<<"\t\t" << tempBody.getPosition()(0) << "\t\t" << tempBody.getPosition()(1)
               << "\t\t" << tempBody.getPosition()(2) << "\t\t" << ttotalEnergy << "\t\t"; //<<  angularMomentum << "\t\t";
         if (i == size() - 1)
         {
            fout << "\n";
         }
      }
   }
   fout << "];" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+1 << "),A(:," << 6*j+2 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << "')" << "\n\n";
//       fout << "hold on;" << "\n\n";

//   }

   fout.close();

}




void
OdeSolver::verlet(double time,
                   int nSteps,
                   std::string filename)
{
   n = nSteps;
   delta_t = time/(n+1);


   std::ofstream fout(filename.c_str());
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(8);
//   fout << "describe = '[t x y z Ek Ep AngMom]' Verlet method" << "\n";
//   fout << "numberOfObjects = " << mysolarsystem.numberOfObject << "\n";
//   fout << "timestep = " << delta_t << "\n";
//   fout << "timelength = " << time << "\n";
   fout << "A = [ " ;

   //set initial acceleration
   for (int i = 0 ; i < size() ; ++i)
   {
      Object &mainbody = mysolarsystem.objectlist[i];

      mysolarsystem.acceleration(mainbody, i);
   }

   for (int k = 1 ; k<=n ; ++k)
   {

      for ( int i = 0 ; i < size() ; ++i)
      {
         Object &mainbody = mysolarsystem.objectlist[i];
         System tempSystem = mysolarsystem;
         Object &tempBody = tempSystem.objectlist[i];

         tempBody.setPosition(mainbody.getPosition()+mainbody.getVelocity()*delta_t +
                              mainbody.getAcceleration()*delta_t*delta_t*0.5); //pos_i+1

         tempSystem.acceleration(tempBody, i);

         tempBody.setVelocity(mainbody.getVelocity() + 0.5*
                              (mainbody.getAcceleration()+tempBody.getAcceleration())*delta_t); //vel_i+1


         double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
         double potentialEnergy = 0.0;
         arma::Col<double> angularMomentum = mysolarsystem.angularMomentum(mainbody);


         for ( int j = 0 ; j < size() ; ++j)
         {
         if (j != i)
         {
            Object body = mysolarsystem.objectlist[j];
            double energy = mysolarsystem.potentialEnergy(mainbody, body);
            potentialEnergy += energy;
         }

         } //end for
         if (i == 0)
         {
            arma::Col<double> temp(3);
            temp.zeros();
            tempBody.setVelocity( temp);
            tempBody.setVelocity( temp);
         }
         double totalEnergy = potentialEnergy+kineticEnergi;

         fout << k*delta_t << "\t\t" << tempBody.getPosition()(0) << "\t\t" << tempBody.getPosition()(1) << "\t\t"
              << tempBody.getPosition()(2) << "\t\t"<< totalEnergy << "\t\t";// << angularMomentum << "\t\t";
         if (i == size() - 1)
         {
            fout << "\n";
         }

         mainbody.setPosition(tempBody.getPosition());
         mainbody.setVelocity(tempBody.getVelocity());
         mysolarsystem.acceleration(mainbody, i);

      }
   }
   fout << "];" << "\n\n";
   fout.close();
}
