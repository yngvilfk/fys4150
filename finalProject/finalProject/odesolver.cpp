#include "odesolver.h"

OdeSolver::OdeSolver(System &mysystem)
    : mysolarsystem(mysystem)
{
}

OdeSolver::~OdeSolver()
{
}

//OdeSolver::rk4( std::string& fileName)
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


   for ( int k = 0 ; k < n ; ++k )
   {

      for (int i = 0 ; i < mysolarsystem.numberOfObject ; ++i)
      {

         Object &mainbody = mysolarsystem.objectlist[i];
         Object tempBody = mainbody;
         arma::Col<double> force(3);
         force.zeros();

         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {

            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force += F;
            }

         } //end for



         tempBody.position = mainbody.position;


         // Calculate k1
         arma::Col<double> dPosdt, dVeldt;
         dPosdt = tempBody.velocity;
         dVeldt = force;


         k1Pos = dPosdt*delta_t;
         k1Vel = dVeldt*delta_t;


        //er dette det samme som tempBody.position=mainpody.position+k1Pos*0.5
         tempBody.position = mainbody.position + k1Pos*0.5;
         tempBody.velocity = mainbody.velocity + k1Vel*0.5;

         // k2
         force.zeros();

         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {

            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force += F;
            }

         } //end for

         dPosdt = tempBody.velocity;
         dVeldt = force;

         k2Pos = dPosdt*delta_t;
         k2Vel = dVeldt*delta_t;
         tempBody.position = mainbody.position + k2Pos*0.5;
         tempBody.velocity = mainbody.velocity + k2Vel*0.5;

         // k3
         force.zeros();
         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {

            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force += F;
            }
         } //end for


         dPosdt = tempBody.velocity;
         dVeldt = force;

         k3Pos = dPosdt*delta_t;
         k3Vel = dVeldt*delta_t;
         tempBody.position = mainbody.position + k3Pos;
         tempBody.velocity = mainbody.velocity + k3Vel;


         // k4
         force.zeros();
         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {
            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force += F;
            }
         } //end for

         dPosdt = tempBody.velocity;
         dVeldt = force;

         k4Pos = dPosdt*delta_t;
         k4Vel = dVeldt*delta_t;

         tempBody.position = mainbody.position + 1.0/6.0 * (k1Pos + 2.*k2Pos + 2.*k3Pos + k4Pos);
         tempBody.velocity = mainbody.velocity + 1.0/6.0 * (k1Vel + 2.*k2Vel + 2.*k3Vel + k4Vel);

         mainbody.update(tempBody.position, tempBody.velocity);

         double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
         double potentialEnergy = 0.0;
         arma::Col<double> angularMomentum = mysolarsystem.angularMomentum(mainbody);


         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
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
            tempBody.velocity = temp;
            tempBody.position = temp;
         }
         double ttotalEnergy = kineticEnergi+potentialEnergy;

         fout << k*delta_t<<"\t\t" << tempBody.position(0) << "\t\t" << tempBody.position(1) << "\t\t" << tempBody.position(2) << "\t\t" << ttotalEnergy << "\t\t"; //<<  angularMomentum << "\t\t";
         if (i == mysolarsystem.numberOfObject - 1)
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



   for (int k = 1 ; k<=n ; ++k)
   {

       for ( int i = 0 ; i < mysolarsystem.numberOfObject ; ++i)
       {
          arma::Col<double> force(3);
          Object &mainbody = mysolarsystem.objectlist[i];
          Object tempBody = mainbody;
          force.zeros();

          for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
          {
             if (j != i)
             {
                arma::Col<double> F(3);
                Object body = mysolarsystem.objectlist[j];
                mysolarsystem.force(tempBody, body,F);
                force += F;
             }
          } //end for

      tempBody.position = mainbody.position+mainbody.velocity*delta_t +force*delta_t*delta_t*0.5; //pos_i+1

      arma::Col<double> newforce(3);
      newforce.zeros();

      for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
      {

         if (j != i)
         {
            arma::Col<double> F(3);
            Object body = mysolarsystem.objectlist[j];
            mysolarsystem.force(tempBody, body,F);
            newforce += F;
         }
      } //end for

      tempBody.velocity = mainbody.velocity + 0.5*(force+newforce)*delta_t; //vel_i+1


      double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
      double potentialEnergy = 0.0;
      arma::Col<double> angularMomentum = mysolarsystem.angularMomentum(mainbody);


      for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
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
         tempBody.velocity = temp;
         tempBody.position = temp;
      }
    double totalEnergy = potentialEnergy+kineticEnergi;

      fout << k*delta_t << "\t\t" << tempBody.position(0) << "\t\t" << tempBody.position(1) << "\t\t" << tempBody.position(2) << "\t\t"<< totalEnergy << "\t\t";// << angularMomentum << "\t\t";
      if (i == mysolarsystem.numberOfObject - 1)
      {
         fout << "\n";
      }

      mainbody.update(tempBody.position, tempBody.velocity);



    }
   }
   fout << "];" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+2 << "),A(:," << 6*j+3 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << "')" << "\n\n";
//       fout << "hold on;" << "\n\n";
//   }

//   fout << "figure()" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+1 << "),A(:," << 6*j+4 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << " kinetic energy" <<"')" << "\n\n";
//       fout << "hold on;" << "\n\n";
//   }
//    fout << "figure" << "\n\n";

//    for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//    {
//        Object thisobject = mysolarsystem.objectlist[j];
//        fout << "plot (A(:," << 6*j+1 << "),A(:," << 6*j+5 <<"))" << "\n\n";
//        fout << "legend('" << thisobject.name << " potential energy" <<"')" << "\n\n";
//        fout << "hold on;" << "\n\n";
//    }
//     fout << "figure" << "\n\n";

//     for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//     {
//         Object thisobject = mysolarsystem.objectlist[j];
//         fout << "plot (A(:," << 6*j+1 << "),A(:," << 6*j+6 <<"))" << "\n\n";
//         fout << "legend('" << thisobject.name << " angular momentum" <<"')" << "\n\n";
//         fout << "hold on;" << "\n\n";
//     }
//      fout << "figure" << "\n\n";
   fout.close();
}
