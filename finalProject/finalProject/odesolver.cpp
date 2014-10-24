#include "odesolver.h"

OdeSolver::OdeSolver(System &mysystem)
    : mysolarsystem(mysystem)
{
}

OdeSolver::~OdeSolver()
{
}

//Odesolver3::rk4( std::string& fileName)
void
OdeSolver::rk4(double time,
                int nSteps,
                std::string filename)
{
   arma::Col<double> k1x(2), k2x(2), k3x(2), k4x(2), k1y(2), k2y(2), k3y(2), k4y(2), k1z(2), k2z(2), k3z(2), k4z(2);
   arma::Col<double> k1Pos(3), k1Vel(3), k2Pos(3), k2Vel(3), k3Pos(3), k3Vel(3), k4Pos(3), k4Vel(3);
   n = nSteps;
   delta_t = time/(n+1);


   std::ofstream fout(filename.c_str());
//   std::ofstream fout("test.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(16);
   fout << "describe = '[objectNmbr x y vx vy objectNmbr x y vx vy etc]'" << "\n";
   fout << "numberOfObjects = " << mysolarsystem.numberOfObject << "\n";
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
               force(0) += F(0);
               force(1) += F(1);
               force(2) += F(2);
            }

         } //end for



         tempBody.position(0) = mainbody.position(0);
         tempBody.velocity(0) = mainbody.velocity(0);
         tempBody.position(1) = mainbody.position(1);
         tempBody.velocity(1) = mainbody.velocity(1);
         tempBody.position(2) = mainbody.position(2);
         tempBody.velocity(2) = mainbody.velocity(2);


         // Calculate k1

         double dxdt = tempBody.velocity(0);
         double dvxdt = force(0);

         double dydt = tempBody.velocity(1);
         double dvydt = force(1);

         double dzdt = tempBody.velocity(2);
         double dvzdt = force(2);


         k1Pos(0) = dxdt*delta_t;
         k1Vel(0) = dvxdt*delta_t;

         k1Pos(1) = dydt*delta_t;
         k1Vel(1) = dvydt*delta_t;

         k1Pos(2) = dzdt*delta_t;
         k1Vel(2) = dvzdt*delta_t;

        //er dette det samme som tempBody.position=mainpody.position+k1Pos*0.5
         tempBody.position(0) = mainbody.position(0) + k1Pos(0)*0.5;
         tempBody.velocity(0) = mainbody.velocity(0) + k1Vel(0)*0.5;

         tempBody.position(1) = mainbody.position(1) + k1Pos(1)*0.5;
         tempBody.velocity(1) = mainbody.velocity(1) + k1Vel(1)*0.5;

         tempBody.position(2) = mainbody.position(2) + k1Pos(2)*0.5;
         tempBody.velocity(2) = mainbody.velocity(2) + k1Vel(2)*0.5;


         // k2
         force.zeros();

         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {

            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force(0) += F(0);
               force(1) += F(1);
               force(2) += F(2);
            }

         } //end for

         dxdt = tempBody.velocity(0);
         dvxdt = force(0);
         dydt = tempBody.velocity(1);
         dvydt = force(1);

         k2Pos(0) = dxdt*delta_t;
         k2Vel(0) = dvxdt*delta_t;
         k2Pos(1) = dydt*delta_t;
         k2Vel(1) = dvydt*delta_t;
         k2Pos(2) = dydt*delta_t;
         k2Vel(2) = dvydt*delta_t;

         tempBody.position(0) = mainbody.position(0) + k2Pos(0)*0.5;
         tempBody.velocity(0) = mainbody.velocity(0) + k2Vel(0)*0.5;
         tempBody.position(1) = mainbody.position(1) + k2Pos(1)*0.5;
         tempBody.velocity(1) = mainbody.velocity(1) + k2Vel(1)*0.5;
         tempBody.position(2) = mainbody.position(2) + k2Pos(2)*0.5;
         tempBody.velocity(2) = mainbody.velocity(2) + k2Vel(2)*0.5;


         // k3
         force.zeros();
         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {

            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force(0) += F(0);
               force(1) += F(1);
               force(2) += F(2);
            }
         } //end for


         dxdt = tempBody.velocity(0);
         dvxdt = force(0);
         dydt = tempBody.velocity(1);
         dvydt = force(1);
         dzdt = tempBody.velocity(2);
         dvzdt = force(2);

         k3Pos(0) = dxdt*delta_t;
         k3Vel(0) = dvxdt*delta_t;
         k3Pos(1) = dydt*delta_t;
         k3Vel(1) = dvydt*delta_t;
         k3Pos(2) = dzdt*delta_t;
         k3Vel(2) = dvzdt*delta_t;

         tempBody.position(0) = mainbody.position(0) + k3Pos(0);
         tempBody.velocity(0) = mainbody.velocity(0) + k3Vel(0);
         tempBody.position(1) = mainbody.position(1) + k3Pos(1);
         tempBody.velocity(1) = mainbody.velocity(1) + k3Vel(1);
         tempBody.position(2) = mainbody.position(2) + k3Pos(2);
         tempBody.velocity(2) = mainbody.velocity(2) + k3Vel(2);



         // k4
         force.zeros();
         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {
            if (j != i)
            {
               arma::Col<double> F(3);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(tempBody, body,F);
               force(0) += F(0);
               force(1) += F(1);
               force(1) += F(1)
            }
         } //end for

         dxdt = tempBody.velocity(0);
         dvxdt = force(0);
         dydt = tempBody.velocity(1);
         dvydt = force(1);
         dzdt = tempBody.velocity(2);
         dvzdt = force(2);

         k4Pos(0) = dxdt*delta_t;
         k4Vel(0) = dvxdt*delta_t;
         k4Pos(1) = dydt*delta_t;
         k4Vel(1) = dvydt*delta_t;
         k4Pos(2) = dzdt*delta_t;
         k4Vel(2) = dvzdt*delta_t;

         tempBody.position(0) = mainbody.position(0) + 1.0/6.0 * (k1Pos(0) + 2.*k2Pos(0) + 2.*k3Pos(0) + k4Pos(0));
         tempBody.velocity(0) = mainbody.velocity(0) + 1.0/6.0 * (k1Vel(0) + 2.*k2Vel(0) + 2.*k3Vel(0) + k4Vel(0));
         tempBody.position(1) = mainbody.position(1) + 1.0/6.0 * (k1Pos(1) + 2.*k2Pos(1) + 2.*k3Pos(1) + k4Pos(1));
         tempBody.velocity(1) = mainbody.velocity(1) + 1.0/6.0 * (k1Vel(1) + 2.*k2Vel(1) + 2.*k3Vel(1) + k4Vel(1));
         tempBody.position(2) = mainbody.position(2) + 1.0/6.0 * (k1Pos(2) + 2.*k2Pos(2) + 2.*k3Pos(2) + k4Pos(2));
         tempBody.velocity(2) = mainbody.velocity(2) + 1.0/6.0 * (k1Vel(2) + 2.*k2Vel(2) + 2.*k3Vel(2) + k4Vel(2));



         mainbody.update(tempBody.position, tempBody.velocity);


         double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
         double potentialEnergy = 0.0;
         double angularMomentum = mysolarsystem.angularMomentum(mainbody);


         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {
            if (j != i)
            {
               Object body = mysolarsystem.objectlist[j];
               double energy = mysolarsystem.potentialEnergy(mainbody, body);
               potentialEnergy += energy;
            }

         } //end for



         fout <<"\t\t" << tempBody.position(0) << "\t\t" << tempBody.position(1) << "\t\t" << k*delta_t << "\t\t" << kineticEnergi << "\t\t" << potentialEnergy << "\t\t" << angularMomentum << "\t\t";
         if (i == mysolarsystem.numberOfObject - 1)
         {
            fout << "\n";
         }
      }
   }
   fout << "];" << "\n\n";
   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
   {
       Object thisobject = mysolarsystem.objectlist[j];
       fout << "plot (A(:," << 6*j+1 << "),A(:," << 6*j+2 <<"))" << "\n\n";
       fout << "legend('" << thisobject.name << "')" << "\n\n";
       fout << "hold on;" << "\n\n";

   }
//   fout << "figure()" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+3 << "),A(:," << 6*j+4 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << " " << "kinetic energy" << "')" << "\n\n";
//       fout << "hold on;" << "\n\n";
//   }

//   fout << "figure()" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+3 << "),A(:," << 6*j+5 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << " " << "potential energy" << "')" << "\n\n";
//       fout << "hold on;" << "\n\n";
//   }

//   fout << "figure()" << "\n\n";
//   for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
//   {
//       Object thisobject = mysolarsystem.objectlist[j];
//       fout << "plot (A(:," << 6*j+3 << "),A(:," << 6*j+6 <<"))" << "\n\n";
//       fout << "legend('" << thisobject.name << " " << "momentum" << "')" << "\n\n";
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
   fout << "describe = '[t x y]'" << "\n";
   fout << "A = [ " ;



   for (int k = 1 ; k<=n ; ++k)
   {

       for ( int i = 0 ; i < mysolarsystem.numberOfObject ; ++i)
       {
          Object &mainbody = mysolarsystem.objectlist[i];
          Object tempBody = mainbody;
          arma::Col<double> force(2);
          force.zeros();

          for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
          {

             if (j != i)
             {
                arma::Col<double> F(2);
                Object body = mysolarsystem.objectlist[j];
                mysolarsystem.force(tempBody, body,F);
                force(0) += F(0);
                force(1) += F(1);
             }

          } //end for

      tempBody.position(0) = mainbody.position(0)+mainbody.velocity(0)*delta_t +force(0)*delta_t*delta_t*0.5; //x_i+1
      tempBody.position(1) = mainbody.position(1)+mainbody.velocity(1)*delta_t +force(1)*delta_t*delta_t*0.5; //y_i+1

      arma::Col<double> newforce(2);
      newforce.zeros();

      for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
      {

         if (j != i)
         {
            arma::Col<double> F(2);
            Object body = mysolarsystem.objectlist[j];
            mysolarsystem.force(tempBody, body,F);
            newforce(0) += F(0);
            newforce(1) += F(1);
         }

      } //end for

      tempBody.velocity(0) = mainbody.velocity(0) + 0.5*(force(0)+newforce(0))*delta_t; //vx_i+1
      tempBody.velocity(1) = mainbody.velocity(1) + 0.5*(force(1)+newforce(1))*delta_t; //vy_i+1

      double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
      double potentialEnergy = 0.0;
      double angularMomentum = mysolarsystem.angularMomentum(mainbody);


      for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
      {
         if (j != i)
         {
            Object body = mysolarsystem.objectlist[j];
            double energy = mysolarsystem.potentialEnergy(mainbody, body);
            potentialEnergy += energy;
         }

      } //end for

      fout << k*delta_t << "\t\t" << tempBody.position(0) << "\t\t" << tempBody.position(1) << "\t\t"<< kineticEnergi << "\t\t" << potentialEnergy << "\t\t" << angularMomentum << "\t\t";
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
