#include "odesolver2.h"



OdeSolver2::OdeSolver2(System &mysystem)
   : mysolarsystem(mysystem)
{

}


OdeSolver2::~OdeSolver2()
{
}


//OdeSolver2::rk4( std::string& fileName)
void
OdeSolver2::rk4(double time,
                int nSteps)
{
   arma::Col<double> k1x(2), k2x(2), k3x(2), k4x(2), k1y(2), k2y(2), k3y(2), k4y(2);
   n = nSteps;
   delta_t = time/(n+1);


   //std::ofstream fout(fileName.c_str());
   std::ofstream fout("test.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(16);
   fout << "describe = '[objectNmbr x y vx vy objectNmbr x y vx vy etc]'" << "\n";
   fout << "numberOfObjects = " << mysolarsystem.numberOfObject << "\n";
   fout << "A = [ " ;


   for ( int k = 0 ; k < n ; ++k )
   {

      for (int i = 1 ; i < mysolarsystem.numberOfObject ; ++i)
      {

         Object &mainbody = mysolarsystem.objectlist[i];
         arma::Col<double> force(2);
         force.zeros();

         double kineticEnergi = mysolarsystem.kineticEnergi(mainbody);
         double potentialEnergy = 0.0;
         double angularMomentum = mysolarsystem.angularMomentum(mainbody);


         for ( int j = 0 ; j < mysolarsystem.numberOfObject ; ++j)
         {
            if (j != i)
            {
               arma::Col<double> F(2);
               Object body = mysolarsystem.objectlist[j];
               mysolarsystem.force(mainbody, body,F);
               force(0) += F(0);
               force(1) += F(1);
               double energy = mysolarsystem.potentialEnergy(mainbody, body);
               potentialEnergy += energy;
            }

         } //end for





         arma::Col<double> position_k(2), velocity_k(2);


         position_k(0) = mainbody.position(0);
         velocity_k(0) = mainbody.velocity(0);
         position_k(1)= mainbody.position(1);
         velocity_k(1) = mainbody.velocity(1);



         // k1
         double dxdt = velocity_k(0);
         double dvxdt = force(0);

         double dydt = velocity_k(1);
         double dvydt = force(1);


         k1x(0) = dxdt*delta_t;
         k1x(1) = dvxdt*delta_t;

         k1y(0) = dydt*delta_t;
         k1y(1) = dvydt*delta_t;


         position_k(0) = mainbody.position(0) + k1x(0)*0.5;
         velocity_k(0) = mainbody.velocity(0) + k1x(1)*0.5;



         position_k(1) = mainbody.position(1) + k1y(0)*0.5;
         velocity_k(1) = mainbody.velocity(1) + k1y(1)*0.5;


         // k2
         dxdt = velocity_k(0);
         dvxdt = force(0);
         dydt = velocity_k(1);
         dvydt = force(1);

         k2x(0) = dxdt*delta_t;
         k2x(1) = dvxdt*delta_t;
         k2y(0) = dydt*delta_t;
         k2y(1) = dvydt*delta_t;

         position_k(0) = mainbody.position(0) + k2x(0)*0.5;
         velocity_k(0) = mainbody.velocity(0) + k2x(1)*0.5;
         position_k(1) = mainbody.position(1) + k2y(0)*0.5;
         velocity_k(1) = mainbody.velocity(1) + k2y(1)*0.5;



         // k3
         dxdt = velocity_k(0);
         dvxdt = force(0);
         dydt = velocity_k(1);
         dvydt = force(1);

         k3x(0) = dxdt*delta_t;
         k3x(1) = dvxdt*delta_t;
         k3y(0) = dydt*delta_t;
         k3y(1) = dvydt*delta_t;

         position_k(0) = mainbody.position(0) + k3x(0)*0.5;
         velocity_k(0) = mainbody.velocity(0) + k3x(1)*0.5;
         position_k(1) = mainbody.position(1) + k3y(0)*0.5;
         velocity_k(1) = mainbody.velocity(1) + k3y(1)*0.5;




         // k4
         dxdt = velocity_k(0);
         dvxdt = force(0);
         dydt = velocity_k(1);
         dvydt = force(1);

         k4x(0) = dxdt*delta_t;
         k4x(1) = dvxdt*delta_t;
         k4y(0) = dydt*delta_t;
         k4y(1) = dvydt*delta_t;

         position_k(0) = mainbody.position(0) + 1.0/6.0 * (k1x(0) + 2.*k2x(0) + 2.*k3x(0) + k4x(0));
         velocity_k(0) = mainbody.velocity(0) + 1.0/6.0 * (k1x(1) + 2.*k2x(1) + 2.*k3x(1) + k4x(1));
         position_k(1) = mainbody.position(1) + 1.0/6.0 * (k1y(0) + 2.*k2y(0) + 2.*k3y(0) + k4y(0));
         velocity_k(1) = mainbody.velocity(1) + 1.0/6.0 * (k1y(1) + 2.*k2y(1) + 2.*k3y(1) + k4y(1));





         mainbody.update(position_k, velocity_k);

         fout << i << "\t\t" << position_k(0) << "\t\t" << position_k(1) << "\t\t" << k*delta_t << "\t\t" << kineticEnergi << "\t\t" << potentialEnergy << "\t\t" << angularMomentum << "\t\t";
         if (i == mysolarsystem.numberOfObject - 1)
         {
            fout << "\n";
         }


      }


   }
   fout << "];" << "\n\n";
   fout << "plot (A(:,2),A(:,3))" << "\n\n";
//   fout << "figure();" << "\n\n";
//   fout << "plot (A(:,4),A(:,5), '-c')" << "\n\n";
//   fout << "hold on;" << "\n\n";
   //fout << "plot (A(:,4),A(:,6), '-g')" << "\n\n";
   //fout << "hold on;" << "\n\n";
   //fout << "plot (A(:,4),A(:,7), '-r')" << "\n\n";
//   fout << "legend ('sun','earth','jupiter')" << "\n\n";
   fout.close();



}



void
OdeSolver2::verlet()
{
   double PI = 4*std::atan(1.0);
   arma::Col<double> y(2), x(2), vy(2), vx(2);
   //startpos
   x(0) = x0;
   y(0) = y0;
   vx(0) = v0x;
   vy(0) = v0y;


   arma::Col<double> earth(2), sun(2);
   sun(0) = 0.0;
   sun(1) = 0.0;
   earth(0) = x(0);
   earth(1) = y(0);


   std::ofstream fout("verlet1.m");
   //fout.setf(std::ios::scientific);
   fout.precision(16);
   //fout.width(8);
   fout << "describe = '[t x y]'" << "\n";
   fout << "A = [ " ;
   fout << 0.0 << "\t\t" << x(0) << "\t\t" << y(0) << "\n";



   for (int i = 1 ; i<=n ; ++i)
   {
       earth(0) = x(0);
       earth(1) = y(0);
       Distance d;
      const double R = d.twoObjects(earth,sun);
      std::cout << R << std::endl;

      x(1) = x(0)+vx(0)*delta_t - 0.5*delta_t*delta_t*4*PI*PI*x(0)/(R*R*R); //x_i+1
      y(1) = y(0)+vy(0)*delta_t - 0.5*delta_t*delta_t*4*PI*PI*y(0)/(R*R*R); //y_i+1
      vx(1) = vx(0) - 0.5*(4*PI*PI*x(0)/(R*R*R)+4*PI*PI*x(1)/(R*R*R))*delta_t; //vx_i+1
      vy(1) = vy(0) - 0.5*(4*PI*PI*y(0)/(R*R*R)+4*PI*PI*y(1)/(R*R*R))*delta_t; //vy_i+1

      fout << i*delta_t << "\t\t" << x(1) << "\t\t" << y(1) << "\n";
      x(0) = x(1);
      vx(0) = vx(1);
      y(0) = y(1);
      vy(0) = vy(1);


   }
}
