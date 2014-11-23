#include "object.h"

Object::Object(double x,
               double y,
               double z,
               double vx,
               double vy,
               double vz,
               const double objectMass)
    : mass_(objectMass)
{
   arma::Col<double> createColDim(3);

   position_ = createColDim;
   velocity_ = createColDim;
   position_(0) = x;
   position_(1) = y;
   position_(2) = z;
   velocity_(0) = vx;
   velocity_(1) = vy;
   velocity_(2) = vz;
   acceleration_ = arma::zeros<arma::vec>(3);
}


Object::Object(arma::Col<double> pos,
               arma::Col<double> vel,
               const double objectMass)
    : mass_(objectMass)
    , position_(pos)
    , velocity_(vel)
{
    acceleration_ = arma::zeros<arma::vec>(3);
}

double
Object::maxTimestep()
{
   timestep_ = 1./std::sqrt(arma::dot(acceleration_,acceleration_));
   return timestep_;
}

void
Object::newFile()
{
   std::string filename = "object_.m";
   std::ofstream fout(filename.c_str());
   fout << "A = [";
   fout.close();
}



void
Object::addToFile()
{
   std::string filename = "object_.m";
   std::ofstream fout(filename.c_str(), std::ios::app);
   fout << position_(0) << "\t\t" << position_(1) << "\t\t" << position_(2) << "\n";
   fout.close();
}



void
Object::closeFile()
{
   std::string filename = "object_.m";
   std::ofstream fout(filename.c_str(), std::ios::app);
   fout << "] \n";
//   fout << "plot3(A(:,1),A(:,2),A(:,3))" ;
   fout << "plot(A(:,1),A(:,2))" ;
   fout.close();
}
