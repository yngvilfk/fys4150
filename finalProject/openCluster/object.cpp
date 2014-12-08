#include "object.h"


Object::Object(arma::Col<double> pos,
               arma::Col<double> vel,
               const double objectMass)
    : mass_(objectMass)
    , position_(pos)
    , velocity_(vel)
{
    acceleration_ = arma::zeros<arma::vec>(3);
}

Object::Object(arma::Col<double> pos,
               arma::Col<double> vel,
               const double objectMass,
               std::string objectName)
    : name_(objectName)
    , mass_(objectMass)
    , position_(pos)
    , velocity_(vel)
{
    acceleration_ = arma::zeros<arma::vec>(3);
}

double
Object::maxTimestep()
{
   /*returns a suggestion for maximum timestep for the given
    * object*/
   timestep_ = 1./std::sqrt(arma::dot(acceleration_,acceleration_));
   return timestep_;
}

void
Object::newFile(std::string addition)
{
   /*This function creates a file
    * 'object_<name><addittion>.m'
    * where the position of the
    * object can be stored to.*/
   std::string filename = "object_";
   filename.append(name_);
   filename.append(addition);
   filename.append(".m");
   std::ofstream fout(filename.c_str());
   fout << "A = [";
   fout.close();
}



void
Object::addToFile(std::string addition)
{
   /*addToFile adds the present position
    * to the file 'object_<name><addition>.m*/
   std::string filename = "object_";
   filename.append(name_);
   filename.append(addition);
   filename.append(".m");
   std::ofstream fout(filename.c_str(), std::ios::app);
   fout << position_(0) << "\t\t" << position_(1)
        << "\t\t" << position_(2) << "\n";
   fout.close();
}



void
Object::closeFile(std::string addition)
{
   /*closeFile close the file
    * 'object_<name><addition>.m' after
    * all positions are saved*/
   std::string filename = "object_";
   filename.append(name_);
   filename.append(addition);
   filename.append(".m");
   std::ofstream fout(filename.c_str(), std::ios::app);
   fout << "] \n";
   fout.close();
}
