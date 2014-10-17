#include "object.h"

Object::Object(double x,
               double y,
               double vx,
               double vy,
               const double objectMass,
               std::string objectName)
    : name(objectName)
    , mass(objectMass)
{
   arma::Col<double> createColDim(2);

   position = createColDim;
   velocity = createColDim;
   position(0) = x;
   position(1) = y;
   velocity(0) = vx;
   velocity(1) = vy;
}


Object::Object(arma::Col<double> pos,
               arma::Col<double> vel,
               const double objectMass,
               std::string objectName)
    : name(objectName)
    , mass(objectMass)
    , position(pos)
    , velocity(vel)
{
}

Object::~Object()
{
}

void
Object::update(arma::Col<double> pos,
               arma::Col<double> vel)
{
   position = pos;
   velocity = vel;
}
