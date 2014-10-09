#include "object.h"

Object::Object()
{
}

Object::~Object()
{
}

Object::Object(arma::Col<double> pos,
               arma::Col<double> vel,
               double objectMass,
               char* objectName)
    : position(pos)
    , velocity(vel)
    , name(objectName)
    , mass(objectMass)
{
}
