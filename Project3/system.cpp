#include "system.h"


System::System()
   : numberOfObject(0)
{
}

System::~System()
{
}


void
System::addObject(Object &newobject)
{
    PI = 4*std::atan(1.0);
    objectlist.push_back(newobject);
    numberOfObject = objectlist.size();
}



void
System::force(Object &object1,
              Object &object2,
              arma::Col<double>& F)
{
    //pos1: position object 1, pos 2: position object 2, GM: gravitational constant times mass of object 2
    double PI = 4*std::atan(1.0);
    Distance d;
    double R = d.twoObjects(object1.position, object2.position);
    F(0) = -4*PI*PI*object2.mass*(object1.position(0)-object2.position(0))/(R*R*R);
    F(1) = -4*PI*PI*object2.mass*(object1.position(1)-object2.position(1))/(R*R*R);
}

double
System::kineticEnergi(Object movingObject)
{
   double energy = 0.5 * movingObject.mass * movingObject.mass * movingObject.mass; //* solarmass * velocitySI * velocitySI ;
   return energy;
}


double
System::potentialEnergy(Object movingObject,
                        Object otherObject)
{
   Distance d;
   double R = d.twoObjects(movingObject.position, otherObject.position);
   double energy = - 4 * PI * PI * otherObject.mass * movingObject.mass /(R*R);//* solarmass  * velocitySI * velocitySI/ ( AU);
   return energy;
}


double
System::angularMomentum(Object movingObject)
{
   double momentum = (movingObject.position(0)*movingObject.velocity(1)-movingObject.position(1)*movingObject.velocity(0));// * AU * velocitySI;
   return momentum;
}
