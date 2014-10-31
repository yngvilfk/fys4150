#include "system.h"

System::System()
    : numberOfObject(0)
    , numberOfSystems(0)
 {
    PI = 4*std::atan(1.0);
 }

 System::~System()
 {
 }

//System is build so that every system needs a center object.
//the fist added object is the center object, object and systems
//are orbiting around
 void
 System::addObject(Object &newobject)
 {
     objectlist.push_back(newobject);
     numberOfObject = objectlist.size();
 }

 void
 System::addSystem(System &newSubSystem)
 {
     systemlist.push_back(newSubSystem);
     numberOfSystems = systemlist.size();
 }



 void
 System::force(Object &object1,
               Object &object2,
               arma::Col<double>& F)
 {
     //pos1: position object 1, pos 2: position object 2,
     Distance d;
     double R = d.twoObjects(object1.position, object2.position);
//     F(0) = -4*PI*PI*object2.mass*(object1.position(0)-object2.position(0))/(R*R*R);
//     F(1) = -4*PI*PI*object2.mass*(object1.position(1)-object2.position(1))/(R*R*R);
//     F(2) = -4*PI*PI*object2.mass*(object1.position(2)-object2.position(2))/(R*R*R);
       F = -4*PI*PI*object2.mass*(object1.position-object2.position)/(R*R*R);
 }


 double
 System::kineticEnergi(Object movingObject)
 {
    double energy = 0.5 * movingObject.mass *arma::dot(movingObject.velocity, movingObject.velocity) ; //* solarmass * velocity * velocity ;
    return energy;
 }


 double
 System::potentialEnergy(Object movingObject,
                         Object otherObject)
 {
    Distance d;
    double R = d.twoObjects(movingObject.position, otherObject.position);
    double energy =  4 * PI * PI * otherObject.mass * movingObject.mass /R;//* solarmass  * velocitySI * velocitySI/ ( AU);
    return energy;
 }


 arma::Col<double>
 System::angularMomentum(Object movingObject)
 {

    arma::Col<double> momentum = arma::cross(movingObject.position, movingObject.mass*movingObject.velocity);// * AU * velocitySI;
    return momentum;
 }
 double
 System::maxTimestep(Object movingObject)
 {
     Object centerObject = objectlist[0];
     Distance d;
     double R = d.twoObjects(movingObject.position, centerObject.position);
     arma::Col<double> v = movingObject.velocity-centerObject.velocity;
     double V = std::sqrt(arma::dot(v,v));
     double dt = R/V;
     return dt;
 }

 double
 System::maxTimestep(System movingSystem)
 {
     Object centerObject = objectlist[0];
     Object movingObject = movingSystem.objectlist[0];
     Distance d;
     double R = d.twoObjects(movingObject.position, centerObject.position);
     arma::Col<double> v = movingObject.velocity-centerObject.velocity;
     double V = std::sqrt(arma::dot(v,v));
     double dt = R/V;
     return dt;
 }
