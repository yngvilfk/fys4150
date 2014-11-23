#include "system.h"

System::System()
    : numberOfObject(0)
    , numberOfSystems(0)
 {
    PI = 4*std::atan(1.0);
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



// void
// System::acceleration(Object &object1,
//               Object &object2,
//               arma::Col<double>& A)
// {
//     //pos1: position object 1, pos 2: position object 2,
//     Distance d;
//     arma::Col<double> position1, position2;
//     position1 = object1.getPosition();
//     position2 = object2.getPosition();
////     double R = d.twoObjects(object1.getPosition(), );
//     double R = d.twoObjects(position1, position2);
//     A = -4*PI*PI*object2.getMass()*(object1.getPosition()-object2.getPosition())/(R*R*R);
// }

 void
 System::acceleration(Object &mainObject,
                      int i)
 {
    arma::Col<double> tempAcceleration(3);
    arma::Col<double> acceleration(3);
    acceleration.zeros();
    for (int j = 0; j < numberOfObject; ++j)
    {
       if (j!=i)
       {
          Object tempObject = objectlist[j];
          Distance d;
          arma::Col<double> position1, position2;
          position1 = mainObject.getPosition();
          position2 = tempObject.getPosition();
//          double R = d.twoObjects(position1, position2);
          double R = Distance::twoObjects(position1, position2);
          tempAcceleration = -4*PI*PI*tempObject.getMass()*(mainObject.getPosition()-tempObject.getPosition())/(R*R*R);//
          acceleration += tempAcceleration;
       } //end if
    } //end for
    mainObject.setAcceleration(acceleration);
 }


 double
 System::kineticEnergi(Object movingObject)
 {
    double energy = 0.5 * movingObject.getMass() *arma::dot(movingObject.getVelocity(), movingObject.getVelocity()) ; //* solarmass * velocity * velocity ;
    return energy;
 }


 double
 System::potentialEnergy(Object movingObject,
                         Object otherObject)
 {
    Distance d;
    arma::Col<double> position1, position2;
    position1 = movingObject.getPosition();
    position2 = otherObject.getPosition();
    double R = d.twoObjects(position1, position2);
//    double R = d.twoObjects(movingObject.getPosition(), otherObject.getPosition());
    double energy =  4 * PI * PI * otherObject.getMass() * movingObject.getMass() /R;//* solarmass  * velocitySI * velocitySI/ ( AU);
    return energy;
 }


 arma::Col<double>
 System::angularMomentum(Object movingObject)
 {

    arma::Col<double> momentum = arma::cross(movingObject.getPosition(), movingObject.getMass()*movingObject.getVelocity());// * AU * velocitySI;
    return momentum;
 }


 void
 System::solveRK4(double time,
               double timestep)
 {
//    for (int i < numberOfObject ; ++i)
//    {
//       Object &mainbody = objectlist[i];
//       acceleration(mainbody, i);
//       SolveStep::rk4(delta_t, i, mainbody, mainbody);
//       acceleration(mainbody, i);
//    }
 }

 void
 System::solveVerlet(double time,
               double timestep)
 {

 }
