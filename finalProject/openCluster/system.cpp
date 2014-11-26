#include "system.h"

System::System()
    : numberOfObject(0)
    , numberOfSystems(0)
    , PI(4.*std::atan(1.0))
    , epsilon_(0.0)
 {
 }


 void
 System::addObject(Object &newobject)
 {
     objectlist.push_back(newobject);
     numberOfObject = static_cast<int>(objectlist.size());
 }


 void
 System::removeObject(int i)
 {
    objectlist.erase(objectlist.begin()+i);
    numberOfObject = static_cast<int>(objectlist.size());
 }


 void
 System::acceleration(Object &mainObject,
                      int i,
                      double time,
                      std::string length)
 {
    arma::Col<double> tempAcceleration(3);
    arma::Col<double> acceleration(3);
    acceleration.zeros();
    for (int j = 0; j < numberOfObject; ++j)
    {
       if (j!=i)
       {
          const Object tempObject = objectlist[j];
          Distance d;
          arma::Col<double> position1, position2;
          position1 = mainObject.getPosition();
          position2 = tempObject.getPosition()+ time*tempObject.getVelocity();
          double R = d.twoObjects(position1, position2);
          tempAcceleration = -4*PI*PI*tempObject.getMass()*(mainObject.getPosition()-tempObject.getPosition())/
                               (R*(R*R+epsilon_*epsilon_));// /6.3241e4);//[ly/yr]
          if (length == "ly")
          {
             tempAcceleration = tempAcceleration/6.3241e4;
          }

          acceleration += tempAcceleration;
       } //end if
    } //end for
    mainObject.setAcceleration(acceleration);
 }

 double
 System::maxTimestep(int i)
 {
    Distance d;
    const Object mainObject = objectlist[i];
    int k;
    if (i == 0)
    {
       k=1;
    }
    else
    {
       k=0;
    }
    Object tempObject = objectlist[k];
    arma::Col<double> position1(3), position2(3);
    arma::Col<double> mainVel = mainObject.getVelocity();

    position1 = mainObject.getPosition();
    position2 = tempObject.getPosition();//+ time*tempObject.getVelocity();

    double R = d.twoObjects(position1, position2);
    arma::Col<double> tempVel = tempObject.getVelocity()-mainVel;
    double V = std::sqrt(arma::dot(tempVel,tempVel));
    double maxTime;
    if (V>1e-5)
    {
       maxTime = R/V;
    }
    else
    {
       maxTime = R/1e-5;
    }
    double temptime;

    for (int j = 1; j < numberOfObject; ++j)
    {
       if (j!=i)
       {
          tempObject = objectlist[j];
          position2 = tempObject.getPosition();//+ time*tempObject.getVelocity();
          R = d.twoObjects(position1, position2);
          tempVel = tempObject.getVelocity()-mainVel;
          V = std::sqrt(arma::dot(tempVel,tempVel));
          temptime = R/V;
          if (temptime<maxTime && V<1e-5)
          {
             maxTime = temptime;
          }
       } //end if
    } //end for
    return maxTime;
 }


 double
 System::kineticEnergi() const
 {
    double energy = 0;
    for (int i = 0 ; i < numberOfObject ; ++i)
    {
       const Object movingObject = objectlist[i];
       energy += 0.5 * movingObject.getMass() *
             arma::dot(movingObject.getVelocity(), movingObject.getVelocity()) ; //* solarmass * velocity * velocity ;
    }

    return energy;
 }


 double
 System::potentialEnergy() const
 {
    Distance d;
    double energy = 0.0;
    for (int i = 0 ; i < numberOfObject ; ++i)
    {
       for (int j = 0 ; j < numberOfObject ; ++j)
       {
          if (i!=j)
          {
             const Object movingObject = objectlist[i];
             const Object otherObject  = objectlist[j];
             const double R = d.twoObjects(movingObject.getPosition(), otherObject.getPosition());
             energy +=  4 * PI * PI * otherObject.getMass() *
                   movingObject.getMass() /R;
          }
       }
    }

    return energy;
 }


 double
 System::boundKineticEnergi(double limit) const
 {
    double energy = 0;
    arma::Col<double> center(3), position(3);
    center.zeros();
    Distance d;
    for (int i = 0 ; i < numberOfObject ; ++i)
    {
       Object movingObject = objectlist[i];
       position = movingObject.getPosition();
       if (d.twoObjects(center, position)<limit)
       {
          energy += 0.5 * movingObject.getMass() *
                arma::dot(movingObject.getVelocity(), movingObject.getVelocity()) ; //* solarmass * velocity * velocity ;
       }
    }

    return energy;
 }


 double
 System::boundPotentialEnergy(double limit) const
 {
    Distance d;
    double energy = 0.0;
    arma::Col<double> position1, position2, center(3);
    for (int i = 0 ; i < numberOfObject ; ++i)
    {
       for (int j = 0 ; j < numberOfObject ; ++j)
       {
          if (i!=j)
          {
             Object movingObject = objectlist[i];
             Object otherObject = objectlist[j];
             position1 = movingObject.getPosition();
             position2 = otherObject.getPosition();
             if (d.twoObjects(center, position1)<limit)
             {
                double R = d.twoObjects(position1, position2);
                energy +=  4 * PI * PI * otherObject.getMass() *
                      movingObject.getMass() /R;
             }
          }
       }
    }
    return energy;
 }
