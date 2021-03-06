#ifndef SYSTEM_H
#define SYSTEM_H

#include<armadillo>
#include<object.h>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>
#include<distance.h>

class System
{
public:
    System();


    void addObject(Object &newobject);

    void addSystem(System &newSubSystem);

//    void acceleration(Object &object1,
//               Object &object2,
//               arma::Col<double>& A);
    void acceleration(Object &mainObject,
                      int i);

    double kineticEnergi(Object movingObject);

    double potentialEnergy(Object movingObject,
                           Object otherObject);

    arma::Col<double> angularMomentum(Object movingObject);

    void solveRK4(double time,
                  double timestep);
    void solveVerlet(double time,
                  double timestep);

    int numberOfObject;
    std::vector<Object> objectlist;

    int numberOfSystems;
    std::vector<System> systemlist;

private:
   double PI;
};

#endif // SYSTEM_H
