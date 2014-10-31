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
    ~System();

    void addObject(Object &newobject);

    void addSystem(System &newSubSystem);

    void force(Object &object1,
               Object &object2,
               arma::Col<double>& F);

    double kineticEnergi(Object movingObject);

    double potentialEnergy(Object movingObject,
                           Object otherObject);

    double maxTimestep(Object movingObject);

    double maxTimestep(System movingSystem);


    arma::Col<double> angularMomentum(Object movingObject);


    int numberOfObject;
    std::vector<Object> objectlist;

    int numberOfSystems;
    std::vector<System> systemlist;

private:
   double PI;
};

#endif // SYSTEM_H
