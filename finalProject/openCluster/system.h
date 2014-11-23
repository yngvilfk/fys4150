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
    void removeObject(int i);

    void setEpsilon(const double& epsilon) { epsilon_ = epsilon ;}

    void acceleration(Object &mainObject,
                      int i);

    double kineticEnergi() const;

    double potentialEnergy() const;

    double totalEnergy() const { return kineticEnergi()+potentialEnergy(); }

    double boundKineticEnergi(double limit) const;

    double boundPotentialEnergy(double limit) const;

    double boundTotalEnergy(double limit) const { return boundKineticEnergi(limit)+boundPotentialEnergy(limit); }

    int numberOfObject;
    std::vector<Object> objectlist;

    int numberOfSystems;
    std::vector<System> systemlist;

private:
   double PI;
   double epsilon_;
};

#endif // SYSTEM_H
