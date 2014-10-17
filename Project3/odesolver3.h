#ifndef ODESOLVER3_H
#define ODESOLVER3_H

#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<distance.h>
#include<iomanip>
#include<system.h>
#include<object.h>

class Odesolver3
{
public:
    Odesolver3(System &mysystem);
    ~Odesolver3();

    void rk4(double time,
             int nSteps);
    void verlet(double time,
                int nSteps);

    System mysolarsystem;

 protected:


 private:
    double delta_t;
    int n;




};

#endif // ODESOLVER3_H
