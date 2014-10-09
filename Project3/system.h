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

class System
{
public:
    System();
    ~System();

    void addObject(arma::Col<double> pos,
                   arma::Col<double> vel,
                   double m,
                   char* name);

    void solve(const std::string filename);
    arma::Col<double> force(arma::Col<double> pos1,
                            arma::Col<double> pos2,
                            double GM);


    Object objectlist[20];
    int numberOfObject;


private:



};

#endif // SYSTEM_H
