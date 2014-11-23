#ifndef DISTANCE_H
#define DISTANCE_H
#include<armadillo>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
class Distance
{
public:
    double twoObjects(const arma::Col<double>&, const arma::Col<double>&) const;
};

#endif // DISTANCE_H
