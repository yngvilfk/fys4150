#ifndef OBJECT_H
#define OBJECT_H

#include<armadillo>


class Object
{
public:
    Object(double x,
           double y,
           double vx,
           double vy,
           const double objectMass,
           std::string objectName);

    Object(arma::Col<double> pos,
           arma::Col<double> vel,
           const double objectMass,
           std::string objectName);
    ~Object();
    void update(arma::Col<double> pos,
                arma::Col<double> vel);


    std::string name;
    double      mass;   //in solar masses

    arma::Col<double> position;   //in astronomical units (AU)
    arma::Col<double> velocity;   //in AU/year
};

#endif // OBJECT_H
