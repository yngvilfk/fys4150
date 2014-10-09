#ifndef OBJECT_H
#define OBJECT_H

#include<armadillo>


class Object
{
public:
    Object();
    Object(arma::Col<double> pos,
           arma::Col<double> vel,
           double objectMass,
           char* objectName);
    ~Object();

private:
    char* name;
    double mass;
    arma::Col<double> position;
    arma::Col<double> velocity;


};

#endif // OBJECT_H
