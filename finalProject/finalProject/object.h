#ifndef OBJECT_H
#define OBJECT_H

#include<armadillo>
#include<fstream>
#include<iostream>

class Object
{
public:
    Object(double x,
           double y,
           double z,
           double vx,
           double vy,
           double vz,
           const double objectMass,
           std::string objectName);

    Object(arma::Col<double> pos,
           arma::Col<double> vel,
           const double objectMass,
           std::string objectName);

    void update(arma::Col<double> pos,
                arma::Col<double> vel);
    void                     setPosition(const arma::Col<double>& pos) { position_ = pos;}
    const arma::Col<double>& getPosition() const {return position_;}
    void                     setVelocity(const arma::Col<double>& vel) { velocity_ = vel; }
    const arma::Col<double>& getVelocity() const {return velocity_;}
    void                     setAcceleration(const arma::Col<double>& acceleration) { acceleration_ = acceleration ;}
    const arma::Col<double>& getAcceleration() const {return acceleration_;}
    const double& getMass() const {return mass_;}
    const std::string& getName() const {return name_;}
    double maxTimestep();
    void newFile();
    void addToFile();
    void closeFile();

protected:

private:
    std::string name_;
    double      mass_;   //in solar masse

    arma::Col<double> position_;   //in astronomical units (AU)
    arma::Col<double> velocity_;   //in AU/year
    arma::Col<double> acceleration_;
    double timestep_;

};

#endif // OBJECT_H

