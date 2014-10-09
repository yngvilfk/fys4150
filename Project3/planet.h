#ifndef PLANET_H
#define PLANET_H

#include<armadillo>

class Planet
{
public:
    Planet();
    ~Planet();

    void earth(double x, double y, double vx, double vy);


private:
    double sunx;
    double suny;
    double earthx;
    double earthy;
    double earthvx;
    double earthvy;
    double massEarth;
    double jupiterx;
    double jupitery;
    double jupitervx;
    double jupitervy;
    double massJupiter;


};

#endif // PLANET_H
