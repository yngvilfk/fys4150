#include "distance.h"

double Distance::twoObjects(const arma::Col<double>& object1,
                            const arma::Col<double>& object2) const
{
    /*function takes in the coordinates of two objects, object 1 and object 2:
    object 1 have the coordinates (x,y): (object1(0), object1(1))
    object 2 have the coordinates (x,y): (object2(0), object2(1))
     */
   const arma::Col<double> tempPos = object1-object2;
    const double R = std::sqrt(arma::dot(tempPos, tempPos));
    return R;
}

