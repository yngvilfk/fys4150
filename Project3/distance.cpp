#include "distance.h"

Distance::Distance()

{
}

Distance::~Distance()
{
}

double Distance::twoObjects(arma::Col<double> object1,
                            arma::Col<double> object2)
{
    /*function takes in the coordinates of two objects, object 1 and object 2:
    object 1 have the coordinates (x,y): (object1(0), object1(1))
    object 2 have the coordinates (x,y): (object2(0), object2(1))
     */
    double x1 = object1(0);
    double y1 = object1(1);
    double x2 = object2(0);
    double y2 = object2(1);
    double R = std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return R;
}
