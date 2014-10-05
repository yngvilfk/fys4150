#include "distance.h"

Distance::Distance()

{
}

Distance::~Distance()
{
}

double Distance::twoObjects(const arma::Col<double>& object1,
                            const arma::Col<double>& object2)
{
    /*function takes in the coordinates of two objects, object 1 and object 2:
    object 1 have the coordinates (x,y): (object1(0), object1(1))
    object 2 have the coordinates (x,y): (object2(0), object2(1))
     */
    const double x1 = object1(0);
    const double y1 = object1(1);
    const double x2 = object2(0);
    const double y2 = object2(1);
    const double R = std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return R;
}
