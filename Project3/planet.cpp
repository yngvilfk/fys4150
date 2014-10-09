#include "planet.h"

Planet::Planet()
{
   suny = 0.0;
   sunx = 0.0;
}

Planet::~Planet()
{
}

void
Planet::earth(double x, double y, double vx, double vy)
{
   earthx = x;
   earthy = y;
   earthvx = vx;
   earthvy = vy;
}

//double
//Planet::distance(int object1, int object2)
//{
//   double R = 0.0;
//   return R;
//}

