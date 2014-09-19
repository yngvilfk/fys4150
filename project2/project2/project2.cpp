#include <cmath>
#include <armadillo>
#include <math.h>



int main()
{

}



implementJacobi(double ** A, double ** R, int n)
{
//function to implement the Jacobis rotation algoritm
    int k, l;
    double epsilon = 1.0e-8;
    double maks_offdiag = maxOffDiag(A, &k, &l, n);
    int maxIterations = n * n * n;
    int iterations = 0;
    mat B(n,n);
    mat S(n,n);
    while ( fabs(maks_offdiag) > epsilon && iterations < maxIterations)
    {
        S.fill(0.0);
        double tau = (A(l,l)-A(k,k))/2*A(k,l); //cot(2*theta)
        double t = -tau + sqrt(1.0 + tau*tau);    // tan(theta)
        //or: double t = -tau - sqrt(1. + tau*tau);
        double c = 1.0/(sqrt(1.0+ t*t));          //cos(theta)
        double s = t*c;                           //sin(theta)
        for (int i = 0 ; i<n ; ++i)
        {
            S(i,i) = 1.0;
        }
        S(k,k) = c;
        S(l,l) = c;
        S(k,l) = -s;
        S(l,k) = s;

        double maks_offdiag = maxOffDiag(A, &k, &l, n);
    }

}




double maxOffDiag(double ** A, int * k, int * l, int n )
{
//find max, non-diagonal matrix element in a symmetric matrix and save line and colomn-number
    double max = 0.0;
            for ( int i = 0 ; i<n ; ++i)
            {
                for ( int j = i+1 ; j<n ; ++j)
                {
                    if( fabs(A(i,j)) > max)
                    {
                        max = fabs(A(i,j));
                        *k = i;
                        *l = j;
                    }
                }
            }
    return max
}

int S
