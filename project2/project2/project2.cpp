#include <string>
#include <iostream> //i standardbiblioteket
#include <vector> //ligger i standardbiblioteket
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <fstream>


using namespace arma;
using namespace std;



namespace
{
   namespace Local
   {
      void jacobi(arma::Mat<double>& A,
                  arma::Mat<double>& B,
                  const int&         n);

      double maxOffDiag(const arma::Mat<double>& A,
                      int&                     k,
                      int&                     l,
                      const int&               n);

      void rotate(arma::Mat<double>& A,
                  arma::Mat<double>& B,
                  const int&         k,
                  const int&         l,
                  const int&         n);
   }
}



int main()
{
   double h, rho, rhoMaks, rhoMin,  e, V, E, hbar, kons, m, d;
   int n, nsteps;
   const int nsteps = 100;
   const int n = nsteps-1;
   const double hbar=1;
   const double alpha = pow((hbar*hbar/m*kons),(0.25));
   const double lambda = (2.*m*alpha*alpha*E)/(hbar*hbar);
   rhoMin = 0.;
   rhoMaks = 100.;
   h = (rhoMaks-rhoMin)/static_cast<double>(nsteps);
   e = -1./(h*h);

   arma::Mat<double> A(n,n), B(n,n);
   A.zeros();
   //fill matrix A
   //first row:
   double rho1 = rhoMin + h;
   V = rho1*rho1;
   d = 2./(h*h) + V;
   A(1,1) = d;
   A(1,2) = e;
   //row n (n_steps-1)
   double rho_n = rhoMin + static_cast<double>(n-1)*h;
   V = rho_n*rho_n;
   d = 2./(h*h) + V;
   A(n-1,n-1) = d;
   A(n-1,n-2) = e;
   //row 2 to n-1 :
   for (int i = 1 ; i<n-1 ; ++i)
   {
      rho = rhoMin + static_cast<double>(i)*h;
      V = rho*rho;
      d = 2./(h*h) + V;
      A(i,i) = d;
      A(i,i+1) = e;
      A(i,i-1) = e;
      i+=1;

   }
   Local::jacobi(A,B,n);
   cout << A << endl;
}



void
Local::jacobi(arma::Mat<double>& A,
              arma::Mat<double>& B,
              const int& n)
{
   //setting up eigenvector-matrix:
   B.eye();
   int k, l;
   double epsilon = 1.0e-8;
   double aMax = maxOffDiag(A, k, l, n);
   int maxIterations = n*n*n;
   int iterations = 0;
   while (fabs(aMax) > epsilon && iterations < maxIterations)
   {
      rotate(A, B, k, l, n);
      A = B;
      aMax = maxOffDiag(A, k, l, n);
      ++ iterations;
   }
   return;
}



double
Local::maxOffDiag(const arma::Mat<double>& A,
                 int& k,
                 int& l,
                 const int& n)
{
   double max = 0.0;
   for (int i = 0 ; i<n ; ++i)
   {
      for (int j = i+1 ; j<n ; ++j)
      {
         if (fabs(A(i,j)) > max)
         {
            max = fabs(A(i,j));
            k = i;
            l = j;
         }
      }
   }
   return max;
}



void
Local::rotate(mat& A,
            mat& B,
            const int& k,
            const int& l,
            const int& n)
{
double c, s;
   if(A(k,l) != 0.0)
   {
      double t;
      const double tau = (A(l,l)-A(k,k))/(2*A(k,l));
      if (tau < 0.0)
      {
            t = -tau - sqrt(1. + tau*tau);
      }
      else
      {
         t = -tau + sqrt(1. + tau*tau);
      }
      c = 1.0/(sqrt(1.0+ t*t));
      s = t*c;
   }
   else
   {
      c = 1.0;
      s = 0.0;
   }
   mat S;
   S.eye(n,n);
   S(k,k) = c;
   S(l,l) = c;
   S(k,l) = -s;
   S(l,k) = s;
   B = S.t()*A*S;
}


