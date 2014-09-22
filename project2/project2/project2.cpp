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



int main()
{
   double h, rho, rhoMaks, rhoMin, alpha, e, V, lambda, E, hbar, kons, m, d;
   int n, nsteps;

   alpha = pow((hbar*hbar/m*kons),(0.25));
   lambda = (2.*m*alpha*alpha*E)/(hbar*hbar);
   rhoMin = 0.;
   rhoMaks = 100.;
   h = (rhoMaks-rhoMin)/static_cast<double>(nsteps);
   e = -1./(h*h);
   n = nsteps-1;
   mat A(n,n), B(n,n);
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
   Jacobi(A,B,n);
   cout << A << endl;
}



void
Local::Jacobi(const mat& A,
            mat& B,
            const int& n)
{
   //setting up eigenvector-matrix:
   B.eye();
   int k, l;
   double epsilon = 1.0e-8;
   maxOffDiag(A, k, l, n);
   double aMax = A(k,l);
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



void maxOffDiag(mat& A,
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
}



void rotate(mat& A,
            mat& B,
            const int& k,
            const int& l,
            const int& n)
{
   if(A(k,l) != 0.0)
   {
      double tau = (A(l,l)-A(k,k))/(2*A(k,l));
      if (tau < 0.0)
      {
            double t = -tau - sqrt(1. + tau*tau);
      }
      else
      {
            double t = -tau + sqrt(1. + tau*tau);
      }
      double c = 1.0/(sqrt(1.0+ t*t));
      double s = t*c;
   }
   else
   {
      double c = 1.0;
      double s = 0.0;
   }
   mat S;
   S.eye(n,n);
   S(k,k) = c;
   S(l,l) = c;
   S(k,l) = -s;
   S(l,k) = s;
   B = S.t()*A*S;
}


