
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

      double maxOffDiag(arma::Mat<double>& A,
                      int&                     k,
                      int&                     l,
                      const int&               n);

      void rotate(arma::Mat<double>& A,
                  arma::Mat<double>& B,
                  int&         k,
                  int&         l,
                  const int&         n);
      void writeMatlabVec(std::ofstream& write,
                       arma::Col<double>&,
                       const std::string& name);
      void writeMatlabMat(std::ofstream& write,
                       arma::Mat<double>&,
                       const std::string& name);
   }
}



int main()
{
   int number = 1; //one(1) or two(2) electron system
   int method = 1; //armadillosystem(1) or jacobifunction(0)


   double V = 0.;
   double d= 0.;
   const double rhoMaks = 10.;
   const int nsteps =50;
   const int n = nsteps-1;
   const int rhoMin = 0.;
   const double h = (rhoMaks-rhoMin)/static_cast<double>(nsteps);
   std::cout << h << std::endl;
   const double e = -1./(h*h);
   std::cout << "e: " << e << std::endl;
   arma::Col<double> rho(n+1);
   rho.zeros();

   arma::Mat<double> A(n,n), B(n,n);
   A.zeros();
   //fill matrix A
   if (number == 1)
   {

      //first row:
      rho(1) = rhoMin + h;
      V = rho(1)*rho(1);
      d = 2./(h*h) + V;
      A(0,0) = d;
      A(0,1) = e;
      //row n (n_steps-1)
      rho(n) = rhoMin + static_cast<double>(n)*h;
      V = rho(n)*rho(n);
      d = 2./(h*h) + V;
      std::cout << "V: " << V << ", d: " << d << std::endl;
      A(n-1,n-1) = d;
      A(n-1,n-2) = e;
      //row 2 to n-1 :
      for (int i = 1 ; i<n-1 ; ++i)
      {
         rho(i+1) = rhoMin + static_cast<double>(i+1)*h;
         V = rho(i+1)*rho(i+1);
         d = 2./(h*h) + V;
         A(i,i) = d;
         A(i,i+1) = e;
         A(i,i-1) = e;
      }
   }
   else
   {
      double omega = 0.01;
      //first row:
      rho(1) = rhoMin + h;
      V = rho(1)*rho(1)*omega*omega + 1/rho(1);
      d = 2./(h*h) + V;
      A(0,0) = d;
      A(0,1) = e;
      //row n (n_steps-1)
      rho(n) = rhoMin + static_cast<double>(n)*h;
      V = rho(n)*rho(n)*omega*omega + 1/rho(n);
      d = 2./(h*h) + V;
      std::cout << "V: " << V << ", d: " << d << std::endl;
      A(n-1,n-1) = d;
      A(n-1,n-2) = e;
      //row 2 to n-1 :
      for (int i = 1 ; i<n-1 ; ++i)
      {
         rho(i+1) = rhoMin + static_cast<double>(i+1)*h;
         V = rho(i+1)*rho(i+1)*omega*omega + 1/rho(i+1);
         d = 2./(h*h) + V;
         A(i,i) = d;
         A(i,i+1) = e;
         A(i,i-1) = e;
      }
   }
   std::cout << "n_steps: " << nsteps  << ", rho_max: " << rhoMaks << std::endl;


   if (method == 0)
   {
      Local::jacobi(A,B,n);
      for (int i = 0; i<16 ; ++i)
      {
         cout << A(i,i) << endl;
      }
      cout << "min: " << min(A.diag()) << endl;
   }


   else if (method == 1)
   {
      vec eigval;
      arma::Mat<double> eigvec;
      eig_sym(eigval, eigvec, A);
      std::cout << "eigval: " << eigval(0) << std::endl;
      std::cout << "eigval: " << eigval(1) << std::endl;
      std::cout << "eigval: " << eigval(2) << std::endl;
      //std::cout << "eigvec: " << eigvec << std::endl;

      std::string fileName("test.m");
      std::ofstream write(fileName.c_str());
      //write to file for plotting
      const std::string fileName_val("eigval");
      Local::writeMatlabVec(write, eigval, fileName_val);
      const std::string fileName_vec("eigvec");
      Local::writeMatlabMat(write, eigvec, fileName_vec);
      write.close();
   }

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
   //std::cout << "Amaks out of loop: " << aMax << std::endl;
   int maxIterations = n*n*n;
   int iterations = 0;
   while ( aMax > epsilon && iterations < maxIterations)
   {
      rotate(A, B, k, l, n);
      aMax = maxOffDiag(A, k, l, n);
      //std::cout << "Amaks: " << aMax << std::endl;
      ++ iterations;
   }
   std::cout << "number of iterations needed: " << iterations << std::endl;

   return;
}



double
Local::maxOffDiag(arma::Mat<double>& A,
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
   //std::cout << "maks: " << max << std::endl;

   return max;
}



/*
void
Local::rotate(mat& A,
            mat& B,
            int& k,
            int& l,
            const int& n)
{
      //std::cout << "rotate" << std::endl;
double c, s;
   if(A(k,l) != 0.0)
   {
      double t;
      double tau = (A(l,l)-A(k,k))/(2*A(k,l));
      if (tau < 0.0)
      {
            //t = -tau - sqrt(1. + tau*tau);
            t = -1./(-tau + sqrt(1. + tau*tau) );
      }
      else
      {
         //t = -tau + sqrt(1. + tau*tau);
         t = 1./(tau +sqrt(1. + tau*tau) );
      }
      c = 1.0/(sqrt(1.0 + t*t));
      s = t*c;
   }
   else
   {
      c = 1.0;
      s = 0.0;
   }

   double a_kk, a_ll, a_ik, a_il;
   a_kk = A(k,k);
   a_ll = A(l,l);
   A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
   A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
   A(k,l) = 0.0;
   A(l,k) = 0.0;
   for (int i=0 ; i<n ; ++i)
   {
      if ( i != k && i != l)
      {
         a_ik = A(i,k);
         a_il = A(i,l);
         A(i,k) = c*a_ik - s*a_il;
         A(k,i) = A(i,k);
         A(i,l) = c*a_il + s*a_ik;
         A(l,i) = A(i,l);
      }
//      b_ik = B(i,k);
//      b_il = B(i,l);
//      B(i,k) = c*b_ik - s*b_il;
//      B(i,l) = c*b_il + s*b_ik;

   }
}
*/



void
Local::rotate(mat& A,
            mat& B,
            int& k,
            int& l,
            const int& n)
{
double c, s;
   if(A(k,l) != 0.0)
   {
      double t;
      double tau = (A(l,l)-A(k,k))/(2*A(k,l));
      if (tau < 0.0)
      {
            //t = -tau - sqrt(1. + tau*tau);
            t = -1./(-tau +sqrt(1.+tau*tau) );
      }
      else
      {
         //t = -tau + sqrt(1. + tau*tau);
         t = 1./(tau +sqrt(1.+tau*tau) );
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
   A = B;
}

void
Local::writeMatlabVec(std::ofstream& write,
                   arma::Col<double>& vec,
                   const std::string& name)
{
   write << name << " = [" << "\n";
   for (std::size_t iVec = 0; iVec < vec.size(); ++iVec)
   {
      write << vec(iVec) << "\n";
   }

   write << "]" << "\n\n";
}


void
Local::writeMatlabMat(std::ofstream& write,
                   arma::Mat<double>& mat,
                   const std::string& name)
{
   write << name << " = [" << "\n";
   for (std::size_t iVec = 0; iVec < mat.size(); ++iVec)
   {
      write << mat(iVec) << "\n";
   }

   write << "]" << "\n\n";
}
