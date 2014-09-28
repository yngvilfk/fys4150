
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
                  const int&         n,
                  arma::Col<double>& eigval);

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
                       const arma::Col<double>&,
                       const std::string& name,
                          const int& n);

      void writeMatlabMat(std::ofstream& write,
                       arma::Mat<double>&,
                       const std::string& name,
                          const int& n);

      void plot(const arma::Mat<double>&   eigvec,
                const arma::Col<double>&   rho,
                const std::string&         fileName,
                const int&                 n,
                const int&                 lambda);

   }
}



int main()
{
   int number = 1; //one(1) or two(2) electron system
   int method = 0; //armadillosystem(1) or jacobifunction(0)
   const int nsteps =500;
   const double rhoMaks = 5.0 ;
   int lambda = 1; //select eigenvalue: lambda 1(1), lambda 2(2) etc
   const std::string fileName("omega5_1.m");
   double omega = 0.01;


   double V = 0.;
   double d= 0.;
   const int n = nsteps-1;
   const double rhoMin = 0.;
   const double h = (rhoMaks-rhoMin)/static_cast<double>(nsteps-1);
   std::cout << h << std::endl;
   const double e = -1./(h*h);
   std::cout << "e: " << e << std::endl;
   arma::Col<double> rho(n+2);
   rho.zeros();
   rho(0) = rhoMin;
   rho(n+1)=rhoMaks;

   arma::Mat<double> A(n,n), B(n,n);
   A.zeros();
   clock_t start, finish;

   //fill matrix A

   if (number == 1)   //look at one electron
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



   else   //look at a system of two electron
   {
      //first row:
      rho(1) = rhoMin + h;
      V = rho(1)*rho(1)*omega*omega; // + 1/rho(1);
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
      //method is the local "jacobi" method
      start = clock();         //start timer
      arma::Col<double> eigval;
      Local::jacobi(A,B,n,eigval);
      finish = clock();        //stop timer
      std::cout << "time jacobi algorithm: " << static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << std::endl;

      arma::Mat<double> eigvec = B;
      arma::Col<double> eigvalSorted = sort(eigval);
      std::cout << "lambda0: " << eigvalSorted(0) << std::endl;
      std::cout << "lambda1: " << eigvalSorted(1) << std::endl;
      std::cout << "lambda2: " << eigvalSorted(2) << std::endl;
   }


   else if (method == 1)
   {
       //method is the armadillo "eig_sum" function
 start = clock();         //start timer
      vec eigval;
      arma::Mat<double> eigvec;
      eig_sym(eigval, eigvec, A);
      finish = clock();        //stop timer
      std::cout << "time armadillo eig_sym algorithm: " << static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << std::endl;
      Local::plot(eigvec, rho, fileName,n,lambda);

   }

}




void
Local::jacobi(arma::Mat<double>& A,
              arma::Mat<double>& B,
              const int& n,
              arma::Col<double>& eigval)
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
      rotate(A, B,k, l, n);
      aMax = maxOffDiag(A, k, l, n);
      //std::cout << "Amaks: " << aMax << std::endl;
      ++ iterations;
   }
   std::cout << "number of iterations needed: " << iterations << std::endl;
   //eigval = sort(A.diag());
   eigval = A.diag();
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
            l = i;
            k = j;
         }

      }
   }

   return max;
}




void
Local::rotate(mat& A,
              mat& B,
            int& k,
            int& l,
            const int& n)
{
    /*
     perform one rotation of the matrix and sets A(k,l)=A(l,k)=0
     A -- symmetric input matrix
     k -- index (line) of max element of A
     l -- index (coloumn) of max element of A
     n -- dimension of matrix A
     */
double c, s;
   if(A(k,l) != 0.0)
   {
      double t;
      double tau = (A(l,l)-A(k,k))/(2*A(k,l));
      if (tau < 0.0)
      {
            t = -1./(-tau + sqrt(1. + tau*tau) );
      }
      else
      {
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

   double a_kk, a_ll, a_ik, a_il, b_ik, b_il;
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
      b_ik = B(i,k);
      b_il = B(i,l);
      B(i,k) = c*b_ik - s*b_il;
      B(i,l) = c*b_il + s*b_ik;

   }
}



/*
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
   B = S*A*S.t();
   A = B;
   std::cout << "A(k,l): " << A(k,l) << "A(l,k): " << A(l,k) << std::endl;
}
*/
void
Local::writeMatlabVec(std::ofstream& write,
                      const arma::Col<double>& vec,
                      const std::string& name,
                      const int& n)
{
   write << name << " = [" << "\n";
   for (int iVec = 0; iVec < n; ++iVec)
   {
      write << vec(iVec) << "\n";
   }

   write << "]" << "\n\n";
}


void
Local::writeMatlabMat(std::ofstream& write,
                      arma::Mat<double>& mat,
                      const std::string& name,
                      const int& n)
{
   write << name << " = [" << "\n";
   for (int iVec = 0; iVec < n; ++iVec)
   {
       for (int jVec = 0; jVec < n; ++jVec)
       {
           if (jVec < n-1)
           {
               write << mat(iVec,jVec) << " ";
           }
           else
           {
               write << mat(iVec,jVec) << "\n";
           }
       }

   }

   write << "]" << "\n\n";
}


void
Local::plot(const arma::Mat<double>& eigvec,
            const arma::Col<double>& rho,
            const std::string& fileName,
            const int& n,
            const int& lambda)
{
   arma::Col<double> u(n+2);
   u(0) = 0.;
   u(n+1) = 0.;
   arma::Col<double> v=eigvec.col(lambda-1);

   for ( int i= 1; i<=n ; ++i)
   {
      u(i)=v(i-1)*v(i-1);
   }

   std::ofstream write(fileName.c_str());
   //write to file for plotting
   const std::string fileName_val("rho");
   Local::writeMatlabVec(write, rho, fileName_val,n);
   const std::string fileName_vec("u_squared");
   Local::writeMatlabVec(write, u, fileName_vec,n);
   write.close();
}



