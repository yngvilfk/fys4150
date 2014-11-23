#include<randomgenerator.h>
#include<system.h>
#include<object.h>
#include<solvestep.h>

int main()
{
   RandomGenerator generate;
   double R0 = 20;
   System mysystem = generate.randomSystem(100, R0);
   mysystem.setEpsilon(0.2);
   SolveStep solver(mysystem);
   solver.solve(100, 1000, "solveR_", "rk4");

   std::ofstream fout("start.m");
   fout << "A = [";
   for (int i = 0 ; i < mysystem.numberOfObject ; ++i)
   {
      Object &mainbody = mysystem.objectlist[i];

      fout << mainbody.getPosition()(0) << "\t\t" << mainbody.getPosition()(1)
           << "\t\t" << mainbody.getPosition()(2) << "\n";
   }
   fout << "] \n";
   fout << "plot3(A(:,1), A(:,2),A(:,3), 'o')";
   fout.close();
   std::cout << mysystem.numberOfObject << std::endl;
}
