#include "system.h"

System::System()
    :numberOfObject(0)
{
}

System::~System()
{
}


void
System::addObject(arma::Col<double> pos,
                       arma::Col<double> vel,
                       double m,
                       char* name)
{
    Object x = Object(pos, vel, m, name);
    objectlist[numberOfObject] = x;
    this->numberOfObject+=1;
}



arma::Col<double>
System::force(arma::Col<double> pos1,
              arma::Col<double> pos2,
              double GM)
{
    //pos1: position object 1, pos 2: position object 2, GM: gravitational constant times mass of object 2
    double R = sqrt((pos1(0)-pos2(0))*(pos1(0)-pos2(0))+(pos1(1)-pos2(1))*(pos1(1)-pos2(1)));
    arma::Col<double> F(2);
    F(0) = GM*(pos2(0)-pos1(0))/(R*R*R);
    F(1) = GM*(pos2(1)-pos1(1))/(R*R*R);
    return F;
}



void
System::solve(const std::string filename)
{
    std::ofstream fout(filename.c_str());

//    //fout.setf(std::ios::scientific);
//    fout.precision(8);
//    //fout.width(8);
//    fout << "describe = '[t x v_x y v_y]'" << "\n";
//    fout << "A = [ " ;

//    fout << i*delta_t << "\t\t" << xout(0) << "\t\t" << xout(1) << "\t\t" << yout(0) << "\t\t" << yout(1) << "\n";
//    fout << "];" << "\n\n";
//    fout << "plot (A(:,2),A(:,4))" << "\n\n";
    fout.close();
}
