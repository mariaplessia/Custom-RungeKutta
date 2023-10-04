// RK4.cpp
//
// This is a customizable Runge-Kutta implementation that solves 1st order ODEs for simple & custom mechanical systems
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdio.h>
using namespace std;


double f1(double t, double x1, double x2, double x3, double x4, double c);    
double f2(double t, double x1, double x2, double x3, double x4, double c);
double f3(double t, double x1, double x2, double x3, double x4, double c);
double f4(double t, double x1, double x2, double x3, double x4, double c);

ofstream data_out("results.txt"); 

int main() {
    int i, j;
    double ttotal = 20.; 
    double delt = 0.001;
    
    const int n = int(ttotal/delt);

    double k1[4], k2[4], k3[4], k4[4];
    double x[4] = {0.0, 0.0, 0.0, 0.0}; // Initial conditions
    double c_values[3] = {0.1, 5.0, 10.0}; // Different damping coefficients

    for (int c_index = 0; c_index < 4; c_index++) {
        double t = 0.;
        double c = c_values[c_index];
        data_out << "Damping Coefficient = : " << c << endl;

        // We need to reset initial conditions for each iteration
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 1.0;
        x[3] = 0.0;
        
        data_out << t << '\t' << x[0] << '\t' << x[1] << '\n' endl;

        for(i = 1; i <= n; i++) {

            k1[0] = delt * f1(t, x[0], x[1], x[2], x[3], c);
            k1[1] = delt * f2(t, x[0], x[1], x[2], x[3], c);
            k1[2] = delt * f3(t, x[0], x[1], x[2], x[3], c);
            k1[3] = delt * f4(t, x[0], x[1], x[2], x[3], c);

            k2[0] = delt * f1(t + 0.5*delt, x[0] + 0.5*k1[0], x[1] + 0.5*k1[1], x[2] + 0.5*k1[2], x[3] + 0.5*k1[3], c);
            k2[1] = delt * f2(t + 0.5*delt, x[0] + 0.5*k1[0], x[1] + 0.5*k1[1], x[2] + 0.5*k1[2], x[3] + 0.5*k1[3], c);
            k2[2] = delt * f3(t + 0.5*delt, x[0] + 0.5*k1[0], x[1] + 0.5*k1[1], x[2] + 0.5*k1[2], x[3] + 0.5*k1[3], c);
            k2[3] = delt * f4(t + 0.5*delt, x[0] + 0.5*k1[0], x[1] + 0.5*k1[1], x[2] + 0.5*k1[2], x[3] + 0.5*k1[3], c);

            k3[0] = delt * f1(t + 0.5*delt, x[0] + 0.5*k2[0], x[1] + 0.5*k2[1], x[2] + 0.5*k2[2], x[3] + 0.5*k2[3], c);
            k3[1] = delt * f2(t + 0.5*delt, x[0] + 0.5*k2[0], x[1] + 0.5*k2[1], x[2] + 0.5*k2[2], x[3] + 0.5*k2[3], c);
            k3[2] = delt * f3(t + 0.5*delt, x[0] + 0.5*k2[0], x[1] + 0.5*k2[1], x[2] + 0.5*k2[2], x[3] + 0.5*k2[3], c);
            k3[3] = delt * f4(t + 0.5*delt, x[0] + 0.5*k2[0], x[1] + 0.5*k2[1], x[2] + 0.5*k2[2], x[3] + 0.5*k2[3], c);

            k4[0] = delt * f1(t + delt, x[0] + k3[0], x[1] + k3[1], x[2] + k3[2], x[3] + k3[3], c);
            k4[1] = delt * f2(t + delt, x[0] + k3[0], x[1] + k3[1], x[2] + k3[2], x[3] + k3[3], c);
            k4[2] = delt * f3(t + delt, x[0] + k3[0], x[1] + k3[1], x[2] + k3[2], x[3] + k3[3], c);
            k4[3] = delt * f4(t + delt, x[0] + k3[0], x[1] + k3[1], x[2] + k3[2], x[3] + k3[3], c);

            for (j = 0; j < 4; j++) {
                x[j] += (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]) / 6.0;
            }

            t += delt;
            data_out << t << '\t' << x[0] << '\t' << x[1] << endl;
            
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double f1(double t, double x1, double x2, double x3, double x4, double c) {
    return x3; // f1 = x3
}

double f2(double t, double x1, double x2, double x3, double x4, double c) {
    return x4; // f2 = x4
}

double f3(double t, double x1, double x2, double x3, double x4, double c) {
    return (-3.0 * x1 + x2) - c*x4; // f3 = -3*x1 + x2 + c*x3
}

double f4(double t, double x1, double x2, double x3, double x4, double c) {
    return (1.0/3.0)*x1 - x2 - c*x4; // f4 = (1/3)*x1 - x2 + c*x4
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
