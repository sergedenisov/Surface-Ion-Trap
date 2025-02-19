#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include </home/serge/Desktop/data/calculas.h>


// Function to find h distance from hole to electrodes
double holeHeight(const std::vector<std::vector<double>>& field, double lz, int nz){
    std::vector<int> pp_hole = findLocalMinimum(field);
    return pp_hole[1]*lz/(nz);
}

// Function to find the depth in a 2D scalar field
double findDepth(const std::vector<std::vector<double>>& field) {
    int rows = field.size();
    if (rows == 0) throw std::invalid_argument("Field is empty");
    int cols = field[0].size();

    std::vector<int> pp_hole = findLocalMinimum(field);
    int hole_x = pp_hole[0],hole_z = pp_hole[1];
    
    std::vector<double> hole_line(cols, 0.0);

    int square = 10;
    int escape = 0;
    for(int i=0; i<cols; ++i){
        hole_line[i] = field[hole_x][i];
    }
    for(int i=hole_z+square; i<cols-square; ++i){
        bool check = true;

        for(int c=1; c<square; ++c){
            if(hole_line[i-square] > hole_line[i]){
                check = false;
                break;
            }
        }
        if(check){
            escape = i;
        }
    }
    double depth = (field[hole_x][escape] - field[hole_x][hole_z]) *1000 / (1.6*std::pow(10, -19));//Depth in meV

    return depth;

}

// Function to find secular frequency axis X
double secFreqX(const std::vector<std::vector<double>>& field,double lx, double lz, double M){
    std::vector<int> pp_hole = findLocalMinimum(field);
    double k_x = computeCurvaturesX(field, lx, lz, pp_hole[0], pp_hole[1]);
    return std::sqrt(k_x/M);
}
// Function to find secular frequency axis Z
double secFreqZ(const std::vector<std::vector<double>>& field,double lx, double lz, double M){
    std::vector<int> pp_hole = findLocalMinimum(field);
    double k_z = computeCurvaturesZ(field, lx, lz, pp_hole[0], pp_hole[1]);
    return std::sqrt(k_z/M);
}
// Function to find Efficiency parameter
double effParam(double freq_sec_z, double Freq, double M, double R0, double Q, double v){
    double e = freq_sec_z * M * Freq * R0*R0 * std::sqrt(2) / (Q* v);
    return e;
}
// Function to find q parameter
double qParam(double freq_sec_z, double Freq, double M, double R0, double Q, double v){
    double e = effParam(freq_sec_z, Freq, M, R0, Q, v);
    double q = e*2*Q*v/(M*R0*R0*Freq*Freq);
    return q;
}
