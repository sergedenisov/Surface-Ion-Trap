#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>

// Function to find the pseudo potential field in 3D
std::vector<std::vector<std::vector<double>>> generatePseudoField3D(
    const std::vector<std::vector<std::vector<double>>>& electric_field_x,
    const std::vector<std::vector<std::vector<double>>>& electric_field_y,
    const std::vector<std::vector<std::vector<double>>>& electric_field_z,
    const double Q, const double M, const double Freq,
    int nx = 80, int ny = 80, int nz = 80) {
    
    std::vector<std::vector<std::vector<double>>> pseudo_potential_field(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                pseudo_potential_field[i][j][k] = (Q * Q) * (electric_field_x[i][j][k] * electric_field_x[i][j][k] +
                                                             electric_field_y[i][j][k] * electric_field_y[i][j][k] +
                                                             electric_field_z[i][j][k] * electric_field_z[i][j][k]) / (4 * M * Freq * Freq);
                /*if(pseudo_potential_field[i][j][k]>0.000002){
                    pseudo_potential_field[i][j][k] = 0.000002;
                }*/
            }
        }
    }
    return pseudo_potential_field;
}

// New Function to slice the pseudo potential field in y_slice
std::vector<std::vector<double>> slice(const std::vector<std::vector<std::vector<double>>>& pseudo_potential_field, int y_slice = 40, int nx = 80, int nz = 80) {

    std::vector<std::vector<double>> sliced_pseudo_potential_field(nx, std::vector<double> (nz, 0.0));
    for(int i=0; i<nx; ++i){
        for(int j=0; j<nz; ++j){
            sliced_pseudo_potential_field[i][j] = pseudo_potential_field[i][y_slice][j];
        }
    }

    return sliced_pseudo_potential_field;
}