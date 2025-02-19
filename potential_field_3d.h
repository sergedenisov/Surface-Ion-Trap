#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>

double computeSolidAngle(const std::vector<std::vector<double>>& rectangle, double obsX, double obsY, double obsZ) {
    // Ensure the rectangle has exactly 4 vertices
    if (rectangle.size() != 4) {
        std::cerr << "Rectangle must have exactly 4 vertices." << std::endl;
        return 0.0;
    }

    // Translate rectangle vertices relative to the observer and normalize
    std::vector<std::vector<double>> unitVertices;
    for (const auto& vertex : rectangle) {
        double dx = vertex[0] - obsX;
        double dy = vertex[1] - obsY;
        double dz = 0.0 - obsZ; // Rectangle is on the xy-plane (z = 0)
        double mag = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (mag == 0.0) {
            std::cerr << "Observer is too close to the rectangle." << std::endl;
            //return 2*3.14;
        }
        unitVertices.push_back({dx / mag, dy / mag, dz / mag});
    }

    // Apply the Oosterom and Strackee formula
    double solidAngle=0;
    std::vector<double> v1 = unitVertices[0];
    std::vector<double> v2 = unitVertices[1];
    std::vector<double> v3 = unitVertices[2];

    //First solid angle v1, v2, v4
    double numerator = 0.0;
    double denominator = 1.0;

    numerator += v1[0]*(v2[1]*v3[2] - v2[2]*v3[1]) - v1[1]*(v2[0]*v3[2] - v2[2]*v3[0]) + v1[2]*(v2[0]*v3[1] - v2[1]*v3[0]);
    denominator += (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) + (v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2]) + (v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    if(numerator>0 && denominator<0){
        solidAngle += std::abs(std::atan2(numerator, denominator) + 3.14159);
        std::cout<<1<<std::endl;
    } else{
        solidAngle += std::abs(std::atan2(numerator, denominator));
    }
    //Second solid angle v2, v3, v4
    numerator = 0.0;
    denominator = 1.0;
    v1 = unitVertices[0];
    v2 = unitVertices[2];
    v3 = unitVertices[3];

    numerator += v1[0]*(v2[1]*v3[2] - v2[2]*v3[1]) - v1[1]*(v2[0]*v3[2] - v2[2]*v3[0]) + v1[2]*(v2[0]*v3[1] - v2[1]*v3[0]);
    denominator += (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) + (v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2]) + (v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2]);

    if(numerator>0 && denominator<0){
        solidAngle += std::abs(std::atan2(std::abs(numerator), denominator) + 3.14159);
                std::cout<<1<<std::endl;

    } else{
        solidAngle += std::abs(std::atan2(std::abs(numerator), denominator));
    }
    //solidAngle += 2*std::atan2(numerator, denominator);
    double sA = std::abs(solidAngle)/(3.14159);
    if(sA > 1.0){
        return 1.0;
    }
    return sA;
}



struct Point2D {
    double x, y;
};


// New Function to generate the potential field in 3D
std::vector<std::vector<std::vector<double>>> generatePotentialField3D_new(
    const std::vector<double>& e_positions,
    const std::vector<double>& e_voltages,
    const std::vector<double>& e_widths, double lx, double ly, double lz,
    int nx = 80, int ny = 80, int nz = 80) {

    // Find vertices positions of electrodes
    
    std::vector<std::vector<double>> e_rect0 = {
        {e_positions[0] - e_widths[0] / 2, 0},
        {e_positions[0] + e_widths[0] / 2, 0},
        {e_positions[0] + e_widths[0] / 2, ly},
        {e_positions[0] - e_widths[0] / 2, ly}
    };
    std::vector<std::vector<double>> e_rect1 = {
        {e_positions[1] - e_widths[1] / 2, 0},
        {e_positions[1] + e_widths[1] / 2, 0},
        {e_positions[1] + e_widths[1] / 2, ly},
        {e_positions[1] - e_widths[1] / 2, ly}
    };
    std::vector<std::vector<double>> e_rect2 = {
        {e_positions[2] - e_widths[2] / 2, 0},
        {e_positions[2] + e_widths[2] / 2, 0},
        {e_positions[2] + e_widths[2] / 2, ly},
        {e_positions[2] - e_widths[2] / 2, ly}
    };
    std::vector<std::vector<double>> e_rect3 = {
        {e_positions[3] - e_widths[3] / 2, 0},
        {e_positions[3] + e_widths[3] / 2, 0},
        {e_positions[3] + e_widths[3] / 2, ly},
        {e_positions[3] - e_widths[3] / 2, ly}
    };
    std::vector<std::vector<double>> e_rect4 = {
        {e_positions[4] - e_widths[4] / 2, 0},
        {e_positions[4] + e_widths[4] / 2, 0},
        {e_positions[4] + e_widths[4] / 2, ly},
        {e_positions[4] - e_widths[4] / 2, ly}
    };

    double dx = lx / (nx), dy = ly / (ny), dz = lz / (nz);
    std::vector<std::vector<std::vector<double>>> V(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    // Initialize boundary conditions for electrodes
    int c=0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double x = dx*i;
                double y = dy*j;
                double z = dz*k;
                V[i][j][k] = (e_voltages[0] * computeSolidAngle(e_rect0, x, y, z) + e_voltages[1] * computeSolidAngle(e_rect1, x, y, z) + e_voltages[2] * computeSolidAngle(e_rect2, x, y, z) + e_voltages[3] * computeSolidAngle(e_rect3, x, y, z)+ e_voltages[4] * computeSolidAngle(e_rect4, x, y, z));
                
            }
        }
        
    }
    //std::cout<<computeSolidAngle(e_rect_test, 320.513, 320.513 ,12.8205)/(2*3.14)<<std::endl;
    //std::cout<<computeSolidAngle(e_rect1, 0, 0, 10)<<std::endl;;
    //std::cout<<computeSolidAngle(e_rect2, 0, 0, 10)<<std::endl;;
    //std::cout<<computeSolidAngle(e_rect3, 0, 0, 10)<<std::endl;;

    return V;
}

 

// Function to find the electric field in the x-direction (3D)
std::vector<std::vector<std::vector<double>>> generateElectricFieldX3D(
    const std::vector<std::vector<std::vector<double>>>& potential_field, double lx = 1200e-6,
    int nx = 80, int ny = 80, int nz = 80) {
    
    double dx = lx / (nx);
    std::vector<std::vector<std::vector<double>>> Ex(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    // Calculate Ex using finite differences
    for (int i = 0; i < nx - 1; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Ex[i][j][k] = -(potential_field[i+1][j][k] - potential_field[i][j][k]) / (dx);
            }
        }
    }
    return Ex;
}

// Function to find the electric field in the y-direction (3D)
std::vector<std::vector<std::vector<double>>> generateElectricFieldY3D(
    const std::vector<std::vector<std::vector<double>>>& potential_field, double ly = 1200e-6,
    int nx = 80, int ny = 80, int nz = 80) {
    
    double dy = ly / (ny);
    std::vector<std::vector<std::vector<double>>> Ey(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    // Calculate Ey using finite differences
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny - 1; ++j) {
            for (int k = 0; k < nz; ++k) {
                Ey[i][j][k] = -(potential_field[i][j+1][k] - potential_field[i][j][k]) / (dy);
            }
        }
    }
    return Ey;
}

// Function to find the electric field in the z-direction (3D)
std::vector<std::vector<std::vector<double>>> generateElectricFieldZ3D(
    const std::vector<std::vector<std::vector<double>>>& potential_field, double lz = 1200e-6,
    int nx = 80, int ny = 80, int nz = 80) {
    
    double dz = lz / (nz);
    std::vector<std::vector<std::vector<double>>> Ez(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));

    // Calculate Ez using finite differences
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz - 1; ++k) {
                Ez[i][j][k] = -(potential_field[i][j][k+1] - potential_field[i][j][k]) / (dz);
            }
        }
    }
    return Ez;
}

