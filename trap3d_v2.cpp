#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include </home/serge/Desktop/data/potential_field_3d.h>
#include </home/serge/Desktop/data/pseudo_potential_field_3d.h>
//#include </home/serge/Desktop/data/find_local_min.h>
#include </home/serge/Desktop/data/parameter_finder.h>



// Main function to generate dataset in 3D
int main() {
    const double Freq = 2*3.14*22*std::pow(10, 6);//MHz as in radial
    const double Q = 1.6*std::pow(10.0, -19.0);//Ion charge
    const double M = 28.8*std::pow(10, -26);//Mass of Yb atom

    double lx = 600e-6;  // Domain width in meters (1000 micrometers)
    double ly = 10000e-6;  // Domain height in meters (1000 micrometers)
    double lz = 200e-6;   // Domain depth in meters (350 micrometers)
    /*
    double lx = 1200e-6;  // Domain width in meters (1000 micrometers)
    double ly = 1200e-6;  // Domain height in meters (1000 micrometers)
    double lz = 500e-6;   // Domain depth in meters (350 micrometers)*/

    int nx = 150, ny = 150, nz = 150;

    std::ofstream inputsFile("/home/serge/Desktop/data/data3d/inputs_trap_3d.csv");
    std::ofstream outputsFile_p("/home/serge/Desktop/data/data3d/outputs_trap_p_3d.csv");
    std::ofstream outputsFile_e_x("/home/serge/Desktop/data/data3d/outputs_trap_e_x_3d.csv");
    std::ofstream outputsFile_e_y("/home/serge/Desktop/data/data3d/outputs_trap_e_y_3d.csv");
    std::ofstream outputsFile_e_z("/home/serge/Desktop/data/data3d/outputs_trap_e_z_3d.csv");
    std::ofstream outputsFile_pp("/home/serge/Desktop/data/data3d/outputs_trap_pp_3d.csv");

    if (!inputsFile.is_open() || !outputsFile_p.is_open() || !outputsFile_e_x.is_open() || 
        !outputsFile_e_y.is_open() || !outputsFile_e_z.is_open() || !outputsFile_pp.is_open()) {
        std::cerr << "Error opening output files!" << std::endl;
        return 1;
    }

    // Electrode positions (in micrometers), voltages, and widths (in micrometers)
    std::vector<double> positions = {0e-6, 80e-6, 180e-6, 400e-6, 600e-6};  // Positions along x-axis in micrometers
    std::vector<double> voltages = {0, 100, 0, 100, 0};        // Voltages in volts
    std::vector<double> widths = {0e-6, 160e-6, 40e-6, 400e-6, 0e-6};       // Widths in micrometers
    //std::vector<double> positions = {120e-6, 360e-6, 600e-6, 840e-6, 1080e-6};//In micrometers
    //std::vector<double> widths = {240e-6, 240e-6, 240e-6, 240e-6, 240e-6};//In micrometers
    
    // Generate potential field
    std::vector<std::vector<std::vector<double>>> potential_field = generatePotentialField3D_new(positions, voltages, widths, lx, ly, lz, nx, ny, nz);


    // Generate electric fields
    std::vector<std::vector<std::vector<double>>> electric_field_x = generateElectricFieldX3D(potential_field, lx, nx, ny, nz);
    std::vector<std::vector<std::vector<double>>> electric_field_y = generateElectricFieldY3D(potential_field, ly, nx, ny, nz);
    std::vector<std::vector<std::vector<double>>> electric_field_z = generateElectricFieldZ3D(potential_field, lz, nx, ny, nz);

    // Generate pseudo potential field
    std::vector<std::vector<std::vector<double>>> pseudo_potential_field = generatePseudoField3D(electric_field_x, electric_field_y, electric_field_z, Q, M, Freq, nx, ny, nz);
    
    

    //Slice pseudo potential at y_slice
    int y_slice = ny/2;
    std::vector<std::vector<double>> sliced_pseudo_potential_field = slice(pseudo_potential_field, y_slice, nx, nz);


    //R_0 - min distance from hole to electrodes
    double R0 = holeHeight(sliced_pseudo_potential_field, lz, nz);
    std::cout<<"R_0 = "<<R0*1000000<<" um"<<std::endl;


    //Find depth
    double depth = findDepth(sliced_pseudo_potential_field);
    std::cout<<"Depth = "<<depth<<" meV"<<std::endl;

    //Find secular frequency axis Z
    double freq_sec_z = secFreqZ(sliced_pseudo_potential_field, lx, lz, M);
    std::cout<<"secular frequency = "<<freq_sec_z<<std::endl;

    //Find q parameter
    double q = qParam(freq_sec_z, Freq, M, R0, Q,  voltages[1]);
    std::cout<<"q = "<<q<<std::endl;
    std::cout<<"e = "<<effParam(freq_sec_z, Freq, M, R0, Q,  voltages[1])<<std::endl;

  /*
    
    //Find secular frequency axis Z and efficiency parameter

    double dz = lz/(nz-1);
    double dPhi_1 = (pseudo_potential_field[pp_hole[0]][y_slice][pp_hole[1]] - pseudo_potential_field[pp_hole[0]][y_slice][pp_hole[1]-1])/dz;
    double dPhi_2 = (pseudo_potential_field[pp_hole[0]][y_slice][pp_hole[1]+1] - pseudo_potential_field[pp_hole[0]][y_slice][pp_hole[1]])/dz;
    double kz = (std::abs((dPhi_2 - dPhi_1)/dz)) / (std::pow((1 + dPhi_1*dPhi_1),1.5));//Curvature
    double freq_sec_z = std::sqrt(kz/M);//Secular frequency axis Z 
    double e = freq_sec_z * M * Freq * R0*R0 * std::sqrt(2) / (Q* voltages[1]);//Efficiency parameter
    double q = e*2*Q*voltages[1]/(M*R0*R0*Freq*Freq);//q parameter
    
    double q0 = 2*Q*voltages[1]/(M*R0*R0*Freq*Freq);//q parameter for ideal trap
    
    std::cout<<"Secular frequency: "<<freq_sec_z<<std::endl;
    std::cout<<e<<" "<<2*std::sqrt(2) * freq_sec_z/Freq<<std::endl;
    std::cout<<q<<" "<<q0<<std::endl;
    */
    // Write inputs to file (positions + voltages + widths)
    inputsFile << Freq << ",";
    
    for (size_t i = 0; i < positions.size(); ++i) {
        inputsFile << positions[i] << ",";
    }
    for (size_t i = 0; i < voltages.size(); ++i) {
        inputsFile << voltages[i] << ",";
    }
    for (size_t i = 0; i < widths.size(); ++i) {
        inputsFile << widths[i] << (i == widths.size() - 1 ? "\n" : ",");
    }
    // Write potential field (flattened) to file
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                outputsFile_p << potential_field[i][j][k] << (k == nz - 1 ? "" : ",");
            }
            outputsFile_p << (j == ny - 1 ? "" : ",");
        }
        outputsFile_p << (i == nx - 1 ? "\n" : ",");
    }

    // Write electric fields (flattened) to file
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                outputsFile_e_x << electric_field_x[i][j][k] << (k == nz - 1 ? "" : ",");
                outputsFile_e_y << electric_field_y[i][j][k] << (k == nz - 1 ? "" : ",");
                outputsFile_e_z << electric_field_z[i][j][k] << (k == nz - 1 ? "" : ",");
            }
            outputsFile_e_x << (j == ny - 1 ? "" : ",");
            outputsFile_e_y << (j == ny - 1 ? "" : ",");
            outputsFile_e_z << (j == ny - 1 ? "" : ",");
        }
        outputsFile_e_x << (i == nx - 1 ? "\n" : ",");
        outputsFile_e_y << (i == nx - 1 ? "\n" : ",");
        outputsFile_e_z << (i == nx - 1 ? "\n" : ",");
    }

    // Write pseudo potential field (flattened) to file
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                outputsFile_pp << pseudo_potential_field[i][j][k] << (k == nz - 1 ? "" : ",");
            }
            outputsFile_pp << (j == ny - 1 ? "" : ",");
        }
        outputsFile_pp << (i == nx - 1 ? "\n" : ",");
    }

    inputsFile.close();
    outputsFile_p.close();
    outputsFile_e_x.close();
    outputsFile_e_y.close();
    outputsFile_e_z.close();
    outputsFile_pp.close();

    std::cout << "3D dataset generation complete. Files saved as inputs_trap_3d.csv, outputs_trap_p_3d.csv, etc." << std::endl;
    return 0;
}