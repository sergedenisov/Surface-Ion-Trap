#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include <random>
#include </home/serge/Desktop/data/potential_field_3d.h>
#include </home/serge/Desktop/data/pseudo_potential_field_3d.h>
#include </home/serge/Desktop/data/parameter_finder.h>

int main() {
    int numsamples = 4; // Number of samples to generate
    const double Freq = 2*3.14*22*1e6;
    const double Q = 1.6e-19;
    const double M = 28.8e-26;

    // Domain settings
    double lx = 600e-6, ly = 10000e-6, lz = 200e-6;
    int nx = 150, ny = 150, nz = 150;

    // Random number setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> volt_dist(50.0, 150.0);
    std::uniform_real_distribution<double> pos_dist(0, lx);
    // Open output files (truncate mode)
    std::ofstream inputsFile("/home/serge/Desktop/data/data3d/inputs_trap_3d.csv");
    std::ofstream outputsFile_p("/home/serge/Desktop/data/data3d/outputs_trap_p_3d.csv");
    std::ofstream outputsFile_e_x("/home/serge/Desktop/data/data3d/outputs_trap_e_x_3d.csv");
    std::ofstream outputsFile_e_y("/home/serge/Desktop/data/data3d/outputs_trap_e_y_3d.csv");
    std::ofstream outputsFile_e_z("/home/serge/Desktop/data/data3d/outputs_trap_e_z_3d.csv");
    std::ofstream outputsFile_pp("/home/serge/Desktop/data/data3d/outputs_trap_pp_3d.csv");

    // Write CSV headers
    inputsFile << "pos0,pos1,pos2,pos3,pos4,volt0,volt1,volt2,volt3,volt4,width0,width1,width2,width3,width4\n";

    for(int sample = 0; sample < numsamples; ++sample) {
        // Generate random parameters
        std::vector<double> positions = {pos_dist(gen), pos_dist(gen), pos_dist(gen), pos_dist(gen), pos_dist(gen)};
        std::vector<double> widths = {0e-6, 160e-6, 40e-6, 400e-6, 0e-6};
        //std::vector<double> voltages = {0, volt_dist(gen), 0, volt_dist(gen), 0};
        std::vector<double> voltages = {volt_dist(gen), volt_dist(gen), volt_dist(gen), volt_dist(gen), volt_dist(gen)};

        // Generate fields
        auto potential_field = generatePotentialField3D_new(positions, voltages, widths, lx, ly, lz, nx, ny, nz);
        auto electric_field_x = generateElectricFieldX3D(potential_field, lx, nx, ny, nz);
        auto electric_field_y = generateElectricFieldY3D(potential_field, ly, nx, ny, nz);
        auto electric_field_z = generateElectricFieldZ3D(potential_field, lz, nx, ny, nz);
        auto pseudo_potential_field = generatePseudoField3D(electric_field_x, electric_field_y, electric_field_z, Q, M, Freq, nx, ny, nz);

        // Parameters (not written to any files)
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


        // Write input parameters
        inputsFile << ","
                  << positions[0] << "," << positions[1] << "," << positions[2] << "," 
                  << positions[3] << "," << positions[4] << ","
                  << voltages[0] << "," << voltages[1] << "," << voltages[2] << "," 
                  << voltages[3] << "," << voltages[4] << ","
                  << widths[0] << "," << widths[1] << "," << widths[2] << "," 
                  << widths[3] << "," << widths[4] << "\n";

        // Write output data (flattened 3D arrays)
        auto writeField = [](auto& field, auto& file) {
            for (int i = 0; i < field.size(); ++i) {
                for (int j = 0; j < field[i].size(); ++j) {
                    for (int k = 0; k < field[i][j].size(); ++k) {
                        file << field[i][j][k];
                        if (!(i == field.size()-1 && j == field[i].size()-1 && k == field[i][j].size()-1)) {
                            file << ",";
                        }
                    }
                }
            }
            file << "\n";
        };

        writeField(potential_field, outputsFile_p);
        //writeField(electric_field_x, outputsFile_e_x);
        //writeField(electric_field_y, outputsFile_e_y);
        //writeField(electric_field_z, outputsFile_e_z);
        writeField(pseudo_potential_field, outputsFile_pp);
    }

    std::cout << "Generated " << numsamples << " samples.\n";
    return 0;
}