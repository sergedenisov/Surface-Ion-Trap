#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>

// Function to find the local minimum in a 2D scalar field
std::vector<int> findLocalMinimum(const std::vector<std::vector<double>>& field) {
    int rows = field.size();
    if (rows == 0) throw std::invalid_argument("Field is empty");
    int cols = field[0].size();

    // Directions for 4-connected neighbors (up, down, left, right)
    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};
    int i0=-1, j0=-1;
    int square = 10;
    // Iterate through the field (excluding boundaries)
    for (int i = square; i < rows - square; ++i) {
        for (int j = square; j < cols - square; ++j) {
            bool isLocalMin = true;
            // Check all 4 neighbors
            for (int k = 0; k < 4; ++k) {
                for(int c=1; c<=square; ++c){
                    int ni = i + c*dx[k];
                    int nj = j + c*dy[k];
                    if (field[ni][nj] < field[i][j]) {
                        isLocalMin = false;
                    }
                }
                
            }
            // If it's a local minimum, return its coordinates
            if (isLocalMin && i0==-1) {
                i0=i;
                j0=j;
                //return {i, j};
            }
            if (isLocalMin && field[i][j] < field[i0][j0]) {
                i0=i;
                j0=j;
                //return {i, j};
            }
        }
    }
    return {i0, j0};
    // If no local minimum is found, return invalid coordinates
    return {-1, -1};
}
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

// Compute first derivative in x-direction (central/forward/backward differences)
double compute_f_x(const vector<vector<double>>& grid, int i, int j, double dx) {
    int cols = grid[0].size();
    if (cols <= 1) return 0.0;
    if (j == 0) {
        return (grid[i][1] - grid[i][0]) / dx; // Forward difference
    } else if (j == cols - 1) {
        return (grid[i][j] - grid[i][j-1]) / dx; // Backward difference
    } else {
        return (grid[i][j+1] - grid[i][j-1]) / (2.0 * dx); // Central difference
    }
}

// Compute second derivative in x-direction (central/forward/backward differences)
double compute_f_xx(const vector<vector<double>>& grid, int i, int j, double dx) {
    int cols = grid[0].size();
    if (cols < 2) return 0.0;
    if (j == 0) {
        return (grid[i][2] - 2 * grid[i][1] + grid[i][0]) / (dx * dx);
    } else if (j == cols - 1) {
        return (grid[i][j] - 2 * grid[i][j-1] + grid[i][j-2]) / (dx * dx);
    } else {
        return (grid[i][j+1] - 2 * grid[i][j] + grid[i][j-1]) / (dx * dx);
    }
}

// Compute first derivative in y-direction (central/forward/backward differences)
double compute_f_y(const vector<vector<double>>& grid, int i, int j, double dy) {
    int rows = grid.size();
    if (rows <= 1) return 0.0;
    if (i == 0) {
        return (grid[1][j] - grid[0][j]) / dy; // Forward difference
    } else if (i == rows - 1) {
        return (grid[i][j] - grid[i-1][j]) / dy; // Backward difference
    } else {
        return (grid[i+1][j] - grid[i-1][j]) / (2.0 * dy); // Central difference
    }
}

// Compute second derivative in y-direction (central/forward/backward differences)
double compute_f_yy(const vector<vector<double>>& grid, int i, int j, double dy) {
    int rows = grid.size();
    if (rows < 2) return 0.0;
    if (i == 0) {
        return (grid[2][j] - 2 * grid[1][j] + grid[0][j]) / (dy * dy);
    } else if (i == rows - 1) {
        return (grid[i][j] - 2 * grid[i-1][j] + grid[i-2][j]) / (dy * dy);
    } else {
        return (grid[i+1][j] - 2 * grid[i][j] + grid[i-1][j]) / (dy * dy);
    }
}

// Compute curvatures k_x and k_y for the scalar field at point (i, j)
double computeCurvaturesX(
    const vector<vector<double>>& grid,
    double width,
    double height,
    int i,
    int j
) {
    int rows = grid.size();
    if (rows == 0) throw runtime_error("Empty grid");
    int cols = grid[0].size();
    if (cols == 0) throw runtime_error("Empty grid");
    if (i < 0 || i >= rows || j < 0 || j >= cols) throw out_of_range("Invalid indices");

    // Calculate grid spacing
    double dx = (cols > 1) ? (width / (cols - 1)) : 0.0;

    // Compute derivatives for k_x (along x-direction)
    double f_x = compute_f_x(grid, i, j, dx);
    double f_xx = compute_f_xx(grid, i, j, dx);
    double k_x = f_xx / pow(1.0 + f_x * f_x, 1.5);
    return k_x;
}
double computeCurvaturesZ(
    const vector<vector<double>>& grid,
    double width,
    double height,
    int i,
    int j
) {
    int rows = grid.size();
    if (rows == 0) throw runtime_error("Empty grid");
    int cols = grid[0].size();
    if (cols == 0) throw runtime_error("Empty grid");
    if (i < 0 || i >= rows || j < 0 || j >= cols) throw out_of_range("Invalid indices");

    // Calculate grid spacing
    double dy = (rows > 1) ? (height / (rows - 1)) : 0.0;

    // Compute derivatives for k_y (along y-direction)
    double f_y = compute_f_y(grid, i, j, dy);
    double f_yy = compute_f_yy(grid, i, j, dy);
    double k_y = f_yy / pow(1.0 + f_y * f_y, 1.5);
    return k_y;
}