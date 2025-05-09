#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <array>
#include "Eigen/Dense"


// Define the type for real numbers
// typedef double real;

// Define the type for a vector of real numbers
typedef std::vector<double> Vector;

// Define the type for a function pointer that takes a real number and returns a real number
typedef double (*FunctionPtr)(double, double);

// Define the type for a function pointer that takes a real number and a vector and returns a vector of real numbers
typedef Vector (*FunctionPtrVec)(double, Vector);

// Define the type for a solution containing time and solution values
struct Solution {
    Vector time;
    Vector sol;
    std::vector<Vector> sol_set; // For system of equations
};

// Function declarations
// void printVector(Vector a);
Vector vectorMultWithScalar(double a, Vector b);
Vector vectorAddition(Vector a, Vector b);
Solution explicitEuler(FunctionPtr f, double initial_value, double initial_time, double dt, int steps);
Solution RK4_single(FunctionPtr f, double initial_value, double initial_time, double dt, int steps);
Solution RK4_system(FunctionPtrVec f, Vector initial_value, double initial_time, double dt, int steps);
Solution adamsBashforth2(FunctionPtr f, std::array<double, 2>& initial_value, double initial_time, double dt, int steps);
std::pair<Eigen::MatrixXd, Eigen::VectorXd> generate_coeffMat_Rhs(double initial_radius, double final_radius, double source_strength, int num_grid_points, double delta_radius);
Eigen::VectorXd solve_linear_system_hht(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);

#endif
