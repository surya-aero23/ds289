#ifndef HEADER_H
#define HEADER_H

#include "Eigen/Sparse"

#include <string>
#include <vector>
#include <array>
#include <iomanip>
#include <fstream> // For file handling

// Define the constant PI
#define M_PI 3.14159265358979323846

// Define the type for real numbers
typedef double real;

// Define the type for a vector of real numbers
typedef std::vector<real> Vector;

// Define the type for a function pointer that takes a real number and returns a real number
typedef real (*FunctionPtr)(real, real);
typedef real (*FunctionPtrSingle)(real);

// Define the type for a function pointer that takes a real number and a vector and returns a vector of real numbers
typedef Vector (*FunctionPtrVec)(real, Vector);

// Define the type for a solution containing time and solution values
struct Solution {
    Vector time;
    Vector space;
    Vector sol;
    std::vector<Vector> sol_set; // For system of equations
};

// Function declarations
Vector vectorMultWithScalar(real a, Vector b);
Vector vectorAddition(Vector a, Vector b);
Solution explicitEuler(FunctionPtr f, real initial_value, real initial_time, real dt, int steps);
Solution RK4_single(FunctionPtr f, real initial_value, real initial_time, real dt, int steps);
Solution RK4_system(FunctionPtrVec f, Vector initial_value, real initial_time, real dt, int steps);
Solution adamsBashforth2(FunctionPtr f, std::array<real, 2>& initial_value, real initial_time, real dt, int steps);
real CD2_first(std::array<real, 3> f, real dx);
real CD2_second(std::array<real, 3> f, real dx);
real CD4_second(std::array<real, 5> f, real dx);
real CD6_second(std::array<real, 7> f, real dx);
Solution explicitEuler_Diffusion(FunctionPtrSingle initCondition, real initial_time, real final_time, real Nx, real rd, real alpha, int cdChoice, real domain_length, int checkPoint, std::string filename, int setPrecision);
int get_index(int x, int y);
void solve_2dlaplacian(int Nx, int Ny, double dx, double dy);
void Jacobi_Solver(const std::vector<Vector>& A, const std::vector<real>& b, std::vector<real>& x, int N, int max_iter, real tolerance) ;
void solveDiffusion_ImpEuler(Vector& u, Vector& u_new, int N, real r_d, real dx, real diffusion_coefficient, real initial_time, real final_time, int maxIterJacobi, real toleranceJacobi);
void buildMatrix_ImpEuler_Diffusion(std::vector<Vector>& A, int N, real r_d);

#endif
