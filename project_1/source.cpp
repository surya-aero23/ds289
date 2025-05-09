// source.cpp
#include "header.h"
#include <iostream>
#include "Eigen/Dense"

// vector scalar multiplication function
Vector vectorMultWithScalar(double a, Vector b) {
    int size = static_cast<int>(b.size());
    for (int i = 0; i < size; ++i) {
        b[i] *= a;
    }
    return b;
}


// vector addition function
Vector vectorAddition(Vector a, Vector b) {
    int size = static_cast<int>(b.size());
    for (int i = 0; i < size; ++i) {
        a[i] += b[i];
    }
    return a;
}


// void printVector(Vector a){
//     for (auto i : a) {
//         std::cout << i << " ";
//     }
//     std::cout << std::endl;
// }


Solution explicitEuler(FunctionPtr f, double initial_value, double initial_time, double dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);
    result.time[0] = initial_time;
    result.sol[0] = initial_value;
    double current = initial_value;
    double t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        current += dt * f(current, t);
        t += dt;
        result.time[i] = t;
        result.sol[i] = current;
    }

    return result;
}


// RK4 for single equation
Solution RK4_single(FunctionPtr f, double initial_value, double initial_time, double dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);
    result.time[0] = initial_time;
    result.sol[0] = initial_value;
    double current = initial_value;
    double t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        double k1 = dt * f(current, t);
        double k2 = dt * f(current + 0.5 * k1, t + 0.5 * dt);
        double k3 = dt * f(current + 0.5 * k2, t + 0.5 * dt);
        double k4 = dt * f(current + k3, t + dt);

        current += (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        t += dt;
        result.time[i] = t;
        result.sol[i] = current;
    }

    return result;
}


Solution adamsBashforth2(FunctionPtr f, std::array<double, 2>& initial_value, double initial_time, double dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);

    result.time[0] = initial_time;
    result.sol[0] = initial_value[0];
    result.time[1] = dt;
    result.sol[1] = initial_value[1];

    double t = initial_time;
    double current = initial_value[0];
    double next = initial_value[1];

    for (int i = 2; i <= steps; ++i) {
        double temp = next + (3.0 / 2.0) * dt * f(next, t + dt) - (1.0 / 2.0) * dt * f(current, t);
        current = next;
        next = temp;
        t += dt;
        result.time[i] = t;
        result.sol[i] = next;
    }
    
    return result;
}


// RK4 for system of equations
Solution RK4_system(FunctionPtrVec f, Vector initial_value, double initial_time, double dt, int steps) {
    int size = static_cast<int>(initial_value.size());
    Solution result;
    result.time.resize(steps + 1);
    result.sol_set.resize(size);
    
    for (int i = 0; i < size; ++i) {
        result.sol_set[i].resize(steps + 1);
    }

    result.time[0] = initial_time;
    for (int i = 0; i < size; ++i) {
        result.sol_set[i][0] = initial_value[i];
    }

    Vector current = initial_value;
    double t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        Vector k1 = vectorMultWithScalar(dt, f(t, current));
        Vector k2 = vectorMultWithScalar(dt, f(t + 0.5 * dt, vectorAddition(current, vectorMultWithScalar(0.5, k1))));
        Vector k3 = vectorMultWithScalar(dt, f(t + 0.5 * dt, vectorAddition(current, vectorMultWithScalar(0.5, k2))));
        Vector k4 = vectorMultWithScalar(dt, f(t + dt, vectorAddition(current, k3)));
        // printVector(f(t, current));
        // printVector(k2);
        // printVector(k3);
        // printVector(k4);
        // std::cout<<std::endl;

        for (int j = 0; j < size; ++j) {
            current[j] += (1 / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            result.sol_set[j][i] = current[j];
        }
        t += dt;
        result.time[i] = t;
    }

    return result;
}


// Function to generate the coefficient matrix and the right-hand side vector with proper variable naming
std::pair<Eigen::MatrixXd, Eigen::VectorXd> generate_coeffMat_Rhs(double initial_radius, double final_radius, double source_strength, int num_grid_points, double delta_radius) {
    
    // Set up coefficient matrix A and right-hand side vector b
    // Eigen::MatrixXd A(num_grid_points, num_grid_points);
    // Eigen::VectorXd b(num_grid_points);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_grid_points, num_grid_points);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(num_grid_points);
    double radius;

    // Set up the coefficient matrix A (tri-diagonal matrix) and vector b
    for (int i = 0; i < num_grid_points; i++) {
        radius = initial_radius + i * delta_radius;
        if (i == 0) {
            A(i, i) = - 2 / (delta_radius * delta_radius);
            A(i, i + 1) = 2 / (delta_radius * delta_radius);
            b(i) = -source_strength;
        } else if (i == num_grid_points - 1) {
            A(i, i) = 1;
            b(i) = 1;
        } else {
            A(i, i) = - 2 / (delta_radius * delta_radius);
            A(i, i + 1) = 1 / (delta_radius * delta_radius) + 1 / (2 * radius * delta_radius);
            A(i, i - 1) = 1 / (delta_radius * delta_radius) - 1 / (2 * radius * delta_radius);
            b(i) = -source_strength;
        }
    }

    return std::make_pair(A, b);
}


// Function to solve the system of linear equations with proper variable naming
Eigen::VectorXd solve_linear_system_hht(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    return A.colPivHouseholderQr().solve(b);
}