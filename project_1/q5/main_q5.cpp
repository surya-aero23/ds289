// main.cpp for Q5
#include <iostream>
#include "../header.h"
#include <bits/stdc++.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// Main function
int main() {
    // Read inputs from the inputs.txt file
    double r_initial, r_final, grid_size, dT_dr, T_final, delta_r;
    Vector source;
    
    ifstream inputFile("inputs.txt");
    if (!inputFile.is_open()) {
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    string line;
    int index = 0;
    while (getline(inputFile, line)) {
        
        if (line.empty() || line[0] == '@') {// Skip empty lines and lines starting with '@'
            continue;
            } 
        if (index == 0){
            r_initial = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 1){
            r_final = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 2){
            grid_size = stod(line.substr(line.find(" ") + 1));
            delta_r = (r_final - r_initial) / (grid_size - 1);
        }
        if (index >= 3){
            string source_str = line.substr(line.find(" ") + 1);
            size_t pos = 0;
            while ((pos = source_str.find(" ")) != string::npos) {
                source.push_back(stod(source_str.substr(0, pos)));
                source_str.erase(0, pos + 1);
            }
            source.push_back(stod(source_str)); // Add the last source value
        }
        
        index++;
    }
    inputFile.close();

    if (grid_size == 0) {
        std::cerr << "Error: No grid size provided in the input file." << std::endl;
        return 1;
    }


for (double S : source) {
    std::cout << "Grid size (N) = " << grid_size << std::endl;
    std::cout<< "S = " << S << std::endl;

    // Calculating the step size
    double delta_r = (r_final - r_initial) / (grid_size - 1);

    // Generate the coefficient matrix and the right-hand side vector
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> result = generate_coeffMat_Rhs(r_initial, r_final, S, grid_size, delta_r);
    Eigen::MatrixXd A = result.first;
    Eigen::VectorXd b = result.second;

    // Solve the system of linear equations
    Eigen::VectorXd solution = solve_linear_system_hht(A, b);

    // Write the results to a CSV file
    std::vector<double> r_values;
    for (size_t i = 0; i < solution.size(); ++i) {
        r_values.push_back(r_initial + i * delta_r);
    }

    std::ofstream outputFile("output_q5_S" + std::to_string((int) S) + ".csv");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return 1;
    }

    outputFile << "r,T" << std::endl;
    for (size_t i = 0; i < solution.size(); ++i) {
        outputFile << r_values[i] << ", " << solution[i] << std::endl;
    }

    outputFile.close();

    std::cout << "Results written to output_q5_S" << S << ".csv" << std::endl;  

}

    return 0;
}