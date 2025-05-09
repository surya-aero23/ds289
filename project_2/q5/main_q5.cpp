#include "../header.h"
#include <iostream>
#include <fstream> // For file handling
#include <cmath>   // For the exp function

using namespace std;

// Initial condition function
double initial_condition(double x) {
    return sin(x) + sin(4 * x);
}


int main() {
    cout << "Reading inputs from input.txt file." << endl;
    // Open the inputs.txt file
    ifstream inputFile("input.txt");
    
    if (!inputFile.is_open()) { // Check if the file is open    
        cerr << "Unable to open input.txt file." << endl;
        return 1;
    }

    double initial_time, final_time, domain_length, diffusion_coefficient;
    double r_d, toleranceJacobi;
    int maxIterJacobi;
    vector<int> N_values;
    
    // Read the variables from the inputs.txt file
    string line;
    int index = 0;  // Index to keep track of the variable being read
    
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '@') { // Skip empty lines and comments  
            continue;
        }

        if (index == 0) {
            initial_time = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 1) {
            final_time = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 2) {
            domain_length = stod(line.substr(line.find(" ") + 1));
            if (domain_length == -1) {
                domain_length = 2 * M_PI;
            }
        }

        if (index == 3) {
            r_d = stod(line.substr(line.find(" ") + 1));
        }
        
        if (index == 4) {
            diffusion_coefficient = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 5) {
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                N_values.push_back(stoi(line));
                getline(inputFile, line);
            }
        }

        if (index == 6) {
            maxIterJacobi = stoi(line.substr(line.find(" ") + 1));
        }

        if (index == 7) {
            toleranceJacobi = stod(line.substr(line.find(" ") + 1));
        }

        // Increment the index  
        index++;
    }

    // Close the file
    inputFile.close();

    // print the variables
    int num_N_values = static_cast<int>(N_values.size());

    cout << "Initial time: " << initial_time << endl;
    cout << "Final time: " << final_time << endl;
    cout << "Domain length: " << domain_length << endl;
    cout << "Stability parameter: " << r_d << endl;
    cout << "Diffusion coefficient: " << diffusion_coefficient << endl;
    cout << "Grid sizes: ";
    for (int i = 0; i < num_N_values; i++) {
        cout << N_values[i] << ", ";
    }
    cout << endl;
    cout << "Maximum iterations (Jacobi method): " << maxIterJacobi << endl;
    cout << "Tolerance (Jacobi method): " << toleranceJacobi << endl;


    // Solve the diffusion equation for different grid sizes
    for (int i = 0; i < num_N_values; i++) { 

        // Define the grid size
        int N = N_values[i]; 
        cout << "\n\n---Solving the diffusion equation for N = " << N << "---" << endl;  

        // Define the solution parameters
        double dx = domain_length / (N - 1);

        // Initialize the solution vectors
        vector<double> u(N, 0.0);
        vector<double> u_new(N, 0.0);
        vector<double> x;

        // Initialize the x and u vector
        cout << "---Initializing the solution vectors ---" << endl;
        for (int i = 0; i < N; i++) {
            x.push_back(i * dx);
            u[i] = initial_condition(x[i]);
        }

        // Solve the diffusion equation
        cout << "---Solving the diffusion equation using Implicit Euler method---" << endl;
        solveDiffusion_ImpEuler(u, u_new, N, r_d, dx, diffusion_coefficient, initial_time, final_time, maxIterJacobi, toleranceJacobi);

        // Write the solution to a csv file with the format: x, u
        cout << "---Writing the solution to a csv file---" << endl;
        ofstream outputFile("N"+ to_string(N) + "_output.csv");
        if (!outputFile.is_open()) {
            cerr << "Unable to open output.csv file." << endl;
            return 1;
        }

        for (int i = 0; i < N; i++) {
            outputFile << scientific << setprecision(20) << x[i] << ", " << scientific << setprecision(20) << u_new[i] << endl;
        }

        outputFile.close();
    
    }
      
    return 0;
}
