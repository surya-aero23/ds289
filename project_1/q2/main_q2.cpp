#include "../header.h"
#include <iostream>
#include <fstream> // For file handling
#include <cmath> // For exp and pow

using namespace std;

// Define the true solution of the differential equation dy/dt = y * t^2 - 1.1 * y
double trueSolution(double t, double y0) {
    return y0 * exp((pow(t, 3) / 3) - 1.1 * t);
}

// Define the differential equation dy/dt = y * t^2 - 1.1 * y
double f(double y, double t) {
    return y * (( t * t ) - 1.1 ); 
}


int main() {
    ofstream outputFile;
        ifstream inputFile("inputs.txt");
    if (!inputFile.is_open()) {
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    // Read inputs from the inputs.txt file
    double initial_time, initial_value, final_time;
    Vector dt_values;
    string line;
    
    int index = 0;
    // read the line if @ is not present in the first character
    while (getline(inputFile, line)) {
        
        if (line.empty() || line[0] == '@') {// Skip empty lines and lines starting with '@'
            continue;
            } 
        if (index == 0){
            initial_time = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 1){
            initial_value = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 2){
            final_time = stod(line.substr(line.find(" ") + 1));
        }
        if (index >= 3){
            string dt_values_str = line.substr(line.find(" ") + 1);
            size_t pos = 0;
            while ((pos = dt_values_str.find(" ")) != string::npos) {
                dt_values.push_back(stod(dt_values_str.substr(0, pos)));
                dt_values_str.erase(0, pos + 1);
            }
            dt_values.push_back(stod(dt_values_str)); // Add the last dt value
        }
        index++;
    }

    int steps;
    Solution explicitEulerSolution, RK4Solution, adamsBashforth2Solution;

    // Open the output file
    outputFile.open("output.csv");
    if (!outputFile.is_open()){
    	cerr << "Unable to open file." << endl;
    	return 1;
    	}
    
    outputFile << "Time,True Solution,Explicit Euler,RK4,Adams-Bashforth 2nd order,Delta t" << endl;

    // Loop over different time step sizes
    for (double dt : dt_values) {
        steps = static_cast<int>((final_time - initial_time) / dt); // Calculate the number of steps based on the time step size
        array<double, 2> adamsBashforth2_ic = {initial_value, trueSolution(dt, initial_value)}; // Calculate the initial condition for AB2 based on true solution

        explicitEulerSolution = explicitEuler(f, initial_value, initial_time, dt, steps);
        RK4Solution = RK4_single(f, initial_value, initial_time, dt, steps);
        adamsBashforth2Solution = adamsBashforth2(f, adamsBashforth2_ic, initial_time, dt, steps);

        // Write results for each time step to the output file
        for (int i = 0; i <= steps; ++i) {
            double t = i * dt;
            double true_sol = trueSolution(t, initial_value);

            outputFile << t << "," << true_sol << "," << explicitEulerSolution.sol[i] << ",";
            outputFile << RK4Solution.sol[i] << ",";
            outputFile << adamsBashforth2Solution.sol[i];
            if (i == 0)
                outputFile << "," << dt << endl;
            else
                outputFile << endl;
        }
        // outputFile << endl;
    }

    // Close the output file
    outputFile.close();

    cout << "Results written to output.csv" << endl;

    return 0;
}

