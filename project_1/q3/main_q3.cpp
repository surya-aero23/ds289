#include "../header.h"
#include <iostream>
#include <fstream> // For file handling
#include <cmath> // For exp and pow

using namespace std;

double mass, k, a, F_0, omega;

Vector f(double t, Vector y) {
    // y contains y1 = x and y2 = dx/dt
    Vector result(2);
    result[0] = y[1];
    result[1] = (0.5) * ( (2.5 * sin(0.5 * t)) - (6 * y[0]) - (5 * abs(y[1]) * y[1]) );
    return result;
}


int main() {
    ifstream inputFile("inputs.txt");
    if (!inputFile.is_open()) {
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    Vector initial_value; // Initial conditions:
    double initial_time, final_time;
    double dt_start, num_dt;
    Vector dt_values;
    // Read the variables from the inputs.txt file
    string line;
    int index = 0;
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '@') {
            continue;
        }
        if (index == 0) {
            mass = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 1) {
            k = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 2) {
            a = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 3) {
            F_0 = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 4) {
            omega = stod(line.substr(line.find(" ") + 1));
        }
        
        if (index == 5) {
            initial_time = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 6) {
            final_time = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 7) {
            initial_value.push_back(stod(line.substr(line.find(" ") + 1)));
        }
        if (index == 8) {
            initial_value.push_back(stod(line.substr(line.find(" ") + 1)));
        }

        if (index == 9) {
            dt_start = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 10) {
            num_dt = stod(line.substr(line.find(" ") + 1));
        }
        // cout<<line<<endl;
        index++;
    }
    inputFile.close();

    // Calculate the dt values starting from dt_start and halving each time
    dt_values.resize(num_dt);
    dt_values[0] = dt_start;
    for (size_t i = 1; i < num_dt; ++i) {
        dt_values[i] = dt_values[i - 1] / 2.0;
    }

    // Open the output file
    std::ofstream outputFile("output.csv");
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file." << std::endl;
        return 1;
    }

    // Write the header
    outputFile << "Time,Position,Velocity,Delta t" << std::endl;

    // Loop over dt values
    for (size_t k = 0; k < dt_values.size(); ++k) {
        double dt = dt_values[k];
        int steps = static_cast<int>((final_time - initial_time) / dt);
        Vector initial_value = {1.0, 0.0}; // Initial conditions: 1 m and 0 m/s

        // Calculate the solution using RK4 with current dt
        Solution solution = RK4_system(f, initial_value, initial_time, dt, steps);

        // Write the data for the current dt to the output file
        for (int i = 0; i <= steps; ++i) {
            outputFile << solution.time[i] << ",";
            for (size_t j = 0; j < initial_value.size(); ++j) {
                outputFile << solution.sol_set[j][i];
                if (j < initial_value.size() - 1) {
                    outputFile << ",";
                }
            
            }
            if (i == 0) {
                outputFile << "," << dt;
            }
            outputFile << std::endl;
        }
        // Separate results for different dt values with an empty line
        // outputFile << std::endl;
    }

    // Close the output file
    outputFile.close();

    std::cout << "Results written to output.csv" << std::endl;

    return 0;
}