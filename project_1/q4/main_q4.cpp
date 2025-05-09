#include "../header.h"
#include <iostream>
#include <fstream> // For file handling
#include <cmath>   // For the exp function

using namespace std;

// Define the differential equation dy/dt = -2 * 10 ^ 5 (y - e ^ (-x)) - e ^ (-x)
double f(double y, double t) {
    return (-2 * pow(10, 5) * y) - (2 * pow(10, 5) * exp(-t)) - exp(-t); 
}

// Analytical solution
double analyticalSolution(double t) {
    return exp(-t) - exp(-2 * pow(10, 5) * t);
}


Solution implicitEuler(FunctionPtr f, double initial_value, double initial_time, double dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);
    result.time[0] = initial_time;
    result.sol[0] = initial_value;

    double current = initial_value;
    double t = initial_time;
    double k = 2 * pow(10, 5);

    for (int i = 1; i <= steps; ++i) {
        t += dt;
        current = ( current + (dt * exp(-t) * (k - 1)) ) / (1 + (dt * k));
        result.time[i] = t;
        result.sol[i] = current;
    }
    return result;
}

int main() {
    double initial_value, initial_time, final_time, dt;
    ifstream inputFile("inputs.txt");
    if (!inputFile.is_open()) {
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    string line;
    int index = 0;
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '@') {
            continue;
        }
        if (index == 0) {
            initial_value = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 1) {
            initial_time = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 2) {
            final_time = stod(line.substr(line.find(" ") + 1));
        }
        if (index == 3) {
            dt = stod(line.substr(line.find(" ") + 1));
        }
        // cout<<line<<endl;
        index++;
    }
    inputFile.close();

    int steps = static_cast<int>((final_time - initial_time) / dt);

    Solution solution = implicitEuler(f, initial_value, initial_time, dt, steps);

    // Write the output to a CSV file
    ofstream outputFile("output_q4.csv");
    if (outputFile.is_open()) {
        outputFile << "Time,True Solution,Implicit Euler" << endl;
        for (size_t i = 0; i < solution.time.size(); ++i) {
            outputFile << solution.time[i] << "," << analyticalSolution(solution.time[i]) << "," << solution.sol[i] << endl;
        }
        outputFile.close();
        cout << "Output written to output_q4.csv" << endl;
    } else {
        cerr << "Unable to open file." << endl;
        return 1;
    }

    return 0;
}
