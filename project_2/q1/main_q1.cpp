#include "../header.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;


int main() {
    cout << "Reading inputs from inputs.txt file." << endl;
    // Open the inputs.txt file
    ifstream inputFile("input.txt");
    
    if (!inputFile.is_open()) { // Check if the file is open    
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    int Nx, Ny, Lx, Ly;
    
    // Read the variables from the inputs.txt file
    string line;
    int index = 0;  // Index to keep track of the variable being read
    
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '@') { // Skip empty lines and comments  
            continue;
        }

        if (index == 0) {
            Lx = stoi(line.substr(line.find(" ") + 1));
        }

        if (index == 1) {
            Ly = stoi(line.substr(line.find(" ") + 1));
        }

        if (index == 2) {
            Nx = stoi(line.substr(line.find(" ") + 1));
        }

        if (index == 3) {
            Ny = stoi(line.substr(line.find(" ") + 1));
        }

        // Increment the index  
        index++;
    }
    
    // Close the file
    inputFile.close();

    double dx = (double) Lx / (double) (Nx + 1);
    double dy = (double) Ly / (double) (Ny + 1);

    // print the parameters
    cout << "Lx = " << Lx << endl;
    cout << "Ly = " << Ly << endl;
    cout << "Nx = " << Nx << endl;
    cout << "Ny = " << Ny << endl;
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;

    // Solve the 2D Laplacian
    solve_2dlaplacian(Nx, Ny, dx, dy);

    return 0;
}