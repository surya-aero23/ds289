#include "../header.h"
#include <iostream>
#include <cmath>   // For the exp function
#include <sys/stat.h>  // For checking if a directory exists


using namespace std;


real initCondition(real x) {
    return sin(x) + sin(4*x);
}


int main() {
    cout << "Reading inputs from inputs.txt file." << endl;
    // Open the inputs.txt file
    ifstream inputFile("input.txt");
    
    if (!inputFile.is_open()) { // Check if the file is open    
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    real initial_time, final_time, domain_length, diffusion_coefficient;
    int checkPoint, setPrecision;
    Vector rd_values, cdChoices;
    vector <int> grid_sizes;
    
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
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                grid_sizes.push_back(stoi(line));
                getline(inputFile, line);
            }
        }

        if (index == 4) {
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                rd_values.push_back(stod(line));
                getline(inputFile, line);
            }
        }
        
        if (index == 5) {
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                cdChoices.push_back(stod(line));
                getline(inputFile, line);
            }
        }

        if (index == 6) {
            diffusion_coefficient = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 7) {
            checkPoint = stoi(line.substr(line.find(" ") + 1));
        }

        if (index == 8) {
            setPrecision = stoi(line.substr(line.find(" ") + 1));
        }

        // Increment the index  
        index++;
    }
    
    // Close the file
    inputFile.close();

    // Print the inputs
    cout << "Initial time: " << initial_time << endl;
    cout << "Final time: " << final_time << endl;
    cout << "Domain length: " << domain_length << endl;
    cout << "Grid sizes: ";
    for (auto& elem: grid_sizes) {
        cout << elem << " ";
    }
    cout << endl;
    cout << "rd values: ";
    for (auto& elem: rd_values) {
        cout << elem << " ";
    }
    cout << endl;
    cout << "cdChoices: ";
    for (auto& elem: cdChoices) {
        cout << elem << " ";
    }
    cout << endl;
    cout << "Diffusion coefficient: " << diffusion_coefficient << endl;
    cout << "Check point: " << checkPoint << endl;
    cout << "Set precision: " << setPrecision << endl;
    cout << endl;

    // check if outputs directory exists
    struct stat info;
        if (stat("outputs", &info) != 0) {
            cout << "Creating outputs directory." << endl;
            system("mkdir outputs");
        }
    
    // Begin the simulation

    // for all values of grid sizes, rd values and cdChoices, solve the diffusion equation
    int num_grid_sizes = static_cast <int> (grid_sizes.size());
    int num_rd_values = static_cast <int> (rd_values.size());
    int num_cdChoices = static_cast <int> (cdChoices.size());

    for (int cd_index = 0; cd_index < num_cdChoices; cd_index++) {  // Loop over all the values of cdChoices
        int cdChoice = cdChoices[cd_index];
        
        for (int nx_index = 0; nx_index < num_grid_sizes; nx_index++) { // Loop over all the values of grid sizes
        int Nx = grid_sizes[nx_index];
        
            for (int rd_index = 0; rd_index < num_rd_values; rd_index++) {  // Loop over all the values of rd values
                real rd = rd_values[rd_index];
                
                cout << "---Solving for CD = " <<cdChoice << " Nx = " << Nx << " rd = " << rd << endl;

                // Filename for the output
                string filename = "outputs/CD_"+ to_string(cdChoice) +"_Nx_" + to_string((int) Nx) + "_rd_" + to_string((real) rd) +  ".csv";
                
                // Solve the diffusion equation
                Solution sol = explicitEuler_Diffusion(initCondition, initial_time, final_time, Nx, rd, diffusion_coefficient, cdChoice, domain_length, checkPoint, filename, setPrecision);
                cout << "---Solution obtained for this case.\n" << endl;
                
            }
        }
    }

    return 0;
}

