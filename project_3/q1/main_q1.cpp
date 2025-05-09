#include "../header.h"
#include <iostream>
#include <cmath>
#include <iomanip>


using namespace std;


double initCondition(double x) {
    return sin(4 * M_PI * x) +  sin(6 * M_PI * x) +  sin(10 * M_PI * x);
}


int main() {
    
    cout << "\n---Reading inputs from inputs.txt file---\n" << endl;
    // Open the inputs.txt file
    ifstream inputFile("input.txt");
    
    if (!inputFile.is_open()) { // Check if the file is open    
        cerr << "Unable to open inputs.txt file." << endl;
        return 1;
    }

    double t_initial, t_final, delta_t, diffusion_coefficient, domain_length;
    vector<int> N, schemes;
    
    // Read the variables from the inputs.txt file
    string line;
    int index = 0;  // Index to keep track of the variable being read
    
    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '@') { // Skip empty lines and comments  
            continue;
        }

        if (index == 0) {
            t_initial = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 1) {
            t_final = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 2) {
            delta_t = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 3) {
            domain_length = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 6) {
            diffusion_coefficient = stod(line.substr(line.find(" ") + 1));
        }

        if (index == 5) {
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                schemes.push_back(stoi(line));
                getline(inputFile, line);
            }
        }

        if (index == 4) {
            while (true) {
                if (line.empty() || line[0] == '@') {
                    break;
                }
                N.push_back(stoi(line));
                getline(inputFile, line);
            }
        }

        // Increment the index  
        index++;
    }
    
    // Close the file
    inputFile.close();




    // print the parameters
    cout << "t_initial: " << t_initial << endl;
    cout << "t_final: " << t_final << endl;
    cout << "delta_t: " << delta_t << endl;
    cout << "domain_length: " << domain_length << endl;
    cout << "diffusion_coefficient: " << diffusion_coefficient << endl;
    
    cout << "scheme: ";
    int scheme_size = static_cast<int>(schemes.size());
    for (int i = 0; i < scheme_size; i++) {
        cout << schemes[i] << " ";
    }
    cout << endl;
    
    cout << "N: ";
    int N_size = static_cast<int>(N.size());
    for (int i = 0; i < N_size; i++) {
        cout << N[i] << " ";
    }
    cout << "\n" << endl;




    // for each scheme, calculate the solution
    for (int s = 0; s < scheme_size; s++) {
        
        double scheme = schemes[s];
        cout << "\n---Calculating solution for scheme = " << scheme << "---\n" << endl;

        // for each N, calculate the solution
        for (int i = 0; i < N_size; i++) {
            cout << "\t---Calculating solution for N = " << N[i] << "---" << endl;
            
            Vector u(N[i] + 1), u_new(N[i] + 1), x(N[i]);
            double dx = domain_length / (N[i] - 1);
            
            // Calculate the initial condition
            u[0] = t_initial;
            for (int j = 1; j < N[i] + 1; j++) {
                u[j] = initCondition((j - 1) * dx);
                x[j - 1] = (j - 1) * dx;
                }
            
            // Calculate the solution using the upwind Euler scheme
            if (scheme == 1) 
            {
                // open a file to write the solution (name it based on N[i] and scheme)
                string filename = "output_N_" + to_string(N[i]) + "_upwind_scheme" + ".txt";
                ofstream outputFile(filename);

                if (!outputFile.is_open()) { // Check if the file is open    
                    cerr << "Unable to open " << filename << " file." << endl;
                    return 1;
                }

                // set file precision to 10 and notation to scientific
                outputFile << scientific << setprecision(10);
                outputFile << "time" << ", ";
                for (int j = 0; j < N[i]; j++) {
                    outputFile << x[j] << ", ";
                }
                outputFile << "\n";

                // write the solution to the file
                for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                outputFile << "\n";


                // Calculate the solution
                int steps = (int) ((t_final - t_initial) / delta_t);
                double current_time = t_initial;

                for (int k = 1; k <= steps; k++) {
                    
                    upwindScheme(u, u_new, delta_t, dx, N[i]);
                    u = u_new;
                    current_time += delta_t;

                    // write the solution to the file
                    for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                    outputFile << "\n";
                
                }

                // Check if the current time is the final time
                if (current_time < t_final) {
                    upwindScheme(u, u_new, t_final - current_time, dx, N[i]);
                    u = u_new;

                    // write the solution to the file
                    for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                
                }
                outputFile.close();

                cout << "\tSolution written to " << filename << " file.\n" << endl;

            }


            // Calculate the solution using the CD2 scheme with Euler (CD2_Euler_ViscousBurgers)
            if (scheme == 2) {
                // open a file to write the solution (name it based on N[i] and scheme)
                string filename = "output_N_" + to_string(N[i]) + "_CD2_scheme" + ".txt";
                ofstream outputFile(filename);

                if (!outputFile.is_open()) { // Check if the file is open    
                    cerr << "Unable to open " << filename << " file." << endl;
                    return 1;
                }

                // set file precision to 10 and notation to scientific
                outputFile << scientific << setprecision(10);
                outputFile << "time" << ", ";
                for (int j = 0; j < N[i]; j++) {
                    outputFile << x[j] << ", ";
                }
                outputFile << "\n";

                // write the solution to the file
                for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                outputFile << "\n";

                // Calculate the solution
                int steps = (int) ((t_final - t_initial) / delta_t);
                double current_time = t_initial;

                for (int k = 1; k <= steps; k++) {
                    CD2_Euler_ViscousBurgers(u, u_new, diffusion_coefficient, delta_t, dx, N[i]);
                    u = u_new;
                    current_time += delta_t;

                    // write the solution to the file
                    for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                    outputFile << "\n";
                
                }

                // Check if the current time is the final time
                if (current_time < t_final) {
                    CD2_Euler_ViscousBurgers(u, u_new, diffusion_coefficient, t_final - current_time, dx, N[i]);
                    u = u_new;

                    // write the solution to the file
                    for (int j = 0; j <= N[i]; j++) {
                        outputFile << u[j] << ", ";
                    }
                
                }
                outputFile.close();

                cout << "\tSolution written to " << filename << " file.\n" << endl;


            }
        }
    }
    
    cout << "\n---End of program---" << endl;
    return 0;
}
