// source.cpp
#include "header.h"
#include <iostream>

// vector scalar multiplication function
Vector vectorMultWithScalar(real a, Vector b) {
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


Solution explicitEuler(FunctionPtr f, real initial_value, real initial_time, real dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);
    result.time[0] = initial_time;
    result.sol[0] = initial_value;
    real current = initial_value;
    real t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        current += dt * f(current, t);
        t += dt;
        result.time[i] = t;
        result.sol[i] = current;
    }

    return result;
}


// RK4 for single equation
Solution RK4_single(FunctionPtr f, real initial_value, real initial_time, real dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);
    result.time[0] = initial_time;
    result.sol[0] = initial_value;
    real current = initial_value;
    real t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        real k1 = dt * f(current, t);
        real k2 = dt * f(current + 0.5 * k1, t + 0.5 * dt);
        real k3 = dt * f(current + 0.5 * k2, t + 0.5 * dt);
        real k4 = dt * f(current + k3, t + dt);

        current += (1 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        t += dt;
        result.time[i] = t;
        result.sol[i] = current;
    }

    return result;
}


Solution adamsBashforth2(FunctionPtr f, std::array<real, 2>& initial_value, real initial_time, real dt, int steps) {
    Solution result;
    result.time.resize(steps + 1);
    result.sol.resize(steps + 1);

    result.time[0] = initial_time;
    result.sol[0] = initial_value[0];
    result.time[1] = dt;
    result.sol[1] = initial_value[1];

    real t = initial_time;
    real current = initial_value[0];
    real next = initial_value[1];

    for (int i = 2; i <= steps; ++i) {
        real temp = next + (3.0 / 2.0) * dt * f(next, t + dt) - (1.0 / 2.0) * dt * f(current, t);
        current = next;
        next = temp;
        t += dt;
        result.time[i] = t;
        result.sol[i] = next;
    }
    
    return result;
}


// RK4 for system of equations
Solution RK4_system(FunctionPtrVec f, Vector initial_value, real initial_time, real dt, int steps) {
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
    real t = initial_time;

    for (int i = 1; i <= steps; ++i) {
        Vector k1 = vectorMultWithScalar(dt, f(t, current));
        Vector k2 = vectorMultWithScalar(dt, f(t + 0.5 * dt, vectorAddition(current, vectorMultWithScalar(0.5, k1))));
        Vector k3 = vectorMultWithScalar(dt, f(t + 0.5 * dt, vectorAddition(current, vectorMultWithScalar(0.5, k2))));
        Vector k4 = vectorMultWithScalar(dt, f(t + dt, vectorAddition(current, k3)));

        for (int j = 0; j < size; ++j) {
            current[j] += (1 / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            result.sol_set[j][i] = current[j];
        }
        t += dt;
        result.time[i] = t;
    }

    return result;
}


// 2nd order central difference for first derivative (CD2_first) at i
real CD2_first(std::array<real, 3> f, real dx) {
    // {f_i-1, f_i, f_i+1}
    return (f[2] - f[0]) / (2 * dx);
}

// 2nd order central difference for second derivative (CD2_second) at i
real CD2_second(std::array<real, 3> f, real dx) {
    // {f_i-1, f_i, f_i+1}
    return (f[2] - 2 * f[1] + f[0]) / (dx * dx);
}


// 4th order central difference for second derivative (CD4_second) at i 
real CD4_second(std::array<real, 5> f, real dx) {
    // {f_i-2, f_i-1, f_i, f_i+1, f_i+2}
    return (-f[4] + 16 * f[3] - 30 * f[2] + 16 * f[1] - f[0]) / (12 * dx * dx);
}


// 6th order central difference for second derivative (CD6_second) at i
real CD6_second(std::array<real, 7> f, real dx) {
    // {f_i-3, f_i-2, f_i-1, f_i, f_i+1, f_i+2, f_i+3}
    return (2 * f[6] - 27 * f[5] + 270 * f[4] - 490 * f[3] + 270 * f[2] - 27 * f[1] + 2 * f[0]) / (180 * dx * dx);
}


// explicit Euler function for diffusion equation from initial to final time
Solution explicitEuler_Diffusion(FunctionPtrSingle initCondition, real initial_time, real final_time, real Nx, real rd, real alpha, int cdChoice, real domain_length, int checkPoint, std::string filename, int setPrecision) {
    
    // Set the simulation parameters
    real dx = domain_length / (Nx - 1);
    real dt = rd * dx * dx / alpha;
    int Nt = static_cast<int>((final_time - initial_time) / dt);
    
    // Facilitate checkpointing
    std :: cout << "Nt: " << Nt << std :: endl;
    if (Nt < checkPoint) {
        checkPoint = Nt;
    }
    
    // Initialize the solution variables
    Solution result;
    // result.time.resize(t_size);
    result.space.resize(Nx);
    result.sol.resize(Nx + 1);
    result.sol_set.resize(checkPoint);
    
    // Set the initial condition
    result.sol[0] = initial_time; // zeroth index of result.sol is time
    for (int i = 1; i < Nx+1 ; i++) {
        // Set the grid locations
        result.space[i-1] = (i-1) * dx;

        // IC
        result.sol[i] = initCondition(result.space[i-1]);   
    }
    result.sol_set[0] = result.sol;

    // Open the file to write the heading
    std :: ofstream outputFile(filename);
    if (!outputFile.is_open()) {    // Check if the file is open
        std :: cerr << "Unable to open " << filename << " file." << std :: endl;
        return result;
    }
    // Write the first entries of the columns for time and space
    outputFile << "time,";
    for (int i = 0; i < Nx; ++i) {
        outputFile << std :: setprecision(setPrecision) << result.space[i] << ",";
    }
    outputFile << std :: endl;
    outputFile.close();
    
    // Set the values for time looping
    real t = initial_time + dt;

    // 't_index' for time and 'count' for correction of final time step
    int t_index = 1;
    int sol_set_index = 1;
    int count = 0; 

    // To hold the temporary solution to avoid time level clashes
    Vector temp_sol(Nx - 2);
    
    while (t <= final_time) { // time loop
        // time values 
        result.sol[0] = t;
        
        // To avoid calculating this at every space point
        real alpha_dt = alpha * dt;
        
        if (cdChoice == 2) {    // 2nd order central difference
            // space loop
            for (int j = 2; j < Nx; ++j) {
                temp_sol[j - 2] = result.sol[j] + alpha_dt * CD2_second({result.sol[j - 1], result.sol[j], result.sol[j + 1]}, dx);
                }
            
            // Update the solution from temp
            for (int k = 2; k < Nx; ++k) {
                result.sol[k] = temp_sol[k - 2];
            }
        }

        else if (cdChoice == 4) {   // 4th order central difference
            // Using the periodic boundary condition
            temp_sol[0] = result.sol[2] + alpha_dt * CD4_second({result.sol[Nx - 1], result.sol[1], result.sol[2], result.sol[3], result.sol[4]}, dx);
            
            // space loop
            for (int j = 3; j < Nx - 1; ++j) { 
                temp_sol[j - 2] = result.sol[j] + alpha_dt * CD4_second({result.sol[j - 2], result.sol[j - 1], result.sol[j], result.sol[j + 1], result.sol[j + 2]}, dx);    
                } 

            // Using the periodic boundary condition
            temp_sol[Nx - 3] = result.sol[Nx - 1] + alpha_dt * CD4_second({result.sol[Nx - 3], result.sol[Nx - 2], result.sol[Nx - 1], result.sol[Nx], result.sol[2]}, dx);
            
            // Update the solution from temp
            for (int k = 2; k < Nx; ++k) { 
                result.sol[k] = temp_sol[k - 2];
            }
        }

        else if (cdChoice == 6) {   // 6th order central difference
            // Using the periodic boundary condition
            temp_sol[0] = result.sol[2] + alpha_dt * CD6_second({result.sol[Nx - 2], result.sol[Nx - 1], result.sol[1], result.sol[2], result.sol[3], result.sol[4], result.sol[5]}, dx);
            temp_sol[1] = result.sol[3] + alpha_dt * CD6_second({result.sol[Nx - 1], result.sol[1], result.sol[2], result.sol[3], result.sol[4], result.sol[5], result.sol[6]}, dx);
            
            // space loop
            for (int j = 4; j < Nx - 2; ++j) { 
                temp_sol[j - 2] = result.sol[j] + alpha_dt * CD6_second({result.sol[j - 3], result.sol[j - 2], result.sol[j - 1], result.sol[j], result.sol[j + 1], result.sol[j + 2], result.sol[j + 3]}, dx);
            }

            // Using the periodic boundary condition
            temp_sol[Nx-4] = result.sol[Nx-2] + alpha_dt * CD6_second({result.sol[Nx-5], result.sol[Nx-4], result.sol[Nx-3], result.sol[Nx-2], result.sol[Nx-1], result.sol[1], result.sol[2]}, dx);
            temp_sol[Nx-3] = result.sol[Nx-1] + alpha_dt * CD6_second({result.sol[Nx-4], result.sol[Nx-3], result.sol[Nx-2], result.sol[Nx-1], result.sol[1], result.sol[2], result.sol[3]}, dx);
            
            // Update the solution from temp
            for (int k = 2; k < Nx; ++k) {
                result.sol[k] = temp_sol[k - 2];
            }
        }
    
        else {
            std::cout << "Invalid choice of central difference scheme" << std::endl;
            break;
        }
        
        // Set the periodic boundary condition (no need to set it calculations at the boundary are avoided)
        // result.sol[1] = result.sol[Nx];

        // Copy the solution to the solution set
        result.sol_set[sol_set_index] = result.sol;
        sol_set_index++;

        // Correction for the final time step if dt does not divide (final_time - initial_time)
        if (t + dt > final_time) { 
            dt = final_time - t;
            
            // // To avoid infinite loop
            count++; 
            if (count > 1) {
                break;
            }  
        }

        // Update the time and time index
        t += dt;
        t_index++;

        // Write the sol_set with appropriate time values to the csv file with the filename after every checkpoint
        if (t_index % checkPoint == 0) {
            std :: ofstream outputFile;
            outputFile.open(filename, std :: ofstream :: app);

            for (int j = 0; j < checkPoint ; j++) {
                for (int i = 0; i < Nx + 1; ++i) {
                    // write in scientific notation with 8 significant digits 
                    outputFile << std :: scientific << std :: setprecision(setPrecision) << result.sol_set[j][i] << ",";    
                }
                outputFile << std :: endl;
            }

            // Clear the solution set
            sol_set_index = 0;

            // Close the file
            outputFile.close();
        }
    }

    // Write the left over solution set to the csv file
    if (sol_set_index != 0) { 
        std :: ofstream outputFile;
        outputFile.open(filename, std :: ofstream :: app);

        for (int j = 0; j < sol_set_index ; j++) {
            for (int i = 0; i < Nx + 1; ++i) {
                // write in scientific notation with 8 significant digits
                outputFile << std :: scientific << std :: setprecision(setPrecision) << result.sol_set[j][i] << ",";    
            }
            outputFile << std :: endl;
        }
        outputFile.close();        
    }

    return result;
}



// To obtain the index of the matrix given the x and y coordinates for 2d laplacian problem
int find_index(int x, int y, int Nx)
{
    return x + (Nx + 2) * y;
}


// Generate pentadiagonal matrix for the 2d laplacian equation
void solve_2dlaplacian(int Nx, int Ny, double dx, double dy) {

    // Define the coefficients
    real a = 1 / (dx * dx);
    real b = 1 / (dy * dy);
    real c = -2 * (a + b);

    
    // Define the grid (excluding the boundary points)
    Eigen::SparseMatrix<double> A((Nx + 2) * (Ny + 2), (Nx + 2) * (Ny + 2)); // coefficient matrix A
    Eigen::VectorXd B((Nx + 2) * (Ny + 2));                                              // RHS vector b
    Eigen::VectorXd T((Nx + 2) * (Ny + 2));                                              // solution vector T

    // Initialize the matrix and the right-hand side
    T.setZero();
    A.setZero();
    B.setZero();

    // Set the boundary conditions
    for (int x = 0; x < Nx + 2; x++)
    {
        std::cout << "---Populating row x:" << x << std::endl;   
        for (int y = 0; y < Ny + 2; y++)
        {
            // Get the indices of the surrounding points    
            int center = find_index(x, y, Nx);
            int left = find_index(x - 1, y, Nx);
            int right = find_index(x + 1, y, Nx);
            int up = find_index(x, y - 1, Nx);
            int down = find_index(x, y + 1, Nx);

            //Left surface
            if (x == 0)
            {
                A.insert(center, center) = 1;
                B[center] = 30;
            }

            //Right surface
            else if (x == Nx + 1)
            {
                A.insert(center, center) = 1;
                B[center] = 60;
            }
            
            //Top surface
            else if (y == 0)
            {
                A.insert(center, center) = c;
                A.insert(center, down) = 2 * b;
                A.insert(center, right) = a;
                A.insert(center, left) = a;
                B[center] = 0;
            }

            //Bottom surface
            else if (y == Ny + 1)
            {
                
                A.insert(center, center) = c + 2/dy;
                A.insert(center, up) = 2 * b;
                A.insert(center, left) = a;
                A.insert(center, right) = a;
                B[center] = 120/dy;
            }

            //Interior points
            else
            {
                A.insert(center, center) = c;
                A.insert(center, up) = b;
                A.insert(center, down) = b;
                A.insert(center, left) = a;
                A.insert(center, right) = a;
                B[center] = 0;
            }
        
        }
    
    }   

    //Write the matrix A to a csv file under
    std :: cout << "\n---Writing the matrix A to a csv file" << std :: endl;
    std::ofstream file;
    file.open("A.csv");
    file << "x,y,value" << std::endl;
    for (int i = 0; i < A.outerSize(); i++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
        {
            file << it.row() << "," << it.col() << "," << it.value() << std::endl;
        }
    }
    file.close();
    
    // Solve the system using sparse LU decomposition
    std :: cout << "\n---Solving the system using sparse LU decomposition" << std :: endl;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_LU;
    sparse_LU.analyzePattern(A);
    sparse_LU.factorize(A);
    T = sparse_LU.solve(B);

    //Output the solution
    std :: cout << "\n---Outputting the solution to Temperature.csv" << std :: endl;
    file.open("Temperature.csv");
    for (int y = 0; y < Ny + 2; y++)
    {
        for (int x = 0; x < Nx + 2; x++)
        {
            if (x != 0) {
                file << ",";
            }
            file << std::scientific << std::setprecision(20) << T[find_index(x, y, Nx)];
        }
        file << std::endl;
    }
    file.close();

}


// Build matrix A for implicit Euler CD2 diffusion equation
void buildMatrix_ImpEuler_Diffusion(std::vector<Vector>& A, int N, real r_d) {
    A.resize(N, std::vector<real>(N, 0.0));
    A[0][0] = 1 + 2 * r_d;
    A[0][1] = -r_d;
    A[0][N - 2] = -r_d;

    A[N - 1][N - 1] = 1 + 2 * r_d;
    A[N - 1][N - 2] = -r_d;
    A[N - 1][1] = -r_d;

    for (int i = 1; i < N - 1; i++) {
        A[i][i] = 1 + 2 * r_d;
        A[i][i - 1] = -r_d;
        A[i][i + 1] = -r_d;
    }
}


// Jacobi method to solve the system of equations
void Jacobi_Solver(const std::vector<Vector>& A, const std::vector<real>& b, std::vector<real>& x, int N, int max_iter, real tolerance) {
    std::vector<real> x_prev(N, 0.0);
    int iter = 0;
    real error = tolerance + 1.0;

    while (iter < max_iter && error > tolerance) {
        x_prev = x;

        for (int i = 0; i < N; i++) {
            real sum = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    sum += A[i][j] * x_prev[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        error = 0.0;
        for (int i = 0; i < N; i++) {
            error = std::max(error, std::abs(x[i] - x_prev[i]));
        }

        iter++;
    }
}


// Use jacobi method to solve the system of equations using implicit Euler and second order CDS
void solveDiffusion_ImpEuler(Vector& u, Vector& u_new, int N, real r_d, real dx, real diffusion_coefficient, real initial_time, real final_time, int maxIterJacobi, real toleranceJacobi) {
    real dt = r_d * dx * dx / diffusion_coefficient;
    int nsteps = std::floor((final_time - initial_time) / dt);
    real time = initial_time;

    std::vector<Vector> A;
    buildMatrix_ImpEuler_Diffusion(A, N, r_d);

    for (int i = 0; i < nsteps; i++) {
        Jacobi_Solver(A, u, u_new, N, maxIterJacobi, toleranceJacobi);
        u = u_new;
        u_new = Vector(N, 0.0);
        time += dt;
    }

    if (final_time - time > 0) {
        dt = final_time - time;
        nsteps = 1;
        Jacobi_Solver(A, u, u_new, N, maxIterJacobi, toleranceJacobi);
    }
}