The commands to use in the terminal (inside this directory) are explained in this file.

**Commands Explanation**

**1. make**

Compiles the C++ program using `g++` and the rules defined in the Makefile, creating the executable `main_q4` from the source files `main_q4.cpp` and `../source.cpp`.

**2. make plot**

Executes `main_q4` to generate numerical solutions, then runs `plot_q4.py` to create comparison plots for different time step sizes.

**3. make clean**

Removes compiled object files (`*.o`), the executable `main_q4`, output files, and plots, keeping the project directory clean.

**4. make cls**

Clears the terminal screen and deletes compiled object files, executables, output files, and plots, providing a clean workspace.