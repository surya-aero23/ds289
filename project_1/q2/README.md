The commands to use in the terminal (inside this directory) are explained in this file.

**Commands Explanation**

**1. make**

This command invokes the Makefile to compile the C++ program using the GNU Compiler Collection (`g++`). It uses the rules defined in the Makefile to compile source files (`main_q2.cpp` and `../source.cpp`) into object files (`main_q2.o` and `../source.o`). Then, it links these object files to create the executable `main_q2`.

**2. make plot**

This command executes the compiled executable (`main_q2`) to generate the output CSV file containing numerical solutions. It then runs the Python script `plot_q2.py`, which reads the output CSV file and generates comparison plots for different time step sizes. This command requires to run **make** before use.

**3. make clean**

This command removes all compiled object files (`*.o`), the executable (`main_q2`), and the output CSV file (`output.csv`). It ensures that the project directory is cleaned up, removing any generated files that are not necessary for further development.

**4. make cls**

This command clears the terminal screen, while also deleting the compiled object files, executable, output files, and plots.
