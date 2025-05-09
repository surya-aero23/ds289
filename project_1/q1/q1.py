from sympy import Matrix, symbols, pprint

#q1 a
# Defining symbolic variables
h1, h2, h3 = symbols('h1 h2 h3')

# Defining the matrix
_matrix = Matrix([[1, 1, 1, 1], [0, -h1, -h2, -h3],[0, h1**2/2, h2**2/2, h3**2/2],[0, -h1**3/6, -h2**3/6, -h3**3/6]])

# Print the matrix
print("Matrix and Inverse for q1 a:")
pprint(_matrix)
print("\n\n")

# Invert the matrix
_matrix_inv = _matrix.inv()

# Print the inverse matrix
print("Inverse:")
pprint(_matrix_inv)
print("\n\n")


# q1 b
h = symbols('h')
print("Inverse for q1 b:")
_matrix = Matrix([[1, 1, 1, 1, 1], [-h, -h, 0, h, h],[h, -h, 0, -h, h],[h**2, h**2, 0, h**2, h**2], [-h**2, h**2, 0, -h**2, h**2]])
_matrix_inv = _matrix.inv()
pprint(_matrix_inv)
