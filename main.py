# main calculator
# handles purely the Math Logic for Numerical Algorithms.
# returns the solutions and iteration loops back to whoever calls it (e.g. GUI)

import sympy as sp
import numpy as np

# =============================================================================#
# 1. ROOT FINDING ALGORI    THMS
# =============================================================================#

# %%
# bitsection method
def run_bisection(f, xl, xu, tol, max_iter):
    iterations_data = []
    
    # Initial Bracket check
    if f(xl) * f(xu) > 0:
        return None, 0, 0, "Failed", iterations_data
        
    xr, xr_old = 0, 0
    error = 100
    
    for i in range(max_iter):
        xr_old = xr
        xr = (xl + xu) / 2
        
        if i > 0 and xr != 0: 
            error = abs((xr - xr_old) / xr) * 100
            
        iterations_data.append({"Iter": i, "xl": xl, "xu": xu, "xr": xr, "Error %": error})
        
        if f(xl) * f(xr) > 0: 
            xl = xr
        else: 
            xu = xr
            
        if error <= tol and i > 0: 
            return xr, error, i+1, "Converged", iterations_data
            
    return xr, error, max_iter, "Max Iters Reached", iterations_data

# %%
# false position method
def run_false_position(f, xl, xu, tol, max_iter):
    iterations_data = []
    if f(xl) * f(xu) > 0: 
        return None, 0, 0, "Failed", iterations_data
        
    xr, xr_old = 0, 0
    error = 100
    
    for i in range(max_iter):
        xr_old = xr
        fxl, fxu = f(xl), f(xu)
        
        if fxl - fxu == 0: 
            break
            
        xr = xu - (fxu * (xl - xu)) / (fxl - fxu)
        
        if i > 0 and xr != 0: 
            error = abs((xr - xr_old) / xr) * 100
            
        iterations_data.append({"Iter": i, "xl": xl, "xu": xu, "xr": xr, "Error %": error})
        
        if f(xl) * f(xr) > 0: 
            xl = xr
        else: 
            xu = xr
            
        if error <= tol and i > 0: 
            return xr, error, i+1, "Converged", iterations_data
            
    return xr, error, max_iter, "Max Iters Reached", iterations_data

# %%
# fixed point method
def run_fixed_point(g, x0, tol, max_iter):
    iterations_data = []
    xi = x0
    error = 100
    
    for i in range(max_iter):
        xi_plus1 = g(xi)
        if xi_plus1 != 0: 
            error = abs((xi_plus1 - xi) / xi_plus1) * 100
            
        iterations_data.append({"Iter": i, "Xi": xi, "Xi+1": xi_plus1, "Error %": error})
        
        if error <= tol and i > 0: 
            return xi_plus1, error, i+1, "Converged", iterations_data
            
        xi = xi_plus1
        
    return xi, error, max_iter, "Max Iters Reached", iterations_data

# %%
# newton-raphson method
def run_newton(f, f_dash, x0, tol, max_iter):
    iterations_data = []
    xi = x0
    error = 100
    
    for i in range(max_iter):
        fx, fd = f(xi), f_dash(xi)
        if fd == 0: 
            return xi, error, i, "Zero Derivative", iterations_data
            
        xi_plus1 = xi - (fx / fd)
        
        if xi_plus1 != 0: 
            error = abs((xi_plus1 - xi) / xi_plus1) * 100
            
        iterations_data.append({"Iter": i, "Xi": xi, "F(xi)": fx, "Error %": error})
        
        if error <= tol and i > 0: 
            return xi_plus1, error, i+1, "Converged", iterations_data
            
        xi = xi_plus1
        
    return xi, error, max_iter, "Max Iters Reached", iterations_data

# %%
# secant method
def run_secant(f, x0, x1, tol, max_iter):
    iterations_data = []
    error = 100
    
    for i in range(max_iter):
        fx0, fx1 = f(x0), f(x1)
        if fx0 - fx1 == 0: 
            return x1, error, i, "Zero Denominator", iterations_data
            
        x2 = x1 - (fx1 * (x0 - x1)) / (fx0 - fx1)
        
        if x2 != 0: 
            error = abs((x2 - x1) / x2) * 100
            
        iterations_data.append({"Iter": i, "x0": x0, "x1": x1, "x2": x2, "Error %": error})
        
        if error <= tol and i > 0: 
            return x2, error, i+1, "Converged", iterations_data
            
        x0 = x1
        x1 = x2
        
    return x1, error, max_iter, "Max Iters Reached", iterations_data

# =============================================================================
# 2. LINEAR ALGEBRA SOLVERS
# =============================================================================

# %%
# gauss-jordan elimination method
def gauss_jordan_elimination(a):
    rows = 3
    A = [row[:] for row in a]
    for i in range(rows):
        pivot = A[i][i]
        if pivot == 0:
            raise ValueError("Zero pivot encountered in Gauss-Jordan!")
        for j in range(i, rows+1):
            A[i][j] /= pivot
        for k in range(rows):
            if k != i:
                factor = A[k][i]
                for j in range(i, rows+1):
                    A[k][j] -= factor * A[i][j]
    x = [A[i][rows] for i in range(rows)]
    return x

# %%
# lu decomposition method
def lu_decomposition(a):
    rows = 3
    A = [row[:] for row in a]
    b = [A[i][rows] for i in range(rows)]
    L = [[0.0] * rows for _ in range(rows)]
    U = [[0.0] * rows for _ in range(rows)]
    
    # Doolittle method
    for i in range(rows):
        for k in range(i, rows):
            total = sum(L[i][j] * U[j][k] for j in range(i))
            U[i][k] = A[i][k] - total
        for k in range(i, rows):
            if i == k:
                L[i][i] = 1
            else:
                if U[i][i] == 0:
                    raise ValueError("Zero found on U diagonal.")
                total = sum(L[k][j] * U[j][i] for j in range(i))
                L[k][i] = (A[k][i] - total) / U[i][i]
                
    # Lc = b
    c = [0.0] * rows
    for i in range(rows):
        total = sum(L[i][j] * c[j] for j in range(i))
        c[i] = (b[i] - total) / L[i][i]
        
    # Ux = c
    x = [0.0] * rows
    for i in range(rows-1, -1, -1):
        total = sum(U[i][j] * x[j] for j in range(i+1, rows))
        x[i] = (c[i] - total) / U[i][i]
        
    return x

# %%
def cramers_rule(a):
    rows = 3
    A = [[a[i][j] for j in range(rows)] for i in range(rows)]
    b = [a[i][rows] for i in range(rows)]
    
    def determinant(matrix, n):
        if n == 1: return matrix[0][0]
        if n == 2: return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        det = 0
        for c in range(n):
            submatrix = [[matrix[i][j] for j in range(n) if j != c] for i in range(1, n)]
            det += ((-1) ** c) * matrix[0][c] * determinant(submatrix, n - 1)
        return det

    detA = determinant(A, rows)
    if detA == 0:
        raise ValueError("Determinant is 0. Cramer's rule fails.")
        
    x = []
    for i in range(rows):
        Ai = [row[:] for row in A]
        for j in range(rows):
            Ai[j][i] = b[j]
        x.append(determinant(Ai, rows) / detA)
        
    return x

# %%
if __name__ == "__main__":
    import gui
    app = gui.NumericalApp()
    app.mainloop()