import numpy as np
import sympy as sp

xi = -1
xf = 1
n = 10001
h = (xf - xi) / (n - 1)

Ei = 0
Ef = 3  
n1 = 2001 
h1 = (Ef - Ei) / (n1 - 1)

# lambda values (kept small for accurate derivative calculation at lam=0)
alist = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
Energy_levels = []

for a in alist:
    EE = []
    yy = []

    for i in range(n1):
        E = Ei + i * h1
        EE.append(E)
        
        y = 0.0
        z = 1.0 # Initial derivative y'(-1)
        
        # RK4 Integration Loop
        for j in range(n-1):
            x = xi + j * h
            
            # The system: 
            # dy/dx = z
            # dz/dx = -2 * (E - a * x**4) * y
            
            def f_z(x_val, y_val):
                return -2 * (E - a * x_val**4) * y_val

            k1y = h * z
            k1z = h * f_z(x, y)

            k2y = h * (z + 0.5 * k1z)
            k2z = h * f_z(x + 0.5 * h, y + 0.5 * k1y)

            k3y = h * (z + 0.5 * k2z)
            k3z = h * f_z(x + 0.5 * h, y + 0.5 * k2y)

            k4y = h * (z + k3z)
            k4z = h * f_z(x + h, y + k3y)

            y = y + (k1y + 2*k2y + 2*k3y + k4y) / 6
            z = z + (k1z + 2*k2z + 2*k3z + k4z) / 6
            
        yy.append(y)

    # Find the first energy level where y(xf) crosses zero
    for i in range(n1 - 1):
        if yy[i] * yy[i + 1] < 0:
            # Linear Interpolation for root refinement
            precise_E = EE[i] - yy[i] * (EE[i+1] - EE[i]) / (yy[i+1] - yy[i])
            Energy_levels.append(precise_E)
            break

# Symbolic Differentiation using Lagrange Interpolation
lam = sp.symbols('lam')
poly = 0
for i in range(len(alist)):
    term = Energy_levels[i]
    for j in range(len(alist)):
        if i != j:
            term *= (lam - alist[j]) / (alist[i] - alist[j])
    poly += term

degree = 3
coeffs = np.polyfit(alist, Energy_levels, degree)

E3 = coeffs[0]
E2 = coeffs[1]
E1 = coeffs[2]
E0_val = coeffs[3]

print(f"E0 (Unperturbed) = {E0_val:.8f}")
print(f"First order correction (E1): {E1:.8f}")
print(f"Second order correction (E2): {E2:.8f}")
print(f"Third order correction (E3): {E3:.8f}")

b = 0.5   # lambda value for which we want to calculate the total energy
E_total = E0_val + b*E1 + (b*2)*E2 + (b*3)*E3
print("Total corrected ground state energy:", E_total)