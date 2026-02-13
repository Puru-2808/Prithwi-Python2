import numpy as np
import sympy as sp
from scipy import integrate
import matplotlib.pyplot as plt

xi = -1
xf = 1
n = 10001
h = (xf - xi) / (n - 1)

Ei = 0
Ef = 12
n1 = 1001
h1 = (Ef - Ei) / (n1 - 1)

alist = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
Energy_levels = []

# ---------- ENERGY (RK4 already correct) ----------
for a in alist:
    EE = []
    yy = []

    for i in range(n1):
        E = Ei + i * h1
        EE.append(E)
        
        y = 0.0
        z = 1.0
        
        for j in range(n-1):
            x = xi + j * h
            
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

            y += (k1y + 2*k2y + 2*k3y + k4y) / 6
            z += (k1z + 2*k2z + 2*k3z + k4z) / 6
            
        yy.append(y)

    for i in range(n1 - 1):
        if yy[i] * yy[i + 1] < 0:
            precise_E = EE[i] - yy[i] * (EE[i+1] - EE[i]) / (yy[i+1] - yy[i])
            Energy_levels.append(precise_E)
            break

degree = 2
coeffs = np.polyfit(alist, Energy_levels, degree)

E2 = coeffs[0]
E1 = coeffs[1]
E0_val = coeffs[2]

print(f"E0 (Unperturbed) = {E0_val:.8f}")
print(f"First order correction (E1): {E1:.8f}")
print(f"Second order correction (E2): {E2:.8f}")

b = 0.5
E_total = E0_val + b*E1 + (b**2)*E2
print("Total corrected ground state energy:", E_total)

# ---------- GET UNPERTURBED EIGENVALUES ----------
def get_roots(lambda_val):
    EE=[]
    yy=[]
    for i in range(n1):
        yi, zi = 0, 1
        E = Ei + i*h1
        EE.append(E)
        for j in range(n-1):
            x = xi + j*h
            def fz(xv, yv):
                return -2*(E - lambda_val*xv**4)*yv
            k1y = h*zi
            k1z = h*fz(x, yi)
            k2y = h*(zi + 0.5*k1z)
            k2z = h*fz(x + 0.5*h, yi + 0.5*k1y)
            k3y = h*(zi + 0.5*k2z)
            k3z = h*fz(x + 0.5*h, yi + 0.5*k2y)
            k4y = h*(zi + k3z)
            k4z = h*fz(x + h, yi + k3y)
            yi += (k1y + 2*k2y + 2*k3y + k4y)/6
            zi += (k1z + 2*k2z + 2*k3z + k4z)/6
        yy.append(yi)

    roots=[]
    for i in range(n1-1):
        if yy[i]*yy[i+1] < 0:
            E_root = EE[i] - yy[i]*(EE[i+1]-EE[i])/(yy[i+1]-yy[i])
            roots.append(E_root)
    return roots

EE1 = get_roots(0)
EE2 = get_roots(b)

num_states = min(len(EE1), len(EE2), 4)

# ---------- PLOT ----------
for i in range(num_states):

    # -------- UNPERTURBED --------
    E = EE1[i]
    yi, zi = 0, 1
    xx=[xi]
    yy=[yi]

    for j in range(n-1):
        x = xi + j*h

        def fz(xv, yv):
            return -2*E*yv

        k1y = h*zi
        k1z = h*fz(x, yi)
        k2y = h*(zi + 0.5*k1z)
        k2z = h*fz(x + 0.5*h, yi + 0.5*k1y)
        k3y = h*(zi + 0.5*k2z)
        k3z = h*fz(x + 0.5*h, yi + 0.5*k2y)
        k4y = h*(zi + k3z)
        k4z = h*fz(x + h, yi + k3y)

        yi += (k1y + 2*k2y + 2*k3y + k4y)/6
        zi += (k1z + 2*k2z + 2*k3z + k4z)/6

        xx.append(x)
        yy.append(yi)

    xx_arr = np.array(xx)
    yy_arr = np.array(yy)
    N = integrate.simpson(yy_arr**2, xx_arr)
    yy_arr /= np.sqrt(N)
    plt.plot(xx_arr, yy_arr, label=f"Unpert n={i}")

    # -------- PERTURBED --------
    E = EE2[i]
    yi, zi = 0, 1
    xx=[xi]
    yy=[yi]

    for j in range(n-1):
        x = xi + j*h

        def fz(xv, yv):
            return -2*(E - b*xv**4)*yv

        k1y = h*zi
        k1z = h*fz(x, yi)
        k2y = h*(zi + 0.5*k1z)
        k2z = h*fz(x + 0.5*h, yi + 0.5*k1y)
        k3y = h*(zi + 0.5*k2z)
        k3z = h*fz(x + 0.5*h, yi + 0.5*k2y)
        k4y = h*(zi + k3z)
        k4z = h*fz(x + h, yi + k3y)

        yi += (k1y + 2*k2y + 2*k3y + k4y)/6
        zi += (k1z + 2*k2z + 2*k3z + k4z)/6

        xx.append(x)
        yy.append(yi)

    xx_arr = np.array(xx)
    yy_arr = np.array(yy)
    N = integrate.simpson(yy_arr**2, xx_arr)
    yy_arr /= np.sqrt(N)
    plt.plot(xx_arr, yy_arr, '--', label=f"Pert n={i}")

plt.xlabel("x")
plt.ylabel("Wavefunction ψ(x)")
plt.title("Unperturbed vs Quartic Perturbed Eigenstates")
plt.legend()
plt.grid()
plt.show()