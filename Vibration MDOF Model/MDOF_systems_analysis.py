import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Constants as provided
m1 = 1.65
m2 = 0.7
m3 = 0.7
m4 = 0.7 
m5 = 410

k1 = 80000
k2 = 80000/(15-12)
k3 = k2   
k4 = k2   
r = 0.242
r2 = r**2 

M = np.array([
    [m1, 0., 0., 0., 0.],
    [0., m2, 0., 0., 0.],
    [0., 0., m3, 0., 0.],
    [0., 0., 0., m4, 0.],
    [0., 0., 0., 0., m5]
])

K = np.array([
    [k1*r2, -k1*r2, 0, 0, 0],
    [-k1*r2, r2*(k1+k2), -k2*r2, 0, -k2*r],
    [0, -k2*r2, (k2+k3)*r2, -k3*r2, r*(k2-k3)],
    [0, 0, -k3*r2, r2*(k3+k4), r*(k3-k4)],
    [0, -k2*r, r*(k2-k3), r*(k3-k4), k2+k3+k4]
])

# Initial conditions
t0 = 0.
x0 = np.zeros(10)  
x0[4] = 0.05  

def f(t):
    return np.zeros(5)  

def second_order_ode(t, y, M, K, f):
    n = len(y) // 2
    x = y[:n]  # Positions
    v = y[n:]  # Velocities
    dxdt = v
    dvdt = np.linalg.inv(M) @ (f(t) - np.dot(K, x))
    return np.concatenate([dxdt, dvdt])

# Solve the ODE
solution = solve_ivp(lambda t, y: second_order_ode(t, y, M, K, f), [t0, 10.], x0, t_eval=np.linspace(t0, 10., 1000))

# Access the solution
t = solution.t
x = solution.y[:5]  # Extract positions


unsrt_w, unsrt_modes = np.linalg.eig(np.linalg.inv(M) @ K)
sort_indices = np.argsort(unsrt_w)
sorted_w = unsrt_w[sort_indices]
sorted_modes = unsrt_modes[:, sort_indices]


MM = np.eye(5)
KK = np.zeros((5, 5))
PHI = np.zeros((5, 5))
for i in range(5):
    KK[i, i] = sorted_w[i]
    mode = sorted_modes[:, i].reshape(5, 1)
    PHI[:, i] = sorted_modes[:, i] / (mode.T @ M @ mode) ** 0.5


print("Modal Frequencies (Hz):\n\n", ["{:.4f}".format(freq) for freq in np.sqrt(np.abs(sorted_w))])
print("\nMode Shapes (normalized eigenvectors):\n\n", np.array([[f"{num:.4f}" for num in row] for row in PHI]))

# Plotting all displacements
plt.figure(figsize=(10, 10))
for i in range(5):  
    plt.subplot(5, 1, i + 1)  
    plt.plot(t, x[i], label=f'Displacement of Mass {i + 1}')
    plt.xlabel('Time (s)')
    plt.ylabel(f'Displacement (m)')
    plt.title(f'Displacement of Mass {i + 1} vs Time')
    plt.legend()
    plt.tight_layout()

plt.show()
