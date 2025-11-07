import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann

# Constants
mass_of_argon = 39.948  # amu
dt = 0.5 
epsilon = 0.0103  # eV
sigma = 3.405  # A
num_steps = 1500
n = 5
d = 10.229*sigma /864**(1/3)
kb = 8.617333262e-5  # eV/K
T_in = 95
cutoff = 2.5 * sigma
cutoff2 = cutoff ** 2  

N = n**3
L = n * d


def initialize_positions(n, d):
    pos = np.zeros((N, 3))
    count = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                pos[count] = np.array([i * d, j * d, k * d])
                count += 1
    return pos


def initialize_velocities(T, N):
    """
    Intializes random velocitites, removes drift and rescales
    them according to desired temperature.
    """
    vel = np.random.randn(N, 3)  
    vel -= np.mean(vel, axis=0)  
    T_current = mass_of_argon * np.sum(vel**2) / (3 * kb * N)
    vel *= np.sqrt(T / T_current)  
    return vel


def get_acc_pot(X, L):
    """
    Computes the LJ-potential, forces and corresponding accelerations. 
    """
    forces = np.zeros((N, 3))
    potential_energy = np.zeros(N)
    D = X[None, :, :] - X[:, None, :]
    D -= np.rint(D / L) * L  
    D2 = np.sum(D**2, axis=-1)
    D2_safe = np.maximum(D2, 1e-10)  
    DS2 = np.divide(sigma**2, D2_safe)
    valid = D2 < cutoff2  
    DS6 = np.where(valid, DS2**3, 0) 
    
    for i in range(N - 1):
        for j in range(i + 1, N):
            r = D[i, j, :]
            r2 = D2[i, j]
            rs6 = DS6[i, j]
            F = np.where(valid[i, j], 48 * epsilon * rs6 * r * (rs6 - 0.5) / r2, 0)
            forces[i, :] -= F
            forces[j, :] += F
            potential_energy[i] += np.where(valid[i, j], 4 * epsilon * (rs6 * (rs6 - 1)), 0)
            potential_energy[j] += np.where(valid[i, j], 4 * epsilon * (rs6 * (rs6 - 1)), 0)
    
    return forces / mass_of_argon, potential_energy


def update_pos_vel(X, V, A, dt, L):
    """
    Updates positions and veocities using the velocity verlet method, 
    and ensures PBC using Minimum image convention .
    """
    X_new = X + V * dt + 0.5 * A * dt**2
    X_new -= np.rint( X_new / L ) * L
    A_new, potential_energy = get_acc_pot(X_new, L)
    V_new = V + 0.5 * (A + A_new) * dt
    return X_new, V_new, A_new, potential_energy


def K_energy(V):
    kinetic_energy = 0.5 * mass_of_argon * np.sum(V**2, axis=1)
    return kinetic_energy


def run_md(dt, num_steps, T):
    pos = initialize_positions(n, d)
    vel = initialize_velocities(T, N)
    acc, pot_energy = get_acc_pot(pos, L)
    pos_history = np.zeros((num_steps, N, 3))
    vel_history = np.zeros((num_steps, N, 3))
    ke_history = np.zeros((num_steps, N))
    pe_history = np.zeros((num_steps, N))
    
    for step in range(num_steps):
        pos, vel, acc, pot_energy = update_pos_vel(pos, vel, acc, dt, L)
        pos_history[step] = pos
        vel_history[step] = vel
        pe_history[step] = pot_energy
        ke_history[step] = K_energy(vel)
    
    return pos_history, vel_history, ke_history, pe_history


pos_hist, vel_hist, ke_hist, pe_hist = run_md(dt, num_steps, T_in)

def plot_energy(i):
    t_vals = np.linspace(0, dt * num_steps, num_steps)
    plt.figure(figsize=(10, 5))
    plt.plot(t_vals, ke_hist[:, i], label=f'Kinetic energy for atom {i+1}')
    plt.plot(t_vals, pe_hist[:, i], label=f'Potential energy for atom {i+1}')
    plt.plot(t_vals, ke_hist[:, i] + pe_hist[:, i], label=f'Total energy for atom {i+1}')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.title(f'Energy evolution for Atom {i+1}')
    plt.show()

def plot_x(i):
    t_vals = np.linspace(0, dt * num_steps, num_steps)
    plt.figure(figsize=(10, 5))
    #for i in range(pos_hist.shape[1] - 115):
    plt.plot(t_vals, pos_hist[:, i], linewidth=0.8)
    plt.xlabel('Time (ps)')
    plt.ylabel('Position (Å)')
    plt.legend()
    plt.title('x-position over time for atom {i+1}')
    plt.show()

def plot_v(i):
    t_vals = np.linspace(0, dt * num_steps, num_steps)
    plt.figure(figsize=(10, 5))
    # for i in range(pos_hist.shape[1] - 115):
    plt.plot(t_vals, vel_hist[:, i, 0], linewidth=0.8)
    plt.xlabel('Time (ps)')
    plt.ylabel('Velocity (Å/ps)')
    plt.legend()
    plt.title(f'x-Velocity over time for atom {i+1}')
    plt.show()


#plot_energy(2)
plot_x(0)
plot_v(0)

