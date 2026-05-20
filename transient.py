import numpy as np
import matplotlib.pyplot as plt

# INPUTS 
print("\n1D Transient Heat Conduction Solver\n")
L = float(input("Enter Length of slab (L): ") or 1.0)
nodes = int(input("Enter Number of nodes: ") or 6)
beta = float(input("Enter Fourier Number (beta): ") or 0.2)
T_initial_peak = float(input("Enter Initial Peak Temp at x=0: ") or 100.0)
bc_right_temp = float(input("Enter Right Boundary Temp: ") or 0.0)
epsilon = float(input("Enter Convergence Tolerance: ") or 1e-4)

# GRID INFO
dx = L / (nodes - 1)
print(f"\nGrid Info:")
print(f"Total Nodes: {nodes}")
print(f"Number of Divisions in slab: {nodes-1}")
print(f"Spatial Step (dx): {dx:.2f}")
print(" ")

# INITIALIZATION 
x = np.linspace(0, L, nodes)

# Initial Condition: T = 100 * cos(pi/2 * x)
T = T_initial_peak * np.cos((np.pi / 2) * x)
T_new = np.zeros(nodes)
T_history = [np.copy(T)]

print(f"Step n=0: " + ", ".join([f"{val:.2f}" for val in T]))
print("-" * 50)

# SOLVER LOOP
max_steps = 20000
converged = False

for n in range(1, max_steps + 1):
    # Calculate Internal Nodes (i = 1 to nodes-2)
    for i in range(1, nodes - 1):
        T_new[i] = beta * T[i+1] + (1 - 2*beta) * T[i] + beta * T[i-1]
    
    # BC Left (Insulated: dT/dx = 0) 
    # Using formula: T0_new = 2*beta*T1 + (1-2*beta)*T0
    T_new[0] = (2 * beta) * T[1] + (1 - 2 * beta) * T[0]
        
    # BC Right (Fixed Temperature)
    T_new[nodes-1] = bc_right_temp

    # Check Convergence
    diff = np.max(np.abs(T_new - T))
    
    # Print all steps
    #if n % 10 == 0 or n == 1: 
    print(f"Step n={n}: " + ", ".join([f"{val:.2f}" for val in T_new]))
    
    # Store history
    T_history.append(np.copy(T_new))

    if diff < epsilon:
        print(f"\nConverged at step {n}")
        converged = True
        break
    
    T = np.copy(T_new)

# FINAL OUTPUT 
if not converged:
    print("\nReached max steps without full convergence.")
    
    print("\nFinal Temperature Distribution:")
for i, temp in enumerate(T):
    print(f"Node {i} (x={x[i]:.2f}): {temp:.4f}")

# PLOTTING
T_history = np.array(T_history)
steps = np.arange(len(T_history))

num_cols = 2
num_rows = (nodes + 1) // 2
fig, axs = plt.subplots(num_rows, num_cols, figsize=(10, 3 * num_rows))
fig.suptitle('Transient Temperature for each Node', fontsize=14, fontname='Times New Roman', fontweight='bold')
axs = axs.flatten()

# Determine global min and max for consistent y-axis
y_min = np.min(T_history)
y_max = np.max(T_history)
y_margin = (y_max - y_min) * 0.05
if y_margin == 0: y_margin = 5

for i in range(nodes):
    axs[i].plot(steps, T_history[:, i], color='b', linestyle='-')
    axs[i].set_title(f"Node {i} (T{i}) vs Time Step", fontname='Times New Roman', fontweight='bold')
    axs[i].set_xlabel("Time Step (n)", fontname='Times New Roman')
    axs[i].set_ylabel(f"T{i} Temperature", fontname='Times New Roman')
    axs[i].set_ylim(y_min - y_margin, y_max + y_margin)
    axs[i].grid(True, linestyle='--', alpha=0.7)
for j in range(nodes, len(axs)):
    fig.delaxes(axs[j])

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.show()