import numpy as np
import matplotlib.pyplot as plt

data_in = np.load("NBodyTestWRogue.npy", allow_pickle=True)
G = 6.67408E-11

steps_per_plot = 5 # Allows you to plot less points with really large datasets

e_total = []
e_k_total = []
e_p_total = []

for i, step in enumerate(data_in):
    if i % steps_per_plot == 0:
        bodies = step[1:]
        e_k = 0
        e_p = 0
        for body in bodies:
            e_k += 0.5 * body.mass * np.dot(body.velocity, body.velocity) # Convert the Velocity to a Scalar and Calculate the Kinetic Energy.
            
        for i, body_1 in enumerate(bodies): # loop through the list of bodies 
            for body_2 in bodies[i+1:]: # Future Bodies are compared, this stops multiple comparisons between two objects
                r = np.linalg.norm(body_2.position - body_1.position) # Get the Scalar Distance
                e_p -= G * body_1.mass * body_2.mass / r # Calculate the Potential Energy
        # Append all Values to their Respective Lists
        e_total.append(e_k + e_p)
        e_k_total.append(e_k)
        e_p_total.append(e_p)

x = np.arange(0, len(data_in) // steps_per_plot) 

fig, axs = plt.subplots(3)

axs[0].plot(x, e_total)
axs[0].set_title("Total Energy in the System")
axs[1].plot(x, e_k_total)
axs[1].set_title("Total Kinetic Energy")
axs[2].plot(x, e_p_total)
axs[2].set_title("Total Potential Energy")
fig.tight_layout()
plt.show()
    