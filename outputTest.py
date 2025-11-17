import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.load("TwoBodyTest.npy", allow_pickle=True)

# Each entry in data is: [time, particle1, particle2, ...]
times = np.array([entry[0] for entry in data])
particle_names = [p.name for p in data[0][1:]]

# Build dictionary to store trajectories
trajectories = {name: {"x": [], "y": [], "z": [], "vx": [], "vy": [], "vz": [], "m": None}
                for name in particle_names}

# Extract data
for entry in data:
    t = entry[0]
    particles = entry[1:]

    for p in particles:
        name = p.name
        trajectories[name]["x"].append(p.position[0])
        trajectories[name]["y"].append(p.position[1])
        trajectories[name]["z"].append(p.position[2])
        trajectories[name]["vx"].append(p.velocity[0])
        trajectories[name]["vy"].append(p.velocity[1])
        trajectories[name]["vz"].append(p.velocity[2])
        trajectories[name]["m"] = p.mass
        

# -----------------------------
# 3D Plot of Orbits
# -----------------------------
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

for name in trajectories:
    x = trajectories[name]["x"]
    y = trajectories[name]["y"]
    z = trajectories[name]["z"]

    if name == "Sun":
        ax.scatter(x[0], y[0], z[0], label=name, s=60)
    else:
        ax.plot(x, y, z, label=name)

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_title("Orbital Paths from Simulation Output")
ax.legend()
plt.tight_layout()
plt.show()


# -----------------------------
# Energy Conservation Check
# -----------------------------
G = 6.67430e-11

def kinetic(m, vx, vy, vz):
    v2 = vx*vx + vy*vy + vz*vz
    return 0.5 * m * v2

def potential(m1, m2, r):
    return -G * m1 * m2 / r

energies = []

for frame in range(len(data)):
    particles = data[frame][1:]
    KE = 0
    PE = 0

    # Kinetic Energy
    for p in particles:
        KE += kinetic(p.mass, p.velocity[0], p.velocity[1], p.velocity[2])

    # Potential Energy (pairwise)
    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1, p2 = particles[i], particles[j]
            r = np.linalg.norm(p1.position - p2.position)
            PE += potential(p1.mass, p2.mass, r)

    energies.append(KE + PE)

energies = np.array(energies)

plt.figure(figsize=(8,5))
plt.plot(times, energ
