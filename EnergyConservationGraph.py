import numpy as np
import matplotlib.pyplot as plt

data_in = np.load("NBodyTest.npy", allow_pickle=True)
G = 6.67408E-11

e_total = []

for i, step in enumerate(data_in):
    bodies = step[1:]
    e_k = 0
    e_p = 0
    for body in bodies:
        e_k += 0.5 * body.mass * np.dot(body.velocity, body.velocity)
        
    for i, body_1 in enumerate(bodies):
        for body_2 in bodies[i+1:]:
            r = np.linalg.norm(body_2.position - body_1.position)
            e_p -= G * body_1.mass * body_2.mass / r
    e_total.append(e_k + e_p)

x = np.arange(0, len(data_in))

fig, ax = plt.subplots()
ax.plot(x, e_total)
plt.show()
    