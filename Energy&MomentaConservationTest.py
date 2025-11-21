import numpy as np

# Load in the Test Data to Check for Accuracy
DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)
bodies0 = DataIn[0][1:]

G = 6.67408E-11

# E_k = 1/2 * m * v^2
# U = -GMm/r
E0_kinetic = 0
E0_potential = 0
for body in bodies0:
    E0_kinetic += 1/2 * body.mass * body.velocity**2

for i, body1 in enumerate(bodies0):
    for body2 in bodies0[i+1:]:
        r = np.linalg.norm(body2.position - body1.position)
        E0_potential -= G * body1.mass * body2.mass / r
        
E0_total = E0_kinetic + E0_potential
print(E0_total)

sampleStep = 1

bodiesFinal = DataIn[sampleStep][1:]

G = 6.67408E-11

# E_k = 1/2 * m * v^2
# U = -GMm/r
E_kinetic = 0
E_potential = 0
for body in bodiesFinal:
    E_kinetic += 1/2 * body.mass * body.velocity**2

for i, body1 in enumerate(bodiesFinal):
    for body2 in bodiesFinal[i+1:]:
        r = np.linalg.norm(body2.position - body1.position)
        E_potential -= G * body1.mass * body2.mass / r
        
E_total = E_kinetic + E_potential
print(E_total)

energyLost = E0_total - E_total
print(energyLost)
    
    
    
    
